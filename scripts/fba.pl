#
# fba to be used with the psomcs program
#
# Author: Govind Nair
# email: govind.nair@boku.ac.at
#
#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use List::Util qw(max min sum);
use Pod::Usage;
use Getopt::Long;
use constant CONSIDERED_ZERO => 1.E-7;
use Math::CPLEX::OP;
use Math::CPLEX::Env;
use File::Basename;

my ($config_file, $parameters_file, $sfile, $rfile, $mfile, $rvfile, $help, $man, $targlim, $debug, $outfile, $optmz, $g_norm_recs);

GetOptions(	'configfile=s' => \$config_file,
		'parameters=s' => \$parameters_file,
		'sfile=s' => \$sfile,
		'rfile=s' => \$rfile,
		'meta=s' => \$mfile,
		'opti=s' => \$optmz,
                'rev=s' => \$rvfile)
    or exit;

my $norm_rec_rate = 10;
my ($cplex_env,$cplex_stat_opt,$lp);
my $num_reacs_expand;
my ($solutions,$sol_norm);
my $solution_cnt = 0;
my $write_lp_file = 1;
my $num_threads = 1;

my $dir2write = dirname $config_file;
$dir2write .= "/";

my $reacs = read_reactions($rfile);

my $metas = read_metabolites($mfile);
my $mmets = @$metas;
my $exmets = 0;
foreach my $met(@$metas){
    if($met =~ m:bb:){ $exmets++; }
}

my $revrer = read_reversibility($rvfile);

my $stoim = read_smatrix($sfile);
my $mets = scalar keys $stoim;

my $g_rxnMap = getRxnMappingFromRfile(); # R_xxx => i
my %g_rxnMap_rev = reverse %$g_rxnMap;

getParameters(); # only getting the norm rec indices

my ($medium, $constraints, $objective) = parseConfigFile($config_file);

# getting data on number of variables
my $numrecs = @$reacs;
my $numrevrecs = sum(@$revrer);
my $numirrevrecs = @$revrer-$numrevrecs;
my $numMets = @$metas;

# create and solve the LP
do_fba();

################################################################################
# modified from Matthias Gerstl's FBA code
################################################################################
sub parseConfigFile {
    my $file = shift;
    my ($md,$constr,$obj);
    my $confl;
    my $inMd = 0;
    my $inConstr = 0;
    my $inObj = 0;
    open($confl, "<", $file) or die "Coundn't open $file: $!";
    while(<$confl>){
    	if ($_ =~ /\D/){
	    chomp;
	    $_ =~ s/^\s+//;
	    $_ =~ s/\s+$//;
	    if ($_ !~ /^#/){
		if ($inMd){
		    if ($_ =~ /<\/medium>/){
			$inMd = 0;
		    } else {
			push(@{$md}, $_);
		    }
		} elsif ($inConstr) {
		    if ($_ =~ /<\/constraints>/){
			$inConstr = 0;
		    } else {
			if ($_ =~ /=/){
			    my $recid;
			    if ($_ =~ />=/){
				my @spl = split(/>=/, $_);
				$spl[0] =~ s/\s+//g;
				$spl[1] =~ s/\s+//g;
				$recid = $g_rxnMap->{$spl[0]};
				$constr->{$recid}{'relation'} = 'G';
				$constr->{$recid}{'value'} = $spl[1];
			    } elsif ($_ =~ /<=/){
				my @spl = split(/<=/, $_);
				$spl[0] =~ s/\s+//g;
				$spl[1] =~ s/\s+//g;
				$recid = $g_rxnMap->{$spl[0]};
				$constr->{$recid}{'relation'} = 'L';
				$constr->{$recid}{'value'} = $spl[1];
			    } else {
				my @spl = split(/=/, $_);
				$spl[0] =~ s/\s+//g;
				$spl[1] =~ s/\s+//g;
				$recid = $g_rxnMap->{$spl[0]};
				$constr->{$recid}{'relation'} = 'E';
				$constr->{$recid}{'value'} = $spl[1];
			    }
			} else {
			    die ("Config file is not correct: constraint has no upper or lower bound:\n$_\n");
			}
		    }
		} elsif ($inObj){
		    if ($_ =~ /<\/objective>/){
			$inObj = 0;
		    } else {
			my $positive = 1;
			my $name;
			if ($_ =~ /^-/){
			    $positive = -1;
			}
			$_ =~ s/^\+//;
			$_ =~ s/^\-//;
			$name = $_;
			my $objInd = $g_rxnMap->{$name};
			if ($objInd == -1){
			    die ("Could not find index of objective:\n$_\n");
			} else {
			    $obj->{$objInd} = 0;
			}
		    }
		} else {
		    if ($_ =~ /<medium>/){
			$inMd = 1;
		    } elsif ($_ =~ /<constraints>/){
			$inConstr = 1;
		    } elsif ($_ =~ /<objective>/){
			$inObj = 1;
		    }
		}
	    }
	}
    }
    close($confl);
    return ($md, $constr, $obj);
}
################################################################################

################################################################################
################################################################################
sub printArr{
    my $fl = shift;
    my $arr = shift;
    open(FL, ">$fl");
    foreach my $el(@{$arr}){
        print FL "$el\n";
    }
    close(FL);
}
################################################################################

################################################################################
################################################################################
sub do_fba
{
    # get CPLEX environment
    $cplex_env = Math::CPLEX::Env->new();
    die "ERROR: creating CPLEX environment failed!" unless $cplex_env;

    $cplex_stat_opt = Math::CPLEX::Base::CPX_STAT_OPTIMAL();

    # create LP problem
    $lp = $cplex_env->createOP();
    die "ERROR: couldn't create Linear Program\n" unless $lp;

    if((defined $optmz) && ($optmz eq 'min')){
	die "ERROR: minimize() failed\n" unless $lp->minimize();
    } else {
	die "ERROR: maximize() failed\n" unless $lp->maximize();
    }

    $cplex_env->setintparam(&Math::CPLEX::Base::CPX_PARAM_THREADS, $num_threads);
    
    fill_lp();

    #############################################################################
    # solve LP problem
    #############################################################################
    die "ERROR: lpopt() failed\n" unless $lp->lpopt();
    #############################################################################

    #############################################################################
    #############################################################################
    my $sol_name = "/tmp/lp_solution_fba.txt";
    die "ERROR: solwrite() failed\n" unless $lp->solwrite($sol_name);
    #############################################################################

    #############################################################################
    # retrieve computed values
    #############################################################################
    my $sol_recs;
    my ($sol_status, $obj_val, @vals) = $lp->solution();
    if (!defined($sol_status)){
	$obj_val = 0;
    }
    if ($sol_status != $cplex_stat_opt){
	die ("solution status is not optimal: $sol_status\n");
    }

    print "$obj_val\n";

    #############################################################################
    #############################################################################
    die "ERROR: free() failed\n" unless $lp->free();
    die "ERROR: close() failed\n" unless $cplex_env->close();
    #############################################################################
}
################################################################################

################################################################################
################################################################################
sub setLength
{
    my $inp = shift;
    my $len = shift;

    while( length $inp < $len )
    {
	$inp = '0' . $inp;
    }

    return $inp;
}
################################################################################

################################################################################
################################################################################
sub fill_lp
{
    #############################################################################
    # define columns
    #############################################################################
    my $obj_coefs;
    my $lower_bnd;
    my $upper_bnd;
    my $col_names;
    my $col_cnt = 0;

    # all variables
    for( my $i = 0; $i < @$reacs; $i++ )
    {
	(exists $objective->{$i}) ? ($obj_coefs->[$col_cnt] = 1.0):($obj_coefs->[$col_cnt] = 0.0);
	($revrer->[$i] == 1) ? ($lower_bnd->[$col_cnt] = - &Math::CPLEX::Base::CPX_INFBOUND):($lower_bnd->[$col_cnt] = 0.0);
	$upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
	if (exists $constraints->{$i}){
	    if($constraints->{$i}{relation} eq 'G'){
		$lower_bnd->[$col_cnt] = $constraints->{$i}{value};
	    } elsif ($constraints->{$i}{relation} eq 'L'){
		$upper_bnd->[$col_cnt] = $constraints->{$i}{value};
	    } else {
		$lower_bnd->[$col_cnt] = $constraints->{$i}{value};
		$upper_bnd->[$col_cnt] = $constraints->{$i}{value};
	    }
	}
	$col_names->[$col_cnt] = $reacs->[$i];
	$col_cnt++;
    }

    my $cols = { num_cols  => $col_cnt,
		 obj_coefs => $obj_coefs,
		 lower_bnd => $lower_bnd,
		 upper_bnd => $upper_bnd,
		 col_names => $col_names
    };
    die "ERROR: newcols() failed\n" unless $lp->newcols($cols);
    #############################################################################

    #############################################################################
    # add rows
    #############################################################################
    add_rows();
    #############################################################################

    #############################################################################
    # write lp file
    #############################################################################
    if( $write_lp_file )
    {
	my $filename = "$dir2write"."myCPLEX_fba.lp";
	die "ERROR: writeprob() failed\n" unless $lp->writeprob($filename);
    }
    #############################################################################
}
################################################################################

################################################################################
################################################################################
sub add_rows()
{
    my $newRows;
    my $rhs;
    my $sense;
    my $row_names;
    my $row_cnt = 0;

    # S*v = 0
    for( my $m = 0; $m < $mets; $m++ ){
	$rhs->[$row_cnt]                     = 0.0;
	$sense->[$row_cnt]                   = 'E';
	$row_names->[$row_cnt]               = "$metas->[$m]";

	for( my $r = 0; $r < @$reacs; $r++ )
	{
	    exists $stoim->{$m}{$r} ? ($newRows->[$row_cnt][$r] = $stoim->{$m}{$r}):($newRows->[$row_cnt][$r] = 0.0);
	}
	$row_cnt++;
    }

    # add sum norm reacs = 1
    $rhs->[$row_cnt]                     = 1.0; # change to 1 when uptake is positive
    $sense->[$row_cnt]                   = 'E';
    $row_names->[$row_cnt]               = "norm_recs";
    for( my $r = 0; $r < @$reacs; $r++ )
    {
	exists $g_norm_recs->{$r} ? ($newRows->[$row_cnt][$r] = 1.0):($newRows->[$row_cnt][$r] = 0.0);
    }
    $row_cnt++;

    my $rows = {num_rows  => $row_cnt,
		rhs       => $rhs,
		sense     => $sense,
		row_names => $row_names,
		row_coefs => $newRows};

    die "ERROR: addrows() failed\n" unless $lp->addrows($rows);
}
################################################################################

################################################################################
# from Jürgen Zanghellini
################################################################################
sub getRxnMappingFromRfile{
    open(RXNS,"<",$rfile);
    $_= <RXNS>;           # read reaction name
    s/^\s+//g;
    s/"//g;
    my @rxnIDs= split;
    close(RXNS);
    my $ii= 0;
    my %rxnMap = map { $_ => $ii++ } @rxnIDs;
    return \%rxnMap;
}
################################################################################

################################################################################
#read the stoichiometric matrix
################################################################################
sub read_smatrix{
    my $file = shift;
    my $sm;

    open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";

    while( <$fh> )
    {
	my $row_no = $. - 1;
	my @st_facs = split;
	for my $j(0 .. $#st_facs){
	    if($st_facs[$j] != 0){
		$sm->{$row_no}{$j} = $st_facs[$j];
	    }
	}
    }

    close $fh;
    return $sm;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub read_reversibility
{
    my $file = shift;
    my $rvs;

    open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
    $_ = <$fh>;
    @$rvs = split;
    close $fh;

    return $rvs;
}
################################################################################

################################################################################
# metabolite file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
# from Christian Jungreuthmayer
################################################################################
sub read_metabolites
{
    my $file = shift;
    my $mbs;

    open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
    $_ = <$fh>;
    s/"//g;
    s/#//g;
    @$mbs = split;
    close $fh;

    return $mbs;
}
################################################################################

################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
# from Christian Jungreuthmayer
################################################################################
sub read_reactions
{
    my $file = shift;
    my $rcs;

    open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
    $_ = <$fh>;
    s/"//g;
    s/>//g;
    s/#//g;
    @$rcs = split;
    close $fh;

    return $rcs;
}
################################################################################

################################################################################
# from stackoverflow
################################################################################
sub out{
    #to output results from a hash
    my ( $file, $hash_ref ) = @_;
    open my $fh, '>', $file or die "Can't write '$file': $!";
    local $Data::Dumper::Terse = 1;   # no '$VAR1 = '
    local $Data::Dumper::Useqq = 1;   # double quoted strings
    print $fh Dumper $hash_ref;
    close $fh or die "Can't close '$file': $!";
}
################################################################################

################################################################################
################################################################################
sub getParameters{
    # putting norm reactions into a new hash, $g_norm_recs
    open(my $paras, "<", $parameters_file);
    while(<$paras>){
        chomp;
        my @varvals = split(/=/, $_);
	if($varvals[0] =~ m|normalization_reactions|){
            my @norm_recs = split(' ',$varvals[1]);
	    foreach my $rec( @norm_recs ){
		my $idx = $g_rxnMap->{$rec};
		$g_norm_recs->{$idx} = 1;
	    }
        }
    }
    close($paras);
}
################################################################################



