#!/usr/bin/perl
#
# Author: Govind Nair
# email: govind.nair@boku.ac.at
#
#
# to calculate cMCS optimizing multiple reactions in a metabolic network
# using mixed integer linear programming (MILP) with particle swarm optimization.
#
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use List::Util qw(max min sum);
use Scalar::Util qw(looks_like_number);
use File::Basename;
use Math::Trig;
use constant CONSIDER_ZERO => 1e-8;
use Math::CPLEX::Env;
use Math::CPLEX::Base;
use threads;
use threads::shared;
use Thread::Queue;
use Getopt::Std;
use vars qw($opt_s $opt_m $opt_r $opt_v $opt_e $opt_o $opt_p $opt_h );

my ($g_parti_idx_fits, $g_population_size, $g_max_iterations, $g_norm_recs, $g_excluded_recs, $g_medium, $g_objective_reactions, $g_objective_maxflux, $g_cuts, $g_obj_best, $g_objective_minflux, $g_norm_constraints, $parameters_file, $g_population, $g_rxnMap, $g_num_reactions, $g_rxnMap_rev, $g_max_fit, $g_gadatfl, $g_gamcsfl, $g_particle_size, $g_model_name, $g_constraints, $g_decimal_places, $g_num_threads, $g_cmcs_limit, $g_ofile, $sfile, $rfile, $mfile, $rvfile, $efile, $g_stoim, $g_reacs, $g_rever, $g_metas);

read_arguments();

my $g_dir2write = dirname $parameters_file;
$g_dir2write .= "/";
$g_model_name = basename( $sfile );
$g_model_name =~ s/.\w+$//;
if( !defined $g_ofile) {
    $g_ofile = "$g_dir2write"."pso_output".".txt";
}
my $conffile = "$g_dir2write"."$g_model_name".".conf";
my $ga_datafile = "$g_dir2write"."$g_model_name"."pso_data".".dat";
my $ga_mcs_file = "$g_dir2write"."$g_model_name"."_cMCS".".txt";

open( my $psoout, '>', $g_ofile) or die "Couldn't open $g_ofile: $!";
print $psoout "Running psomcs ... \n";

@$g_medium = ();
$g_decimal_places = "%.6f";
$g_parti_idx_fits = {};
$g_cuts = {};
@$g_obj_best = ();
@$g_population = ();
share($g_population);
share($g_parti_idx_fits);
share($g_cuts);
share($g_max_fit);
share($g_obj_best);

getParameters();

# setting threads to (number of cores - 1), if not in parameters
if(!defined $g_num_threads) {
    $g_num_threads = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
    chomp($g_num_threads);
    $g_num_threads -= 1;
}
print $psoout "Using $g_num_threads threads\n";
print $psoout "... $g_population_size particles\n";
print $psoout "... $g_max_iterations iterations\n";
if(defined $g_cmcs_limit) {
    print $psoout "cMCS will be limited to size $g_cmcs_limit\n";
}

# reading the input files
$g_stoim = read_stoichiomat($sfile);
$g_reacs = read_reactions($rfile);
$g_metas = read_metabolites($mfile);
$g_rever = read_reversibility($rvfile);
$g_rxnMap = getRxnMappingFromRfile(); #rxn name => rxn number
$g_num_reactions = scalar keys %$g_rxnMap;
%$g_rxnMap_rev = reverse %$g_rxnMap; # rxn number => rxn name
if(-e $efile){
    $g_excluded_recs = read_reactions_2_hash($efile);
}
my $num_excluded_recs = scalar keys %$g_excluded_recs;

print $psoout "Number of reactions: $g_num_reactions\n";
print $psoout "Number of reactions excluded from knockouts: $num_excluded_recs\n";

for(my $j=0; $j<@$g_norm_constraints; $j++){
    $g_constraints->{$g_norm_recs->[$j]}{'relationship'} = 'L'; # change to 'L/G' when uptake is positive/negative
    $g_constraints->{$g_norm_recs->[$j]}{'val'} = $g_norm_constraints->[$j];
}

# getting the maximum values for objective/s stored in $g_objective_maxflux
get_obj_maxmin_fluxes();

for(my $i=0; $i<@$g_objective_reactions; $i++){
    print $psoout "$g_objective_reactions->[$i] maximum yield: $g_objective_maxflux->[$i]\n";
    print $psoout "$g_objective_reactions->[$i] minimum yield: $g_objective_minflux->[$i]\n";
}

# making the particles
generate_initial_population();

# opening the pso data file 
open($g_gadatfl, '>', $ga_datafile) or die "Couldn't open $ga_datafile: $!";
printf $g_gadatfl ("%10s %10s", "#Iteration", "  MaxFit   ");
foreach my $rec(@$g_objective_reactions){
    my $rec_pri = substr($rec, 0 , 12);
    printf $g_gadatfl ("%12s ", "$rec_pri");
}
printf $g_gadatfl ("%9s %8s", "CutSize", "Time");
print $g_gadatfl "\n";

# the PSO runs
$g_max_fit = 0;
run_pso($psoout);
print $psoout "psomcs has finished running\n";
close($g_gadatfl);

# writing the cutsets to file
my $num_total_mcs = scalar keys %$g_cuts;
print $psoout "A total of $num_total_mcs cMCS were found\n";
print $psoout "Writing cMCS to $ga_mcs_file in descending order of fitness\n";
write_cutsets2file($ga_mcs_file);

close ($psoout);


################################################################################
################################################################################
sub cal_fitness {
    my ($cutsize, $obj_vals) = @_;

    my $fitness = 1;
    if($cutsize) {
	for(my $i = 0; $i<@$g_objective_reactions; $i++){
	    $fitness *= $obj_vals->[$i]/$g_objective_maxflux->[$i];
	}
	$fitness *= (1 - ($cutsize/$g_num_reactions));
    } else {
	$fitness = 0;
    }
    return $fitness;
}
################################################################################

################################################################################
################################################################################
sub write_cutsets2file {
    my $cuts_file = shift;
    my $cuts_fits;

    foreach my $cut_string_bin(keys %$g_cuts) {
    my $cut_size = unpack("%64b*", $cut_string_bin);
    my @cut_vals = @{$g_cuts->{$cut_string_bin}};
    my $fitness = cal_fitness($cut_size, \@cut_vals);
    $cuts_fits->{$cut_string_bin} = $fitness;
    }

    open($g_gamcsfl, '>', $cuts_file) or die "Couldn't open $cuts_file: $!";
    foreach my $cut_string_bin(sort {$cuts_fits->{$b} <=> $cuts_fits->{$a} } keys %$g_cuts) {
	my $cut_string_unmodified = unModifyString($cut_string_bin,$g_num_reactions);
	my @cutvals = split('', $cut_string_unmodified);
	for(my $i = 0; $i < @cutvals; $i++){
	    if($cutvals[$i]) {
		print $g_gamcsfl "$g_rxnMap_rev->{$i} ";
	    }
	}
	print $g_gamcsfl "\n";
    }
    close($g_gamcsfl);
}
################################################################################

################################################################################
# doing a final round of fba on the cuts for getting the exact max/min vals
################################################################################
sub do_final_fba {

    foreach my $cut_string_bin(keys %$g_cuts) {
	my $obj_minimum_vals;
	my $obj_maximum_vals;
	my $constraints;

	my $cut_string_unmodified = unModifyString($cut_string_bin,$g_num_reactions);
	my @cutvals = split('', $cut_string_unmodified);
	my $cutsize = 0;

	foreach my $rec(keys %$g_constraints){
	    $constraints->{$rec}{'relationship'} = $g_constraints->{$rec}{'relationship'};
	    $constraints->{$rec}{'val'} = $g_constraints->{$rec}{'val'};
	}

	my @cutsname;
	for(my $i = 0; $i < @cutvals; $i++) {
	    if($cutvals[$i]){
		$cutsize++;
		push @cutsname, $g_rxnMap_rev->{$i};
		my $rec = $g_rxnMap_rev->{$i};
		$constraints->{$rec}{'relationship'} = 'E';
		$constraints->{$rec}{'val'} = 0;
	    }
	}

	for(my $i=0; $i<@$g_objective_reactions; $i++) {
	    my $conffile = write_conffile($g_model_name,$g_dir2write,$g_medium,$constraints,$g_objective_reactions->[$i]);
	    $obj_minimum_vals->[$i] = do_fba_maxmin($conffile,1);
	    $obj_maximum_vals->[$i] = do_fba_maxmin($conffile);
	}

	# new fitness for the final cuts
	my $fitness = 0;
	print "@cutsname : ";
	print "$obj_minimum_vals->[0] ";
	$fitness = $obj_minimum_vals->[0];
	for(my $i = 1; $i<@$g_objective_reactions; $i++){
	    $fitness *= $obj_minimum_vals->[$i];
	    print "$obj_minimum_vals->[$i] ";
	}

	print ": $fitness\n";
    }
}
################################################################################

################################################################################
################################################################################
sub print_gen_data{
    my ($num_gen, $gen_time) = @_;

    printf  $g_gadatfl ("%10d %10.5f ", $num_gen, $g_max_fit);
    for(my $i = 0; $i < @$g_objective_reactions; $i++){
	printf  $g_gadatfl ("%12f ", $g_obj_best->[$i]);
    }
    printf  $g_gadatfl ("%9d ", $g_obj_best->[-1]);
    printf  $g_gadatfl ("%8d ", $gen_time);
    printf  $g_gadatfl "\n";
}
################################################################################

################################################################################
################################################################################
sub standard_deviation {
    my ($values, $mean) = @_;
    my @values_var = map { ($_ - $mean)**2 } @$values;
    my $variance = sum(@values_var)/(scalar @values_var);
    my $std_dev = sqrt $variance;

    return $std_dev;
}

################################################################################
################################################################################
sub median {
    my $values = shift;
    my $median_fitness;
    my @sorted_fitnesses = sort { $a <=> $b } @$values;
    my $pop_modulo = $#sorted_fitnesses % 2;
    if($pop_modulo == 0){
        $median_fitness = $sorted_fitnesses[($#sorted_fitnesses)/2];
    } else {
        $median_fitness = ($sorted_fitnesses[($#sorted_fitnesses-1)/2] + $sorted_fitnesses[($#sorted_fitnesses+1)/2])/2;
    }
    return $median_fitness;
}
################################################################################

################################################################################
################################################################################
sub generate_initial_population{
    $g_particle_size = @$g_objective_reactions; # setting chromosome length
    for(my $i = 1; $i<= $g_population_size; $i++){
	$g_population->[$i-1] = &share( [] );
	my $particle;
	@$particle = ();
	share($particle);
	$particle = make_particle($g_particle_size);
	@{$g_population->[$i-1]} = @$particle;
    }
}
################################################################################

################################################################################
################################################################################
sub dual_worker{
    my $tid = threads->tid();
    my $parti_idx = shift;
    do_dual($parti_idx);
}
################################################################################

################################################################################
################################################################################
sub run_pso{
    my $out_fl = shift;
    print $out_fl "Processing iteration: ";

    #looping through the genartions
    for(my $num_iter = 1; $num_iter <= $g_max_iterations; $num_iter++){
	print $out_fl "...$num_iter ";
	my $t_bfgen = time;
	my $partiqueue;

	$partiqueue = Thread::Queue->new();

	for (1..$g_num_threads) {
	    async {
		while (defined(my $parti_idx = $partiqueue->dequeue())) {
		    dual_worker($parti_idx);
		}
	    };
	}

	for(my $i = 0; $i < @$g_population; $i++){
	    $partiqueue->enqueue($i);
	}

	$partiqueue->end();

	$_->join() for threads->list();

	# update particles
	update_particles5() unless ($num_iter == $g_max_iterations);

	my $t_afgen = time-$t_bfgen;

	# print generation results to file
	print_gen_data($num_iter,$t_afgen);

    }
    print $out_fl "\n";
}
################################################################################

################################################################################
# with ring topology
################################################################################
sub update_particles3 {
    my $phi1 = 2;
    my $phi2 = 2;
    my $chi = 0.7298; # ClercnKennedy 2oo2
    my $thereshold = 0.997;

    # calculate 
    for(my $i = 0; $i < @$g_population; $i++){
	my $chromosome = $g_population->[$i];
	my $idxs_fits;
	my $max_fit_idx = 0;
	my $max_fit = 0;

	# getting the neighbour particles, 2 now
	# idxs before current
	my $t_idx = $i-1;
	if($i == 0){
	    $t_idx = @$g_population-1;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# idxs after current
	$t_idx = $i+1;
	if($i == $#$g_population){
	    $t_idx = 0;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}  else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# also including the current particle
	$idxs_fits->{$i} = $g_parti_idx_fits->{$i};

	# getting the max fitness and corresponding idx
	foreach my $idx(keys %$idxs_fits){
	    if($idxs_fits->{$idx} > $max_fit){
		$max_fit = $idxs_fits->{$idx};
		$max_fit_idx = $idx;
	    }
	}

	my $n_fit_chromosome = $g_population->[$max_fit_idx];

	for(my $j = 0; $j < @$g_objective_reactions; $j++){
	    my $rand1 = rand();
	    my $rand2 = rand();
	    if($rand1 > $thereshold) {
		$rand1 = 1.0;
	    }
	    if($rand2 > $thereshold) {
		$rand2 = 1.0;
	    }
	    my $old_obj_pos = $j + @$g_objective_reactions + 1;
	    my $old_vel_pos = $j + 2*@$g_objective_reactions + 2;
	    my $nu_vel = $phi1*$rand1*($n_fit_chromosome->[$old_obj_pos] - $chromosome->[$j]) + $phi2*$rand2*($g_obj_best->[$j] - $chromosome->[$j]);
	    $nu_vel *= $chi;
	    $nu_vel = sprintf $g_decimal_places, $nu_vel;
	    $chromosome->[$old_vel_pos] = $nu_vel;
	    $chromosome->[$j] += $nu_vel;

	    # checking if value greater than maximum or less than minimum
	    # if greater, then assigning random value
	    if(($chromosome->[$j] > $g_objective_maxflux->[$j]) || ($chromosome->[$j] < $g_objective_minflux->[$j])){
		my $min_val = 0.0001;
		my $diff = $g_objective_maxflux->[$j] - $min_val;
		my $val = $min_val + rand($diff);
		$chromosome->[$j] = $val;
	    }
	    $chromosome->[$j] = sprintf $g_decimal_places, $chromosome->[$j];
	}
    }
}
################################################################################

################################################################################
# with von Neumann topology
################################################################################
sub update_particles4 {
    my $phi1 = 2;
    my $phi2 = 2;
    my $chi = 0.7298; # ClercnKennedy 2oo2
    my $thereshold = 0.997;

    # calculate 
    for(my $i = 0; $i < @$g_population; $i++){
	my $chromosome = $g_population->[$i];
	my $idxs_fits;
	my $max_fit_idx = 0;
	my $max_fit = 0;

	# getting the neighbour particles, 4 = von Neumann
	# idxs before current 
	my $t_idx = $i-1;
	if($i == 0){
	    $t_idx = @$g_population-1;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx--;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} elsif ($i == 1){
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx = @$g_population-1;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx--;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# idxs after current
	$t_idx = $i+1;
	if($i == $#$g_population){
	    $t_idx = 0;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx++;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} elsif ($i == ($#$g_population-1)){
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx = 0;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx++;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# also including the current particle
	$idxs_fits->{$i} = $g_parti_idx_fits->{$i};

	# getting the max fitness and corresponding idx
	foreach my $idx(keys %$idxs_fits){
	    if($idxs_fits->{$idx} > $max_fit){
		$max_fit = $idxs_fits->{$idx};
		$max_fit_idx = $idx;
	    }
	}

	my $n_fit_chromosome = $g_population->[$max_fit_idx];

	for(my $j = 0; $j < @$g_objective_reactions; $j++){
	    my $rand1 = rand();
	    my $rand2 = rand();
	    if($rand1 > $thereshold) {
		$rand1 = 1.0;
	    }
	    if($rand2 > $thereshold) {
		$rand2 = 1.0;
	    }
	    my $old_obj_pos = $j + @$g_objective_reactions + 1;
	    my $old_vel_pos = $j + 2*@$g_objective_reactions + 2;
	    my $nu_vel = $phi1*$rand1*($n_fit_chromosome->[$old_obj_pos] - $chromosome->[$j]) + $phi2*$rand2*($g_obj_best->[$j] - $chromosome->[$j]);
	    $nu_vel *= $chi;
	    $nu_vel = sprintf $g_decimal_places, $nu_vel;
	    $chromosome->[$old_vel_pos] = $nu_vel;
	    $chromosome->[$j] += $nu_vel;
	    # checking if value greater than maximum or less than minimum
	    # if greater, then assigning random value
	    if(($chromosome->[$j] > $g_objective_maxflux->[$j]) || ($chromosome->[$j] < $g_objective_minflux->[$j])){
		my $min_val = 0.0001;
		my $diff = $g_objective_maxflux->[$j] - $min_val;
		my $val = $min_val + rand($diff);
		$chromosome->[$j] = $val;
	    }

	    $chromosome->[$j] = sprintf $g_decimal_places, $chromosome->[$j];
	}
    }
}
################################################################################

################################################################################
# Gong and Zhang Small World topology
################################################################################
sub update_particles5 {
    my $phi1 = 2;
    my $phi2 = 2;
    my $chi = 0.7298; # ClercnKennedy 2oo2
    my $thereshold = 0.997;

    # calculate 
    for(my $i = 0; $i < @$g_population; $i++){
	my $chromosome = $g_population->[$i];
	my $idxs_fits;
	my $max_fit_idx = 0;
	my $max_fit = 0;

	# getting random idxs
	my $random_idx = int rand (scalar @$g_population);
	$idxs_fits->{$random_idx} = $g_parti_idx_fits->{$random_idx};

	# idxs before current 
	my $t_idx = $i-1;
	if($i == 0){
	    $t_idx = @$g_population-1;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx--;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} elsif ($i == 1){
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx = @$g_population-1;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx--;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# idxs after current
	$t_idx = $i+1;
	if($i == $#$g_population){
	    $t_idx = 0;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx++;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} elsif ($i == ($#$g_population-1)){
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx = 0;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	} else {
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	    $t_idx++;
	    $idxs_fits->{$t_idx} = $g_parti_idx_fits->{$t_idx};
	}

	# getting the max fitness and corresponding idx
	foreach my $idx(keys %$idxs_fits){
	    if($idxs_fits->{$idx} > $max_fit){
		$max_fit = $idxs_fits->{$idx};
		$max_fit_idx = $idx;
	    }
	}

	my $n_fit_chromosome = $g_population->[$max_fit_idx];

	for(my $j = 0; $j < @$g_objective_reactions; $j++){
	    my $rand1 = rand();
	    my $rand2 = rand();
	    if($rand1 > $thereshold) {
		$rand1 = 1.0;
	    }
	    if($rand2 > $thereshold) {
		$rand2 = 1.0;
	    }
	    my $old_obj_pos = $j + @$g_objective_reactions + 1;
	    my $old_vel_pos = $j + 2*@$g_objective_reactions + 2;
	    my $nu_vel = $phi1*$rand1*($n_fit_chromosome->[$old_obj_pos] - $chromosome->[$j]) + $phi2*
$rand2*($g_obj_best->[$j] - $chromosome->[$j]);
	    $nu_vel *= $chi;
	    $nu_vel = sprintf $g_decimal_places, $nu_vel;
	    $chromosome->[$old_vel_pos] = $nu_vel;
	    $chromosome->[$j] += $nu_vel;
	    # checking if value greater than maximum or less than minimum
	    # if greater, then assigning random value
	    if(($chromosome->[$j] > $g_objective_maxflux->[$j]) || ($chromosome->[$j] < $g_objective_minflux->[$j])){
		my $min_val = 0.0001;
		my $diff = $g_objective_maxflux->[$j] - $min_val;
		my $val = $min_val + rand($diff);
		$chromosome->[$j] = $val;
	    }
	    $chromosome->[$j] = sprintf $g_decimal_places, $chromosome->[$j];
	}
    }
}
################################################################################

################################################################################
################################################################################
sub do_dual{
    my $idx = shift;
    my $chromosome = $g_population->[$idx];
    my $new_cut_sz_pos = @$g_objective_reactions; # pos in chrom for cut size
    my $old_cut_sz_pos = 2*$new_cut_sz_pos + 1; # pos corresponding to cut sz of previous best

    my $cutsets;
    $cutsets = run_defig($chromosome);

    if (ref($cutsets) eq "ARRAY") {
	lock($g_cuts);
	foreach my $mcs(@$cutsets) {
	    my $cut_pos;
	    @$cut_pos = map {$g_rxnMap->{$_}} @$mcs;

	    my $cut_string_bin = generateString0($g_num_reactions,$cut_pos);
	    if(exists $g_cuts->{$cut_string_bin}){
		for(my $i = 0; $i < @$g_objective_reactions; $i++){
		    $g_cuts->{$cut_string_bin}[$i] = $chromosome->[$i];
		}
	    } else {
		$g_cuts->{$cut_string_bin} =  &share( [] );
		for(my $i = 0; $i < @$g_objective_reactions; $i++){
		    $g_cuts->{$cut_string_bin}[$i] = $chromosome->[$i];
		}
	    }
	    # inserting new cut size into chromosome
	    $chromosome->[$new_cut_sz_pos] = scalar @$mcs;

	}
    }  else {
	# set the new cut size to 0
	$chromosome->[$new_cut_sz_pos] = 0;
    }

    # calculate fitness
    my $new_fitness = 0; # based on new positions
    my $pbest_fitness = 0; # for previous best
    my $cutsize = 0;
    my $obj_vals;

    # calculating the new fitness
    $cutsize = $chromosome->[$new_cut_sz_pos];
    for(my $i = 0; $i < @$g_objective_reactions; $i++){
	$obj_vals->[$i] = $chromosome->[$i];
    }
    $new_fitness = cal_fitness($cutsize, $obj_vals);

    # calculating the old fitness
    $cutsize = $chromosome->[$old_cut_sz_pos];
    for(my $i = 0; $i < @$g_objective_reactions; $i++){
	my $old_obj_pos = $i + @$g_objective_reactions + 1;
	$obj_vals->[$i] = $chromosome->[$old_obj_pos];
    }
    $pbest_fitness = cal_fitness($cutsize, $obj_vals);

    # update particle best
    if($new_fitness > $pbest_fitness){
	for(my $i = 0; $i < @$g_objective_reactions; $i++){
	    my $old_obj_pos = $i + @$g_objective_reactions + 1;
	    $chromosome->[$old_obj_pos] = $chromosome->[$i];
	}
	$chromosome->[$old_cut_sz_pos] = $chromosome->[$new_cut_sz_pos];
    }

    # update global best
    if($new_fitness > $g_max_fit){
	lock($g_max_fit);
	lock($g_obj_best);
	for(my $i = 0; $i <= @$g_objective_reactions; $i++){
	    $g_obj_best->[$i] = $chromosome->[$i]; # obj vals + cut size
	}
	$g_obj_best->[$old_cut_sz_pos] = $chromosome->[$new_cut_sz_pos];
	$g_max_fit = $new_fitness;
    }

    # updating $g_parti_idx_fits, should be equal to the previous best
    # this info is used in updating the particles
    {
	lock($g_parti_idx_fits);
	$g_parti_idx_fits->{$idx} = $pbest_fitness;
    }

}
################################################################################

################################################################################
################################################################################
sub read_mcs_file {
    my $mcsfile = shift;
    my $cutsets = ();

    open(my $mcsfl, '<', $mcsfile) or die "Couldn't open $mcsfile: $!";
    while(<$mcsfl>){
	if($_){
	    chomp;
	    s/^\s+//g;
	    s/"//g;
	    s/\s+$//g;
	    my @cut= split;

	    push @$cutsets, \@cut;
	}
    }
    return $cutsets;
}
################################################################################

################################################################################
################################################################################
sub getParameters{
    open(my $paras, "<", $parameters_file);
    while(<$paras>){
        chomp;
        my @varvals = split(/=/, $_);
	if($varvals[0] =~ m|objective_reactions|){
            @$g_objective_reactions = split(' ',$varvals[1]);
        }elsif($varvals[0] =~ m|normalization_reaction|){
            @$g_norm_recs = split(' ',$varvals[1]);
        }elsif($varvals[0] =~ m|medium|){
            @$g_medium = split(' ',$varvals[1]);
        }elsif($varvals[0] =~ m|normalization_constraints|){
            @$g_norm_constraints = split(' ',$varvals[1]);
        }elsif($varvals[0] =~ m|population_size|){
            $g_population_size = $varvals[1];
        }elsif($varvals[0] =~ m|max_knockout_size|){
            $g_cmcs_limit = $varvals[1];
        }elsif($varvals[0] =~ m|threads|){
            $g_num_threads = $varvals[1];
        }elsif($varvals[0] =~ m|max_iterations|){
            $g_max_iterations = $varvals[1];
        }
    }
    close($paras);
}
################################################################################

################################################################################
################################################################################
sub unModifyString{
    my $str = shift;
    my $str_length = shift;
    my $org_str = unpack('b*', $str);
    $org_str = substr $org_str,0,$str_length;
    return $org_str;
}
################################################################################

################################################################################
################################################################################
sub getModifiedStr{
    my @string_splits = unpack('(a64)*', shift);
    my $length = 64-(length $string_splits[-1]);
    if($length){
        for(my $i=0;$i<$length;$i++){
	    $string_splits[-1] .= '0';
	}
    }
    my $return_string = join('', @string_splits);
    my $return_string_bin = pack('b*', $return_string);
    return $return_string_bin;
}
################################################################################

################################################################################
# generates a string of 0s of given length with 1s at the given positions
################################################################################
sub generateString0{
    my $length = shift;
    my $positions_list = shift;
    my @string = (0)x$length;
    if(defined $positions_list){
	for(my $i=0;$i<@$positions_list;$i++){
	    my $pos = $positions_list->[$i];
	    $string[$pos] = 1;
	}
    }
    my $string0 = join('', @string);
    $string0 = getModifiedStr($string0);
    return $string0;
}
################################################################################

################################################################################
################################################################################
sub get_obj_maxmin_fluxes{
    for(my $i=0; $i<@$g_objective_reactions; $i++){
	# write config file
	my $conffile = write_conffile($g_model_name,$g_dir2write,$g_medium,$g_constraints,$g_objective_reactions->[$i]);

	# calculate max flux using fba
	$g_objective_maxflux->[$i] = do_fba_maxmin($conffile);
	$g_objective_minflux->[$i] = do_fba_maxmin($conffile,1);

    }
}
################################################################################

################################################################################
################################################################################
sub do_fba_maxmin{
    my ($conffile,$min) = @_;
    my $fba_cmd = '';
    # calculate flux using fba
    my $flux = 0;
    if(defined($min) && ($min == 1)) {
	$fba_cmd = "perl fba.pl --parameters $parameters_file --configfile $conffile --sfile $sfile --rfile $rfile --meta $mfile --rev $rvfile --opti min";
    } else {
	$fba_cmd = "perl fba.pl --parameters $parameters_file --configfile $conffile --sfile $sfile --rfile $rfile --meta $mfile --rev $rvfile";
    }
    my @fba_results = `$fba_cmd`;
    $flux = $fba_results[0];
    chomp($flux);
    if(!looks_like_number($flux)){
	die "Error! Unable to optimize flux: $!";
    }	
    return $flux;
}
################################################################################

################################################################################
################################################################################
sub write_conffile {
    my ($name,$dirname,$medium,$constraints,$objective) = @_;
    # $constraints - rec->{relationship},rec->{val}

    my $conffile = "$dirname"."$name".".conf";

	# write config file
	open(my $confl, ">", $conffile) or die "Couldn't open $conffile: $!";
	print $confl "<medium>\n";
	if(@$g_medium){
	    foreach my $met(@$medium){
		print $confl "$met\n";
	    }
	}
	print $confl "</medium>\n\n";
	print $confl "<constraints>\n";
	foreach my $rec(keys %$constraints){
	    if($constraints->{$rec}{relationship} eq 'E'){
		print $confl "$rec = $constraints->{$rec}{val}\n";
	    } elsif ($constraints->{$rec}{relationship} eq 'L'){
		print $confl "$rec <= $constraints->{$rec}{val}\n";
	    } elsif ($constraints->{$rec}{relationship} eq 'G'){
		print $confl "$rec >= $constraints->{$rec}{val}\n";
	    }
	}
	print $confl "</constraints>\n\n";
	print $confl "<objective>\n";
	print $confl "$objective\n";
	print $confl "</objective>";
	close($confl);

    return $conffile;
}
################################################################################

################################################################################
################################################################################
sub make_particle{
    my $chrom_length = shift;
    if($chrom_length <= 0){
	print "Error: Particle should have at least 1 element!!\n";
	exit 1;
    }
    my $chromosome;

    for(my $i=0; $i<$chrom_length; $i++){
	my ($obj_max_val,$val);
	if($i == 0) {
	    my $min_val = 0.0001;
	    my $diff = $g_objective_maxflux->[$i] - $min_val;
	    $val = $min_val + rand($diff);
	} else {
	    my $min_val = 0.0001;
	    my $diff = $g_objective_maxflux->[$i] - $min_val;
	    $val = $min_val + rand($diff);
	}

	$val = sprintf $g_decimal_places, $val;
	my $realval = 0;
	if($val == 0){
	    print "Problem!! 0 value generated\n";
	}

	# convert $val to a value in the variable range
	$realval = $val;
	$chromosome->[$i] = $realval; # position for storing current and init val
	$chromosome->[$i + $chrom_length + 1] = 0; # position for storing previous best
    }
    $chromosome->[$chrom_length] = 0; # position for storing the current mcs size
    $chromosome->[2*$chrom_length + 1] = 0; # position for storing the previous best mcs size

    # position for velocities
    for(my $i=0; $i<$chrom_length; $i++){
	my $v_max = $g_objective_maxflux->[$i]/2 + 0.001;
	my $v = rand($v_max);
	$v = sprintf $g_decimal_places, $v;
	push @$chromosome, $v; 
    }
    return $chromosome;
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
    #s/_//g;
    my @rxnIDs= split;
    close(RXNS);
    my $ii= 0;
    my %rxnMap = map { $_ => $ii++ } @rxnIDs;
    return \%rxnMap;
}
################################################################################

################################################################################
################################################################################
sub run_defig {
    my $chromosome = shift;
    my($dual_stoim, $dual_metas, $dual_reacs, $dual_rever);

    my($desired, $desired_stoich, $desired_inhomo) = make_DT($chromosome, 0);

    my($targets, $targets_stoich, $targets_inhomo) = make_DT($chromosome, 1);

    my $arr_targets = check_targets($targets, $targets_stoich, $g_reacs);
    my $arr_desired = check_targets($desired, $desired_stoich, $g_reacs);

    # append all to create dual system
    my $tmp1 = transpose($g_stoim);

    my $tmp2 = append_I_right($tmp1);

    my $tmp4 = append_target_cols($tmp2, $arr_targets);

    $dual_stoim = append_desired($tmp4, $g_stoim, $arr_desired);

    ($dual_metas, $dual_reacs, $dual_rever) = write_output_data_defigueiredo($g_reacs, $g_metas, $g_rever, $arr_targets, $arr_desired);

    my $cutsets;
    $cutsets = defigueiredo($dual_stoim, $dual_metas, $dual_reacs, $dual_rever, $targets_inhomo, $desired_inhomo);

    return $cutsets;
}
################################################################################

################################################################################
################################################################################
sub defigueiredo {
    my ($stoim, $metas, $reacs, $rever, $targets_inhomo, $desired_inhomo) = @_;
    my ($sfile,$mfile,$rfile,$rvfile,$ifile,$tfile,$efm_filename,$o_rvfile);
    my ($cplex_env,$lp);
    my ($solutions,$sol_norm);
    my $solution_cnt = 0;
    my $write_lp_file = 0;
    my $num_threads = 5;
    my ($start_networkTrans, $start_vdual, $start_wtarget, $start_desired);
    my ($stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired);
    my ($num_networkTrans, $num_vdual, $num_wtarget, $num_desired);
    my ($start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc);
    my ($stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc);
    my ($num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc);
    my ($num_desired_rev, $desired_map);
    my $cutsets;

    # initialising the above variables
    ($start_networkTrans, $stopp_networkTrans, $start_vdual, $stopp_vdual, $start_wtarget, $stopp_wtarget, $start_desired, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired) = read_reactions_dual($reacs);
    ($start_OrigReac, $stopp_OrigReac, $start_MetaDesiStoich, $stopp_MetaDesiStoich, $start_MetaDesiFunc, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc) = read_metabolites_dual($metas);
    ($num_desired_rev, $desired_map) = read_reversibility_dual($rever, $reacs);

    # doing defiguieredo
    $cutsets = do_figueiredo($stoim, $metas, $reacs, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $solutions,$sol_norm, $solution_cnt, $write_lp_file, $num_threads, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_desired_rev, $desired_map);

    return $cutsets;
}
################################################################################

################################################################################
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub do_figueiredo
{
    my($stoim, $metas, $reacs, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $solutions,$sol_norm, $solution_cnt, $write_lp_file, $num_threads, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_desired_rev, $desired_map) = @_;
    my($num_MIP_columns, $MIP_column_names);
    my $cutsets;

   # get CPLEX environment
   $cplex_env = Math::CPLEX::Env->new();
   die "ERROR: creating CPLEX environment failed!" unless $cplex_env;

   # create MIP problem
   $lp = $cplex_env->createOP();
   die "ERROR: couldn't create Linear Program\n" unless $lp;

   die "ERROR: minimize() failed\n" unless $lp->minimize();

   $cplex_env->setintparam(&Math::CPLEX::Base::CPX_PARAM_THREADS, $num_threads);
   $cplex_env->setdblparam(&Math::CPLEX::Base::CPX_PARAM_EPAGAP,  1e-10);

   fill_init_lp($stoim, $reacs, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $solutions,$sol_norm, $solution_cnt, $write_lp_file, $num_threads, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, \$num_MIP_columns, \$MIP_column_names, $num_desired_rev, $desired_map);
   #############################################################################
   # solve MIP problem
   #############################################################################
   die "ERROR: mipopt() failed\n" unless $lp->mipopt();
   #############################################################################

   #############################################################################
   # retrieve computed values
   #############################################################################
   if( $lp->getstat() == &Math::CPLEX::Base::CPXMIP_OPTIMAL )
   {
      my ($sol_status, $obj_val, @vals) = $lp->solution();
      my $cut;

      for( my $i = 0; $i < $num_MIP_columns; $i++ )
      {
         my $indicator_idx = $i - $stopp_networkTrans - 2*$num_vdual;
         if( $i >=  $stopp_networkTrans + 2*$num_vdual && $i < $stopp_networkTrans + 4*$num_vdual )
         {
            if( $vals[$i] > 0.5 )
            {
               my $idx = $i - 2*$num_vdual;
	       if($i >=  $stopp_networkTrans + 3*$num_vdual) { $idx = $i - 3*$num_vdual; }
               (my $reac_name = $reacs->[$idx]) =~ s/^vdual//;
               push @$cut, $reac_name;
            }
         }
      }
      push @$cutsets, $cut;
      
   }

   #############################################################################
   #############################################################################
   die "ERROR: free() failed\n" unless $lp->free();
   die "ERROR: close() failed\n" unless $cplex_env->close();
   #############################################################################

    # return the cutset
    return $cutsets;
}
################################################################################

################################################################################
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub add_solution
{
    my($stoim, $reacs, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $solutions,$sol_norm, $solution_cnt, $write_lp_file, $num_threads, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_MIP_columns, $MIP_column_names, $num_desired_rev, $desired_map) = @_;
   my $newRows;

   for( my $r = 0; $r < $num_MIP_columns; $r++ )
   {
      if( $solutions->[$solution_cnt][$r] > 0 )
      {
         $newRows->[0][$r] = 1.0;
      }
   }

   my $row = { num_rows  => 1,
               rhs       => [$sol_norm->[$solution_cnt] - 1],
               sense     => ['L'],
               row_names => ["sol_$solution_cnt"],
               row_coefs => $newRows};
               
   die "ERROR: addrows() failed\n" unless $lp->addrows($row);

   #############################################################################
   # write lp file
   #############################################################################
   if( $write_lp_file )
   {
      my $file_num = setLength($solution_cnt,4);
      my $filename = "/tmp/myCPLEX_maha_$file_num.lp";
      die "ERROR: writeprob() failed\n" unless $lp->writeprob($filename);
   }
   #############################################################################
}
################################################################################

################################################################################
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub fill_init_lp
{
    my($stoim, $reacs, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $solutions,$sol_norm, $solution_cnt, $write_lp_file, $num_threads, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_MIP_columns_ptr, $MIP_column_names_ptr, $num_desired_rev, $desired_map) = @_;
    my($num_MIP_columns, $MIP_column_names);
   #############################################################################
   # define columns
   #############################################################################
   my $obj_coefs;
   my $lower_bnd;
   my $upper_bnd;
   my $col_types;
   my $col_names;
   my $col_cnt = 0;

   # networkTrans
   for( my $i = $start_networkTrans; $i < $stopp_networkTrans; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = -1*&Math::CPLEX::Base::CPX_INFBOUND;
      $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
      $col_types->[$col_cnt] = 'C';
      $col_names->[$col_cnt] = $reacs->[$i];
      $col_cnt++;
   }

   # vdual -> positive
   for( my $i = $start_vdual; $i < $stopp_vdual; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = 0;
      $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
      $col_types->[$col_cnt] = 'C';
      $col_names->[$col_cnt] = $reacs->[$i] . 'pos';
      $col_cnt++;
   }

   # vdual -> negative
   for( my $i = $start_vdual; $i < $stopp_vdual; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = 0.0;
      $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
      $col_types->[$col_cnt] = 'C';
      $col_names->[$col_cnt] = $reacs->[$i] . 'neg';
      $col_cnt++;
   }

   # vdual -> indicator positive
   for( my $i = $start_vdual; $i < $stopp_vdual; $i++ )
   {
      $obj_coefs->[$col_cnt] = 1.0;
      $lower_bnd->[$col_cnt] = 0;
      $upper_bnd->[$col_cnt] = 1.0;
      $col_types->[$col_cnt] = 'B';
      $col_names->[$col_cnt] = 'ind' . $reacs->[$i] . 'pos';
      if(exists $g_excluded_recs->{$reacs->[$i]}) {
	  $upper_bnd->[$col_cnt] = 0.0;
      }
      $col_cnt++;
   }

   # vdual -> indicator negative
   for( my $i = $start_vdual; $i < $stopp_vdual; $i++ )
   {
      $obj_coefs->[$col_cnt] = 1.0;
      $lower_bnd->[$col_cnt] = 0.0;
      $upper_bnd->[$col_cnt] = 1.0;
      $col_types->[$col_cnt] = 'B';
      $col_names->[$col_cnt] = 'ind' . $reacs->[$i] . 'neg';
      if(exists $g_excluded_recs->{$reacs->[$i]}) {
	  $upper_bnd->[$col_cnt] = 0.0;
      }
      $col_cnt++;
   }

   # target
   for( my $i = $start_wtarget; $i < $stopp_wtarget; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = 0;
      $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
      $col_types->[$col_cnt] = 'C';
      $col_names->[$col_cnt] = $reacs->[$i];
      $col_cnt++;
   }

   # desired columns
   for( my $i = $start_desired; $i < $stopp_desired; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = 0.0;
      $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
      $col_types->[$col_cnt] = 'C';
      $col_names->[$col_cnt] = $reacs->[$i];
      $col_names->[$col_cnt] .= 'forward' if $rever->[$i] == 1;
      $col_cnt++;

      if( $rever->[$i] == 1 )
      {
         $obj_coefs->[$col_cnt] = 0.0;
         $lower_bnd->[$col_cnt] = 0.0;
         $upper_bnd->[$col_cnt] = &Math::CPLEX::Base::CPX_INFBOUND;
         $col_types->[$col_cnt] = 'C';
         $col_names->[$col_cnt] = $reacs->[$i] . 'backward';
         $col_cnt++;
      }
   }

   # desired columns -> indicator variable
   for( my $i = $start_desired; $i < $stopp_desired; $i++ )
   {
      $obj_coefs->[$col_cnt] = 0.0;
      $lower_bnd->[$col_cnt] = 0.0;
      $upper_bnd->[$col_cnt] = 1.0;
      $col_types->[$col_cnt] = 'B';
      $col_names->[$col_cnt] = 'ind' . $reacs->[$i];
      $col_names->[$col_cnt] .= 'forward' if $rever->[$i] == 1;
      $col_cnt++;

      if( $rever->[$i] == 1 )
      {
         $obj_coefs->[$col_cnt] = 0.0;
         $lower_bnd->[$col_cnt] = 0.0;
         $upper_bnd->[$col_cnt] = 1.0;
         $col_types->[$col_cnt] = 'B';
         $col_names->[$col_cnt] = 'ind' . $reacs->[$i] . 'backward';
         $col_cnt++;
      }
   }

   my $cols = { num_cols  => $col_cnt,
                obj_coefs => $obj_coefs,
                lower_bnd => $lower_bnd,
                upper_bnd => $upper_bnd,
                col_types => $col_types,
                col_names => $col_names
              };
   die "ERROR: newcols() failed\n" unless $lp->newcols($cols);

   $$num_MIP_columns_ptr = $col_cnt;
   $MIP_column_names_ptr = $col_names;
   $num_MIP_columns = $col_cnt;
   $MIP_column_names = $col_names;

   #############################################################################
   # add initial rows
   #############################################################################
   add_inital_rows($stoim, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_MIP_columns, $MIP_column_names, $num_desired_rev, $desired_map);
   #############################################################################

   #############################################################################
   # write lp file
   #############################################################################
   if( $write_lp_file )
   {
      my $filename = "/tmp/myCPLEX_maha.lp";
      die "ERROR: writeprob() failed\n" unless $lp->writeprob($filename);
   }
   #############################################################################
}
################################################################################


################################################################################
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub add_inital_rows()
{
    my($stoim, $rever, $targets_inhomo, $desired_inhomo, $cplex_env, $lp, $start_networkTrans, $start_vdual, $start_wtarget, $start_desired, $stopp_networkTrans, $stopp_vdual, $stopp_wtarget, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired, $start_OrigReac, $start_MetaDesiStoich, $start_MetaDesiFunc, $stopp_OrigReac, $stopp_MetaDesiStoich, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc, $num_MIP_columns, $MIP_column_names, $num_desired_rev, $desired_map) = @_;
   my $newRows;
   my $rhs;
   my $sense;
   my $row_names;
   my $row_cnt = 0;

   for( my $r = 0; $r < $stopp_OrigReac; $r++ )
   {
      $rhs->[$row_cnt]                      = 0.0;
      $sense->[$row_cnt]                    = 'E';
      if($g_rever->[$row_cnt] == 0) { $sense->[$row_cnt] = 'G'; }  #ineq
      $row_names->[$row_cnt]                = "origNet$r";
      for( my $c = $start_networkTrans; $c < $stopp_networkTrans; $c++ )
      {
         if( abs($stoim->[$r][$c]) > CONSIDER_ZERO )
         {
            $newRows->[$row_cnt][$c] = $stoim->[$r][$c];
         }
      }

      for( my $c = $start_vdual; $c < $stopp_vdual; $c++ )
      {
         if( abs($stoim->[$r][$c]) > CONSIDER_ZERO )
         {
            $newRows->[$row_cnt][$c]              = $stoim->[$r][$c];
            $newRows->[$row_cnt][$c + $num_vdual] = -1*$stoim->[$r][$c];
         }
      }

      for( my $c = $start_wtarget; $c < $stopp_wtarget; $c++ )
      {
         if( abs($stoim->[$r][$c]) > CONSIDER_ZERO )
         {
            $newRows->[$row_cnt][$c + 3*$num_vdual] = $stoim->[$r][$c];
         }
      }
      $row_cnt++;
   }

   for( my $r = $start_MetaDesiStoich; $r < $stopp_MetaDesiStoich; $r++ )
   {
      $rhs->[$row_cnt]                      = 0.0;
      $sense->[$row_cnt]                    = 'E';
      $row_names->[$row_cnt]                = "desiStoi$r";
      my $col_cnt = $start_desired;

      for( my $c = $start_desired; $c < $stopp_desired; $c++, $col_cnt++ )
      {
         if( abs($stoim->[$r][$c]) > CONSIDER_ZERO )
         {
            $newRows->[$row_cnt][$col_cnt + 3*$num_vdual] = $stoim->[$r][$c];

            if( $rever->[$c] == 1 )
            {
               $col_cnt++;
               $newRows->[$row_cnt][$col_cnt + 3*$num_vdual] = -1*$stoim->[$r][$c];
            }
         }
         else
         {
            if( $rever->[$c] == 1 )
            {
               $col_cnt++;
            }
         }
      }

      $row_cnt++;
   }

   my $desi_idx = 0;
   for( my $r = $start_MetaDesiFunc; $r < $stopp_MetaDesiFunc; $r++, $desi_idx++ )
   {
      $rhs->[$row_cnt]                      = $desired_inhomo->[$desi_idx];
      $sense->[$row_cnt]                    = 'L';
      $row_names->[$row_cnt]                = "desiFunc$r";
      my $col_cnt = $start_desired;

      for( my $c = $start_desired; $c < $stopp_desired; $c++, $col_cnt++ )
      {
         if( abs($stoim->[$r][$c]) > CONSIDER_ZERO )
         {
            $newRows->[$row_cnt][$col_cnt + 3*$num_vdual] = $stoim->[$r][$c];
            if( $rever->[$c] == 1 )
            {
               $col_cnt++;
               $newRows->[$row_cnt][$col_cnt + 3*$num_vdual] = -1*$stoim->[$r][$c];
            }
         }
         else
         {
            if( $rever->[$c] == 1 )
            {
               $col_cnt++;
            }
         }
      }

      $row_cnt++;
   }

   my $col_cnt = $start_desired;
   for( my $c = $start_desired; $c < $stopp_desired; $c++, $col_cnt++ )
   {
      if( $rever->[$c] == 1 )
      {
         $rhs->[$row_cnt]                      = 1.0;
         $sense->[$row_cnt]                    = 'L';
         $row_names->[$row_cnt]                = "PosPlusNeg$c";
      
         $newRows->[$row_cnt][$col_cnt + 3*$num_vdual  + $num_desired_rev + $num_desired] = 1;
         $col_cnt++;
         $newRows->[$row_cnt][$col_cnt + 3*$num_vdual  + $num_desired_rev + $num_desired] = 1;
         $row_cnt++;
      }
   }

   # sum of desired must be greater than 1
   $rhs->[$row_cnt]                      = 1.0;
   $sense->[$row_cnt]                    = 'G';
   $row_names->[$row_cnt]                = "avoidTrivDesi";
   for( my $c = 0; $c < $num_desired + $num_desired_rev; $c++ )
   {
        my $i = $c + $stopp_desired + $num_desired_rev + 3*$num_vdual;
        $newRows->[$row_cnt][$i] = 1;
   }
   $row_cnt++;

   # support of w
   $rhs->[$row_cnt]                        = -1.00;
   $sense->[$row_cnt]                      = 'L';
   $row_names->[$row_cnt]                  = "support_w";
   for( my $c = $start_wtarget; $c < $stopp_wtarget; $c++ )
   {
       my $t = $c - $start_wtarget;
       $newRows->[$row_cnt][$c + 3*$num_vdual] = $targets_inhomo->[$t];
   }
   $row_cnt++;

   # indicator_pos + indicator_pos <= 1, for all reversible reactions
   my $indicator_start = $stopp_networkTrans + 2*($stopp_vdual - $start_vdual);
   for( my $i = $indicator_start; $i < $indicator_start + $num_vdual; $i++ )
   {
      $rhs->[$row_cnt]                      = 1.0;
      $sense->[$row_cnt]                    = 'L';
      $row_names->[$row_cnt]                = "indposneg$i";
      $newRows->[$row_cnt][$i]              = 1.0;
      $newRows->[$row_cnt][$i + $num_vdual] = 1.0;
      $row_cnt++;
   }

   # sum of indicators <= cmcs_limit
    if(defined $g_cmcs_limit) {
	$rhs->[$row_cnt]                      = $g_cmcs_limit;
	$sense->[$row_cnt]                    = 'L';
	$row_names->[$row_cnt]                = "min_mcs_size";
	$indicator_start = $stopp_networkTrans + 2*($stopp_vdual - $start_vdual);
	for( my $i = $indicator_start; $i < $indicator_start + 2*$num_vdual; $i++ )
	{
	    $newRows->[$row_cnt][$i]              = 1.0;
	}
	$row_cnt++;
    }

   my $rows = {num_rows  => $row_cnt,
               rhs       => $rhs,
               sense     => $sense,
               row_names => $row_names,
               row_coefs => $newRows};

   die "ERROR: addrows() failed\n" unless $lp->addrows($rows);


   # add indicator constraints
   for( my $i = $start_vdual; $i < $start_vdual + 2*$num_vdual; $i++ )
   {
      my $val;
      $val->[$i] = 1.0;
      my $indconstr = {
                         indvar       => $i + 2*$num_vdual,
                         complemented => 1,
                         rhs          => 0.0,
                         sense        => 'L',
                         val          => $val,
                         name         => "indLess$i",
                      };

      die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
   }

   # add indicator constraints
   for( my $i = $start_vdual; $i < $start_vdual + 2*$num_vdual; $i++ )
   {
      my $val;
      $val->[$i] = 1.0;
      my $indconstr = {
                         indvar       => $i + 2*$num_vdual,
                         complemented => 0,
                         rhs          => CONSIDER_ZERO,
                         sense        => 'G',
                         val          => $val,
                         name         => "indGreater$i",
                      };

      die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
   }

   # add indicator constraints for desired network
   for( my $i = 0; $i < $num_desired_rev + $num_desired; $i++ )
   {
      my $val;
      my $idx = $i + $start_desired + 3*$num_vdual;
      $val->[$idx] = 1.0;
      my $indconstr = {
                         indvar       => $idx + $num_desired_rev + $num_desired,
                         complemented => 1,
                         rhs          => 0.0,
                         sense        => 'L',
                         val          => $val,
                         name         => "indLessDesi$i",
                      };

      die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
   }

   # add indicator constraints for desired network
   for( my $i = 0; $i < $num_desired_rev + $num_desired; $i++ )
   {
      my $val;
      my $idx = $i + $start_desired + 3*$num_vdual;
      $val->[$idx] = 1.0;
      my $indconstr = {
                         indvar       => $idx + $num_desired_rev + $num_desired,
                         complemented => 0,
                         rhs          => CONSIDER_ZERO,
                         sense        => 'G',
                         val          => $val,
                         name         => "indGreaterDesi$i",
                      };

      die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
   }

   # add indicator for MCS acting on desired
   for( my $i = 0; $i < $num_vdual; $i++ )
   {
      foreach my $idx ( @{$desired_map->[$i]} )
      {
         my $val;
         $val->[$start_desired + 3*$num_vdual + $idx] = 1.0;
         my $indconstr = {
                            indvar       => $i + $start_vdual + 2*$num_vdual,
                            complemented => 0,
                            rhs          => 0.0,
                            sense        => 'E',
                            val          => $val,
                            name         => "indMCSpos$i$idx",
                         };
         die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
      }

      foreach my $idx ( @{$desired_map->[$i]} )
      {
         my $val;
         $val->[$start_desired + 3*$num_vdual + $idx] = 1.0;
         my $indconstr = {
                            indvar       => $i + $start_vdual + 3*$num_vdual,
                            complemented => 0,
                            rhs          => 0.0,
                            sense        => 'E',
                            val          => $val,
                            name         => "indMCSneg$i$idx",
                         };

            die "ERROR: addindconstr() failed\n" unless $lp->addindconstr($indconstr);
      }
   }
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub write_output_data_defigueiredo
{
   my $rcs      = shift;
   my $mts      = shift;
   my $rvs      = shift;
   my $tgs      = shift;
   my $des      = shift;

   # dual system: reactions are not metabolites!
   my $new_metas;
   foreach my $origReac (@$rcs )
   {
      push @$new_metas, "OrigReac" . $origReac;
   }
   foreach my $meta_des (@$mts )
   {
      push @$new_metas, "MetaDesi" . $meta_des;
   }
   for( my $i = 0; $i < @$des; $i++ )
   {
      push @$new_metas, "Desired" . setLength($i,2);
   }

   # create reaction and reversibility files of dual system
   my $new_reacs;
   my $new_rever;
   for( my $i = 0; $i < @$mts; $i++ )
   {
      push @$new_reacs, 'networkTrans' .$mts->[$i];
      push @$new_rever, '1';
   }
   for( my $i = 0; $i < @$rcs; $i++ )
   {
      push @$new_reacs, 'vdual' . $rcs->[$i];
      push @$new_rever, '1';
   }
   for( my $i = 0; $i < @$tgs; $i++ )
   {
      push @$new_reacs, 'wtarget' . set_length($i,2);
      push @$new_rever, '0';
   }

   for( my $i = 0; $i < @$rcs; $i++ )
   {
      push @$new_reacs, 'desired' . $rcs->[$i];
      if($rvs->[$i] )
      {
         push @$new_rever, '1';
      }
      else
      {
         push @$new_rever, '0';
      }
   }

   return ($new_metas, $new_reacs, $new_rever);
}
################################################################################

################################################################################
# one reaction in target and rest in desired
################################################################################
sub make_DT {
    my($ratios, $if_target) = @_;
    my $desired;
    my $desired_stoich;
    my $desired_inhomo;

    # the stoichiometry for the norm recs
    my $norm_so;
    for(0 .. @$g_norm_recs){
	push @$norm_so, -1;
    }
    push @$desired, $g_norm_recs;
    push @$desired_stoich, $norm_so;
    push @$desired_inhomo, -1;

    if($if_target) {
	my $recs;
	my $so;
	foreach my $norm_rec(@$g_norm_recs){
	    my $neg_rat = -1*$ratios->[0];
	    push @$recs, $norm_rec;
	    push @$so, $neg_rat;
	}
	push @$recs, $g_objective_reactions->[0];
	push @$so, 1;
	push @$desired, $recs;
	push @$desired_stoich, $so;
	push @$desired_inhomo, 0;
    } else {
	for(my $i = 0; $i < @$g_objective_reactions; $i++){
	    my $recs;
	    my $so;
	    foreach my $norm_rec(@$g_norm_recs){
		push @$recs, $norm_rec;
		push @$so, $ratios->[$i];
	    }
	    push @$recs, $g_objective_reactions->[$i];
	    push @$so, -1;
	    push @$desired, $recs;
	    push @$desired_stoich, $so;
	    push @$desired_inhomo, 0;
	}
    }

    return $desired, $desired_stoich, $desired_inhomo;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub append_target_cols
{
   my $in  = shift;
   my $tgs = shift;
   my $out;

   for( my $m = 0; $m < @$in; $m++ )
   {
      for( my $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         $out->[$m][$r] = $in->[$m][$r];
      }
   }

   my $c = @{$in->[0]};

   for( my $t = 0; $t < @$tgs; $t++ )
   {
      for( my $u = 0; $u < @{$tgs->[$t]}; $u++ )
      {
         $out->[$u][$c] = $tgs->[$t][$u];
      }
      $c++;
   }

   return $out;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub append_desired
{
   my $in  = shift;
   my $stm = shift;
   my $ds  = shift;
   my $out;

   for( my $m = 0; $m < @$in; $m++ )
   {
      for( my $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         $out->[$m][$r] = $in->[$m][$r];
      }
   }

   my $c_store = @{$in->[0]};
   my $r_store = @{$in};

   # add zero at bottom for N
   for( my $r = 0; $r < @$stm + @$ds; $r++ )
   {
      for( my $c = 0; $c < $c_store; $c++ )
      {
         $out->[$r_store + $r][$c] = 0;
      }
   }

   for( my $r = 0; $r < $r_store; $r++ )
   {
      for( my $c = 0; $c < @{$stm->[0]}; $c++ )
      {
         $out->[$r][$c + $c_store] = 0;
      }
   }

   for( my $r = 0; $r < @$stm; $r++ )
   {
      for( my $c = 0; $c < @{$stm->[$r]}; $c++ )
      {
         $out->[$r + $r_store][$c + $c_store] = $stm->[$r][$c];
      }
   }

   for( my $r = 0; $r <  @$ds; $r++ )
   {
      for( my $c = 0; $c < @{$ds->[$r]}; $c++ )
      {
         $out->[$r + $r_store + @$stm][$c + $c_store] = $ds->[$r][$c];
      }
   }


   return $out;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub append_I_right
{
   my $in  = shift;
   my $out = [];

   # copy upper part
   my $m = 0;
   my $r = 0;
   for( ; $m < @$in; $m++ )
   {
      # print "m=$m\n";
      for( $r = 0; $r < @{$in->[$m]}; $r++ )
      {
         $out->[$m][$r] = $in->[$m][$r];
      }
   }
   my $last = $r;

   for( $m = 0; $m < @$in; $m++ )
   {
      for( my $r = 0; $r < @$in; $r++ )
      {
         $out->[$m][$r+$last] = 0.0;
         $out->[$m][$r+$last] = 1.0 if $m == $r;
      }
   }

   return $out;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub transpose
{
   my $mat = shift;
   my $trp = [];;

   # number of rows
   my $m = @$mat;
   # number of colums
   my $n = @{$mat->[0]};

   for( my $r = 0; $r < $m; $r++ )
   {
      for( my $c = 0; $c < $n; $c++ )
      {
         $trp->[$c][$r] = $mat->[$r][$c];
      }
   }

   return $trp;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub print_vector
{
   my $mes = shift;
   my $vec = shift;

   print $mes;

   print join(",",@$vec), "\n"
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub print_matrix
{
   my $message = shift;
   my $matrix  = shift;

   warn $message;

   for( my $m = 0; $m < @$matrix; $m++ )
   {
      for( my $r = 0; $r < @{$matrix->[$m]}; $r++ )
      {
         print STDERR "\t$matrix->[$m][$r]";
      }
      warn "\n";
   }
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub check_targets
{
   my $rcs_map = $g_rxnMap;
   my $tgs     = shift;
   my $tgs_sto = shift;
   my $rcs_nm  = shift;
   my $t_arr;

   foreach my $target_line (@$tgs)
   {
      foreach my $target (@$target_line)
      {
         die "FATAL ERROR: $target not found in set of reactions specified in rfile\n" unless defined $rcs_map->{$target};
      }
   }

   for( my $x = 0; $x < @$tgs; $x++ )
   {
      for( my $r = 0; $r < @$rcs_nm; $r++ )
      {
         $t_arr->[$x][$r] = 0;
      }

      for( my $t = 0; $t < @{$tgs->[$x]}; $t++ )
      {
         my $r;
         my $found_reacs = 0;
         for( $r = 0; $r < @$rcs_nm; $r++ )
         {
            if( $rcs_nm->[$r] eq $tgs->[$x][$t] )
            {
                $found_reacs = 1;
                last;
            }
         }

         if( $found_reacs )
         {
             $t_arr->[$x][$r] = $tgs_sto->[$x][$t];
         }
         else
         {
            die "FATAL ERROR: didn't find target '$tgs->[$x][$t]' in array of known reactions\n";
         }
      }
   }

   return $t_arr;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub set_length
{
   my $in  = shift;
   my $num = shift;

   while( $num > length($in) )
   {
      $in = '0' . $in;
   }

   return $in;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
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
# from Christian Jungreuthmayer
################################################################################
sub read_stoichiomat
{
   my $file = shift;
   my $sm;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";

   while( <$fh> )
   {
      s/^\s+//;
      s/\s+$//;
      my @st_facs = split;
      push @$sm, \@st_facs;
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
################################################################################
sub read_reactions_2_hash
{
   my $file = shift;
   my $rcs;
   my $rcs_h;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";
   while(<$fh>) {
       chomp;
       s/"//g;
       s/'//g;
       s/>//g;
       s/#//g;
       push @$rcs, $_;
   }
   close $fh;

   %$rcs_h = map {"vdual".$_ => 0} @$rcs;

   return $rcs_h;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub print_sfile
{
   my $name = shift;
   my $mat  = shift;

   open my $fh, ">$name" or die "ERROR: couldn't open file '$name' for writing: $!\n";

   for( my $m = 0; $m < @$mat; $m++ )
   {
      print $fh "@{$mat->[$m]}\n";
   }

   close $fh;
}
################################################################################

################################################################################
# from Christian Jungreuthmayer
################################################################################
sub print_file
{
   my $name = shift;
   my $rcs  = shift;
   my $add_double_quotes = shift;

   open my $fh, ">$name" or die "ERROR: couldn't open file '$name' for writing: $!\n";

   for( my $m = 0; $m < @$rcs; $m++ )
   {
      if( $add_double_quotes )
      {
         print $fh "\"$rcs->[$m]\"";
      }
      else
      {
         print $fh "$rcs->[$m]";
      }
      print $fh " " unless $m == @$rcs - 1;
   }

   close $fh;
}
################################################################################

################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub read_reversibility_dual
{
   my ($rvs, $reacs) = @_;

   my $des_cnt = 0;
   my $num_desired_rev = 0;
   my $desired_map;

   for( my $r = 0; $r < @$rvs; $r++ )
   {
      if( $reacs->[$r] =~ /^desired/ )
      {
         my $arr_ref;
         push @$arr_ref, $des_cnt;
         $des_cnt++;
         if( $rvs->[$r] == 1 )
         {
            push @$arr_ref, $des_cnt;
            $des_cnt++;
            $num_desired_rev++;
         }

         push @$desired_map, $arr_ref;
      }
   }

   return ($num_desired_rev, $desired_map);
}
################################################################################

################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub read_metabolites_dual
{
   my $found_OrigReac       = 0;
   my $found_MetaDesiStoich = 0;
   my $found_MetaDesiFunc   = 0;
   my $mbs = shift;
   my ($start_OrigReac, $stopp_OrigReac, $start_MetaDesiStoich, $stopp_MetaDesiStoich, $start_MetaDesiFunc, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc);

   for( my $m = 0; $m < @$mbs; $m++ )
   {
      if( $mbs->[$m] =~ /^OrigReac/ && $found_OrigReac == 0)
      {
         $found_OrigReac = 1;
         $start_OrigReac = $m;
      }
      if( $mbs->[$m] =~ /^MetaDesi/ && $found_MetaDesiStoich == 0)
      {
         $found_MetaDesiStoich = 1;
         $start_MetaDesiStoich = $m;
      }
      if( $mbs->[$m] =~ /^Desired/ && $found_MetaDesiFunc == 0)
      {
         $found_MetaDesiFunc = 1;
         $start_MetaDesiFunc = $m;
      }

      if( $found_OrigReac == 1 && $mbs->[$m] !~ /^OrigReac/ )
      {
         $found_OrigReac = 0;
         $stopp_OrigReac = $m;
      }
      if( $found_MetaDesiStoich == 1 && $mbs->[$m] !~ /^MetaDesi/ )
      {
         $found_MetaDesiStoich = 0;
         $stopp_MetaDesiStoich = $m;
      }
      if( $found_MetaDesiFunc == 1 && $mbs->[$m] !~ /^Desired/ )
      {
         die "FATAL ERROR: We expect that the last metabolite is a desired reaction. Hence, we should never end up in here!\n";
      }
   }
   $stopp_MetaDesiFunc = @$mbs;

  $num_OrigReac       = $stopp_OrigReac       - $start_OrigReac;
  $num_MetaDesiStoich = $stopp_MetaDesiStoich - $found_MetaDesiStoich;
  $num_MetaDesiFunc   = $stopp_MetaDesiFunc   - $start_MetaDesiFunc;

   return ($start_OrigReac, $stopp_OrigReac, $start_MetaDesiStoich, $stopp_MetaDesiStoich, $start_MetaDesiFunc, $stopp_MetaDesiFunc, $num_OrigReac, $num_MetaDesiStoich, $num_MetaDesiFunc);
}
################################################################################

################################################################################
# reaction file has one line that contains a list of the reaction names
# the reaction names are separated by white spaces
# modified from Christian Jungreuthmayer's defigueiredo.pl
################################################################################
sub read_reactions_dual
{
   my $found_networkTrans = 0;
   my $found_vdual        = 0;
   my $found_wtarget      = 0;
   my $found_desired      = 0;
   my $rcs = shift;
   my ($start_networkTrans, $stopp_networkTrans, $start_vdual, $stopp_vdual, $start_wtarget, $stopp_wtarget, $start_desired, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired);

   for( my $r = 0; $r < @$rcs; $r++ )
   {
      if( $rcs->[$r] =~ /^networkTrans/ && $found_networkTrans == 0)
      {
         $found_networkTrans = 1;
         $start_networkTrans = $r;
      }
      if( $rcs->[$r] =~ /^vdual/ && $found_vdual == 0)
      {
         $found_vdual = 1;
         $start_vdual = $r;
      }
      if( $rcs->[$r] =~ /^wtarget/ && $found_wtarget == 0)
      {
         $found_wtarget = 1;
         $start_wtarget = $r;
      }
      if( $rcs->[$r] =~ /^desired/ && $found_desired == 0)
      {
         $found_desired = 1;
         $start_desired = $r;
      }

      if( $found_networkTrans == 1 && $rcs->[$r] !~ /^networkTrans/ )
      {
         $found_networkTrans = 0;
         $stopp_networkTrans = $r;
      }
      if( $found_vdual == 1 && $rcs->[$r] !~ /^vdual/ )
      {
         $found_vdual = 0;
         $stopp_vdual = $r;
      }
      if( $found_wtarget == 1 && $rcs->[$r] !~ /^wtarget/ )
      {
         $found_wtarget = 0;
         $stopp_wtarget = $r;
      }
      if( $found_desired == 1 && $rcs->[$r] !~ /^desired/ )
      {
         die "FATAL ERROR: We expect that the last reaction is a desired reaction. Hence, we should never end up in here!\n";
      }
   }
   $stopp_desired = @$rcs;

   $num_networkTrans = $stopp_networkTrans - $start_networkTrans;
   $num_vdual        = $stopp_vdual        - $start_vdual;
   $num_wtarget      = $stopp_wtarget      - $start_wtarget;
   $num_desired      = $stopp_desired      - $start_desired;

   return ($start_networkTrans, $stopp_networkTrans, $start_vdual, $stopp_vdual, $start_wtarget, $stopp_wtarget, $start_desired, $stopp_desired, $num_networkTrans, $num_vdual, $num_wtarget, $num_desired);
}
################################################################################

################################################################################
# read in program options
################################################################################
sub read_arguments
{
   getopts('s:m:r:v:e:o:p:h');

   if( $opt_h )
   {
      usage();
   }

   if( $opt_e )
   {
      $efile = $opt_e;
   }

   if( $opt_p )
   {
      $parameters_file = $opt_p;
   }
   else
   {
      usage('ERROR: please provide an input file containing pso parameters ',-1);
   }

   if( $opt_s )
   {
      $sfile = $opt_s;
   }
   else
   {
      usage('ERROR: please provide an input file containing stoichiometric matrix ',-1);
   }

   if( $opt_m )
   {
      $mfile = $opt_m;
   }
   else
   {
      usage('ERROR: please provide an input file containing metabolite names ',-1);
   }

   if( $opt_r )
   {
      $rfile = $opt_r;
   }
   else
   {
      usage('ERROR: please provide an input file containing reaction names ',-1);
   }

   if( $opt_v )
   {
      $rvfile = $opt_v;
   }
   else
   {
      usage('ERROR: please provide an input file containing reaction reversibility information ',-1);
   }

   if( $opt_o )
   {
      $g_ofile = $opt_o;
   }

}
################################################################################

################################################################################
################################################################################
sub usage
{
   my $message   = shift || '';
   my $exit_code = shift || 0;

   print "$message\n" if $message;

   print "psomcs.pl -s sfile -m file -r rfile -v rvfile -e excluded_reactions -p pso_parameters -o output_file [-h]\n";
   print "\n";
   print "-s ..... name of file containing stoichiometric matrix (input)\n";
   print "-m ..... name of file containing metabolites (input)\n";
   print "-r ..... name of file containing reactions (input)\n";
   print "-v ..... name of file containing reversibility information (input)\n";
   print "-e ..... name of file containing reactions to be excluded from knockouts(input) \n";
   print "-p ..... name of file containing pso parameters (input)\n";
   print "-o ..... program output file\n";
   print "-h ..... print this message\n";
   print "\n";
   print "psomcs.pl finds optimal knockouts in metabolic networks\n";

   exit($exit_code);
}
################################################################################
