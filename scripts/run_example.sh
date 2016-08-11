#!/bin/bash
echo "running psomcs.pl ....."
perl psomcs.pl -s ../examples/trinh.sfile -r ../examples/trinh.rfile -m ../examples/trinh.mfile -v ../examples/trinh.rvfile -e ../examples/trinh.efile -p ../examples/parameters.txt -o ../examples/trinh_pso_out.txt
echo " ..... psomcs.pl has finished calculating cMCSs"
