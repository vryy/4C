#!/bin/sh
#PBS -N ca_test
#PBS -l nodes=1:ppn=2,cput=0:30:00,mem=500mb


WORKDIR=$SCRDIR/gee
cp ./cca_ita1_debg.exe $WORKDIR/gee/.
cp ./dyn_shell_dirich.dat $WORKDIR/gee/.
cd $WORKDIR
#
date
echo  "start application:"
#get_input job
#mpirun -np 2 cca_ita1_debg.exe dyn_shell_dirich.dat xxx
#save_result data
echo "job finished"
date
