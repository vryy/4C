#!/bin/sh
#PBS -N ca_test
#PBS -l nodes=1:ppn=2,cput=0:30:00,mem=500mb

HOMEDIR=$HOME/<myname>/ccarat
SCRATCHDIR=$SCRDIR/gee
cp $HOMEDIR/<exe> $SCRATCHDIR/.
cp $HOMEDIR/<inputfile> $SCRATCHDIR/.
cd $SCRATCHDIR
#
date
echo  "start application:"
#get_input job
/opt/NECmpi/bin/mpirun -np 2 <exe> <inputfile xxx > xxx.stdout
#save_result data
echo "job finished"
date
