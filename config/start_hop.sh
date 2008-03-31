#!/bin/sh

#PBS -m ae
#PBS -M <email_adress>
#PBS -N <job_name>
#PBS -l nodes=2:ppn=2
#mkdir -p /scratch/$USER/benchmark
#cd /scratch/$USER/test


srcdir=$PBS_O_WORKDIR
cd $srcdir

PROCS=4
PREFIX=<outputfile_with_path>
EXE=<exename>
INPUT=<inputfile_with_path>
RESTART=
echo $PBS_NODEFILE
VAPILIB=/usr/local/ibgd/driver/infinihost/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VAPILIB

# For bugfixing use newer version of openmpi!
MPIDIR=/cluster/openmpi-1.1.4
# give path for gnu 4.2.2 libs
export CPPLIB=/cluster/gcc-4.2.2/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CPPLIB

$MPIDIR/bin/mpirun -np $PROCS -hostfile $PBS_NODEFILE $srcdir/$EXE $INPUT $PREFIX $RESTART | tee $PREFIX.log
