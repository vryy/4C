#!/bin/bash
#PBS -o /home/hlrb2/pr32ne/lu43rob/output/baci.out
#PBS -j oe
#PBS -N baci
#PBS -S /bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=4:ncpus=1:mem=8gb
#PBS -M user@lnm.mw.tum.de
#PBS -m abe
#PBS -r n
. /etc/profile.d/modules.sh
BACIDIR=$HOME/baci
INPUTDIR=$HOME/input
OUTPUTDIR=$HOME/output
PREFIX=filename
echo $OPT_TMP
mpiexec $BACIDIR/cca_suse_mpi.altix.fast $INPUTDIR/$PREFIX.dat $OPT_TMP/zzz
cp $OPT_TMP/zzz* $OUTDIR/


