#!/bin/bash 
#$-o $HOME/gee/ccarat/baci.$JOB_ID.out -j y 
#$-N baci
#$-l h_rt=00:30:00 
#$-S /bin/bash 
#$-l march=x86_64
## mpi_8: multiples of 8, machine size is 38 times 8
## mpi_16: multiples of 16, machine size is 15 times 16 (machine owned by TU Mathematics and might go away)
#$-pe mpi_8 8
#$-M <name>@lnm.mw.tum.de
#$-m abe
#$-cwd
#$-r n
. /etc/profile
CURRDIR=$HOME/gee/ccarat
echo $NSLOTS
mpiexec -n $NSLOTS $CURRDIR/cca_amd_parastation.fast $CURRDIR/cylinder_3d_001.dat $OPT_TMP/xxx
cp $OPT_TMP/zzz* $CURRDIR/.
jobperf -j $JOB_ID
