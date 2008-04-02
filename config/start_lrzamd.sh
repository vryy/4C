#!/bin/bash 
#$-o $HOME/gee/ccarat/baci.$JOB_ID.out -j y 
#$-N baci
#$-l h_rt=00:30:00 
#$-S /bin/bash 
#$-l march=x86_64
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
