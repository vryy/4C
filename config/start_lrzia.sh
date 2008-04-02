#!/bin/bash 
#$-o $HOME/gee/ccarat/baci.$JOB_ID.out -j y 
#$-N baci
#$-l h_rt=00:10:00 
#$-S /bin/bash 
#$-l march=ia64 
#$-pe mpi 4
#$-M gee@lnm.mw.tum.de
#$-m abe
#$-cwd
#$-r n
. /etc/profile
CURRDIR=$HOME/gee/ccarat
mpiexec -n $NSLOTS $CURRDIR/cca_ia_parastation.fast $CURRDIR/cylinder_3d_001.dat $OPT_TMP/zzz
cp $OPT_TMP/zzz* $CURRDIR/.
jobperf -j $JOB_ID