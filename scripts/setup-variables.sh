
# default values
MPIBOOT=
MPIHALT=
MPIRUN=mpirun
MPIRUNARGS=

# for filters we need non-mip versions of these
PLAIN_CC=cc
PLAIN_CP=c++
PLAIN_F77=f77
PLAIN_LD=f77

# for the filters we do NOT need the solver libs
PLAIN_LIBS=$LIBS

# filters based on drt will need trilinos and metis as well.
FILTER_LIBS=$LIBS

# most of the time these have to begin changed
CC=cc
CP=CC
F77=f77
LD=f77

# 64 bit exe?
if grep '^SIXTYFOUR' "$definefile" 2>&1 > /dev/null ; then
   BIT=64
else
   BIT=32
fi


# C++ is used => us CP for linking
USES_CPP=0


# check for DEBUG version
IS_DEBUG=0

# check for trilinos usage
IS_TRILINOS=0
if grep '^[[:blank:]]*TRILINOS_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
  IS_TRILINOS=1
fi

# - from definefile
if grep '^[[:blank:]]*DEBUG' "$definefile" 2>&1 > /dev/null ; then
  echo DEBUG from definefile
  IS_DEBUG=1
fi

# - from command line
if [ "x$extraargs" = "xDEBUG" ] ; then
  echo DEBUG from command line
  IS_DEBUG=1
fi

# - from environment variable
if [ x$DEBUG = "xyes" ] ; then
  echo DEBUG from env
  IS_DEBUG=1
fi


# source config file
# Here we get all the variables.
. "$configfile"


# get the define flags
DEFINES=`sed -e 's/#.*//' -e 's/^DEBUG//' "$definefile" | awk 'BEGIN { defs = "" }
END { print defs } { if ( $NF ) { if ( $1 != "" ) { defs = defs" -D"$1 } } }'`

#echo ">>>" $DEFINES "<<<"
#exit

# if we use trilinos, we have to add some more trilinos defines
if [ x$IS_TRILINOS = "x1" ] ; then
  DEFINES="$DEFINES $TRILINOS_DEFINES"
fi

# always look into headers directory of destination
CFLAGS="-I$DEST/src/headers $CFLAGS"
CPFLAGS="-I$DEST/src/headers $CPFLAGS"

# the debug version needs special treatment
if [ x$IS_DEBUG = "x1" ] ; then
  DEBUGFLAG=-g
  DEFINES=" -DDEBUG $DEFINES"
  PROGRAMNAME=$PROGRAMNAME.debg
else
  DEBUGFLAG=
  PROGRAMNAME=$PROGRAMNAME.fast

  CFLAGS="$CFLAGS $CFLAGS_OPT"
  CPFLAGS="$CPFLAGS $CPFLAGS_OPT"
  FFLAGS="$FFLAGS $FFLAGS_OPT"
  LDFLAGS="$LDFLAGS $LDFLAGS_OPT"
fi


# a parallel version needs one more compile time flag
if [ x$PARALLEL = "xyes" ] ; then
    DEFINES=" -DPARALLEL $DEFINES"
fi

# define the length of the loops for fast elements
DEFINES="$DEFINES $LOOPL"

#echo ">>>" $DEFINES "<<<"
#exit
