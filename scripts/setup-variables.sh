
# default values
SRC=.

# for filters we need non-mip versions of these
PLAIN_CC=cc
PLAIN_F77=f77
PLAIN_LD=f77

# for the filters we do NOT nedd the solver libs
PLAIN_LIBS=$LIBS

# most of the time these have to begin changed
CC=cc
F77=f77
LD=f77

# 64 bit exe?
if grep '^SIXTYFOUR' "$definefile" 2>&1 > /dev/null ; then
   BIT=64
else
   BIT=32   
fi


# source config file
# Here we get all the variables.
. "./$configfile"


# get the define flags

DEFINES=`sed -e 's/#.*//' "$definefile" | awk 'BEGIN { defs = "" }
END { print defs }
  { if (length($1) > 0) { defs = defs" -D"$1 } }'`

#echo ">>>" $DEFINES "<<<"
#exit

# add special cflags from command line
if [ "x$specialargs" != "x" ] ; then
    echo "Found special cflags on the command line !!"
    CFLAGS="$CFLAGS $specialargs"
    FFLAGS="$FFLAGS $specialargs"
fi


# the debug version needs special treatment
if grep '^DEBUG' "$definefile" 2>&1 > /dev/null ; then
    DEBUGFLAG=-g
    PROGRAMNAME=$PROGRAMNAME.debg
else
    DEBUGFLAG=
    PROGRAMNAME=$PROGRAMNAME.fast

    CFLAGS="$CFLAGS $CFLAGS_OPT"
    FFLAGS="$FFLAGS $FFLAGS_OPT"
    LDFLAGS="$LDFLAGS $LDFLAGS_OPT"
fi

# a parallel version needs one more compile time flag
if [ x$PARALLEL = "xyes" ] ; then
    DEFINES=" -DPARALLEL $DEFINES"
fi
