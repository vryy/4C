
# build the Makefile
# This is done in three steps

# step 1: store all user supplied variables
#
cat > $makefile <<EOF
#
# Variables
#
SRC=$SRC

PLAIN_CC=$PLAIN_CC
PLAIN_F77=$PLAIN_F77
PLAIN_LD=$PLAIN_LD

CC=$CC
F77=$F77
LD=$LD

PROGRAM=$PROGRAMNAME

CFLAGS=$CFLAGS $DEBUGFLAG -D$PLATFORM $DEFINES
FFLAGS=$FFLAGS $DEBUGFLAG -D$PLATFORM $DEFINES

LDFLAGS=$LDFLAGS
INCLUDES=$INCLUDEDIRS
LIBS=$LIBDIRS $LIBS
PLAIN_LIBS=$LIBDIRS $PLAIN_LIBS

VISUAL2_LIB=$VISUAL2_LIB
#VISUAL2_INC=$VISUAL2_INC

#----------------------- binaries -----------------------------------
include ./Makefile.objects
#--------------------------------------------------------------------
#
# The main rule called when no arguments are given
ccarat: \$(PROGRAM)

# Build (nearly) everything.
# Some filters (like the visual ones) are very specific and system dependent
# and thus not included here. This rule is supposed to work everywhere.
all: \$(PROGRAM) post_gid_txt post_out post_monitor

\$(PROGRAM): \\
		$OBJECTS
		@echo "Linking \$(LD) \$(LDFLAGS) \$(LIBS) \$(INCLUDES) -o \$(PROGRAM)"
		@\$(LD) \$(LDFLAGS) \\
		$OBJECTS \\
		\$(LIBS)  -o \$(PROGRAM)
		@echo "\$(PROGRAM) successfully built"
EOF


# step 2: copy Makefile.in
cat Makefile.in >> $makefile
#awk '$1 !~ /include/ { print $0 } $1 ~ /include/ { system("cat " $2)}' Makefile.in >> Makefile

# step 3: build dependencies if gcc can be found
#
# you can skip this step by saying
# $ NODEPS=yes ./configure ...
#
if [ x$NODEPS != "xyes" ] ; then
  # hopefully nobody translates which's messages...
  if which gcc | grep '^no ' 2>&1 > /dev/null ; then
     echo $0: gcc not found. No dependencies generated. Use Makefile with care.
  else
    for file in `find src -name "*.c"` ; do
      echo "build deps for" $file
      gcc -D$PLATFORM $DEFINES -MM -MT `echo $file|sed -e 's/c$/o/'` $INCLUDEDIRS $file >> $makefile
    done
  fi
fi
