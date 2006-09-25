
# check the different solvers

# trilinos
if grep '^[ \t]*TRILINOS_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$TRILINOS_LIB" = "x" ] ; then
        echo $0: Warning: Variable TRILINOS_LIB undefined but TRILINOS_PACKAGE requested.
    fi
    LIBS="-L$TRILINOS_LIBDIR $TRILINOS_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS -I$TRILINOS_INCDIR"
fi

# aztec (only if no trilinos, bcause trilinos contains aztec anyway)
if grep '^[ \t]*TRILINOS_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
  IS_TRILINOS=1
else
  if grep '^[ \t]*AZTEC_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
      if [ "x$AZTEC_LIB" = "x" ] ; then
          echo $0: Warning: Variable AZTEC_LIB undefined but AZTEC_PACKAGE requested.
      fi
      LIBS="$AZTEC_LIB $LIBS"
      INCLUDEDIRS="$INCLUDEDIRS $AZTEC_INC"
  fi
fi

# spooles
if grep '^[ \t]*SPOOLES_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$SPOOLES_LIB" = "x" ] ; then
        echo $0: Warning: Variable SPOOLES_LIB undefined but SPOOLES_PACKAGE requested.
    fi
    LIBS="$SPOOLES_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS $SPOOLES_INC"
fi

# umfpack
if grep '^[ \t]*UMFPACK' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$UMFPACK_LIB" = "x" ] ; then
        echo $0: Warning: Variable UMFPACK_LIB undefined but UMFPACK_PACKAGE requested.
    fi
    LIBS="$UMFPACK_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS $UMFPACK_INC"
fi

# nurbs
if grep '^[ \t]*NURBS' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$NURBS_LIB" = "x" ] ; then
        echo $0: Warning: Variable NURBS_LIB undefined but NURBS_PACKAGE requested.
    fi
    LIBS="$NURBS_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS $NURBS_INC"
fi


# visual2
if grep '^[ \t]*VISUAL2_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$VISUAL2_LIB" = "x" ] ; then
        echo $0: Error: Variable VISUAL2_LIB undefined but VISUAL2_PACKAGE requested.
        exit 1
    fi
    if [ x$BIT = "x64" ] ; then
        echo $0: Error: VISUAL2 not possible in 64 bit executable.
        exit 1
    fi
    LIBS="$VISUAL2_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS $VISUAL2_INC"
fi

# visual3
if grep '^[ \t]*VISUAL3_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
    if grep '^VISUAL2_PACKAGE' "$definefile" 2>&1 > /dev/null ; then
       echo $0: Error: VISUAL2 & VISUAL3 cannot linked into the same executable.
       exit 1
    fi
    if [ "x$VISUAL3_LIB" = "x" ] ; then
        echo $0: Error: Variable VISUAL3_LIB undefined but VISUAL3_PACKAGE requested.
        exit 1
    fi
    if [ x$BIT = "x64" ] ; then
        echo $0: Error: VISUAL3 not possible in 64 bit executable.
        exit 1
    fi
    LIBS="$VISUAL3_LIB $LIBS"
    INCLUDEDIRS="$INCLUDEDIRS $VISUAL2_INC"
fi
