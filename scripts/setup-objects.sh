
# fluid objects:
OBJ_FLUID=0

# control-flag, to ensure that OBJS_S2ML get OBJS_WALL,OBJS_WALLGE,OBJS_IF:
OBJ_W1=0
OBJ_WGE=0
OBJ_IF=0


# objects:
OBJECTS="\$(OBJS_MAIN) \$(OBJS_GLOBAL) \$(OBJS_SOLVER) \$(OBJS_PSS) \$(OBJS_INPUT) \
\$(OBJS_PAR) \$(OBJS_MATH) \$(OBJS_OUTPUT) \$(OBJS_FORTRAN) \$(OBJS_VISUAL) \$(OBJS_STRUCTURE) \$(OBJS_IO)"

# ALE
if grep '^D_ALE' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_ALE)"
fi

# AXISHELL
if grep '^D_AXISHELL' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_AXISHELL)"
fi

# BEAM3
if grep '^D_BEAM3' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_BEAM3)"
fi

# BRICK1
if grep '^D_BRICK1' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_BRICK1)"
fi

# FLUID2
if grep '^D_FLUID2' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_FLUID2) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID2)"
    fi
fi

# FLUID2_ML
if grep '^FLUID2_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2ML) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2ML)"
    fi
fi

# FLUID3
if grep '^D_FLUID3' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_FLUID3) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID3)"
    fi
fi

# FLUID3_FAST
if grep '^D_FLUID3_F' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3F) \$(OBJS_F3F_F) \$(OBJS_FLUID)"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3F) \$(OBJS_F3F_F)"
    fi
fi

# FLUID3_ML
if grep '^FLUID3_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3ML) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3ML)"
    fi
fi

# FLUID2_PR
if grep '^D_FLUID2_PR' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2PRO) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2PRO)"
    fi
fi

# SHELL8
if grep '^D_SHELL8' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SHELL8)"
fi

# SHELL9
if grep '^D_SHELL9' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SHELL9)"
fi

# WALL1
if grep '^D_WALL1' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_WALL1)"
    OBJ_W1=1
fi

# INTERF
if grep '^D_INTERF' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_IF)"
    OBJ_IF=1
fi

# WALLGE
if grep '^D_WALLGE' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_WALLGE)"
    OBJ_WGE=1
fi

# MLSTRUCT
if grep '^D_MLSTRUCT' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_W1" = "x1" -a "x$OBJ_WGE" = "x1" -a "x$OBJ_IF" = "x1" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML)"
    fi
    if [ "x$OBJ_W1" = "x1" -a "x$OBJ_WGE" = "x1" -a "x$OBJ_IF" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_IF)"
        DEFINES="$DEFINES -DD_INTERF"
    fi
    if [ "x$OBJ_W1" = "x1" -a "x$OBJ_WGE" = "x0" -a "x$OBJ_IF" = "x1" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALLGE)"
        DEFINES="$DEFINES -DD_WALLGE"
    fi
    if [ "x$OBJ_W1" = "x0" -a "x$OBJ_WGE" = "x1" -a "x$OBJ_IF" = "x1" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALL1)"
        DEFINES="$DEFINES -DD_WALL1"
    fi
    if [ "x$OBJ_W1" = "x0" -a "x$OBJ_WGE" = "x0" -a "x$OBJ_IF" = "x1" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALL1) \$(OBJS_WALLGE)"
        DEFINES="$DEFINES -DD_WALL1 -DD_WALLGE"
    fi
    if [ "x$OBJ_W1" = "x0" -a "x$OBJ_WGE" = "x1" -a "x$OBJ_IF" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALL1) \$(OBJS_IF)"
        DEFINES="$DEFINES -DD_WALL1 -DD_INTERF"
    fi
    if [ "x$OBJ_W1" = "x1" -a "x$OBJ_WGE" = "x0" -a "x$OBJ_IF" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALLGE) \$(OBJS_IF)"
        DEFINES="$DEFINES -DD_WALLGE -DD_INTERF"
    fi
    if [ "x$OBJ_W1" = "x0" -a "x$OBJ_WGE" = "x0" -a "x$OBJ_IF" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_S2ML) \$(OBJS_WALL1) \$(OBJS_WALLGE) \$(OBJS_IF)"
        DEFINES="$DEFINES -DD_WALL1 -DD_WALLGE -DD_INTERF"
    fi
fi

# FSI
if grep '^D_FSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_FSI)"
fi

# SSI
if grep '^D_SSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SSI)"
fi

# OPTIM
if grep '^D_OPTIM' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_OPT)"
fi

# MAT
if grep '^D_MAT' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_MAT)"
fi

# NURBS
if grep '^NURBS' "$definefile" 2>&1 > /dev/null ; then
  if [ "x$NURBS_LIB" = "x" ] ; then
    echo $0: Warning: NURBS_PACKAGE not available on this platform.
    echo $0: Warning: Using straight interpolation for subdivision!!
    DEFINES1=`echo $DEFINES | sed -e 's/-DNURBS//'`
    DEFINES=$DEFINES1
  else
    OBJECTS="$OBJECTS \$(OBJS_NURBS_CPP)"
    USES_CPP=1
  fi
fi


