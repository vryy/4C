
# fluid objects:
OBJ_FLUID=0

# control-flag, to ensure that OBJS_S2ML get OBJS_WALL,OBJS_WALLGE,OBJS_IF:
OBJ_W1=0
OBJ_WGE=0
OBJ_IF=0


# objects:
OBJECTS="\$(OBJS_MAIN) \$(OBJS_GLOBAL) \$(OBJS_GLOBAL_CPP) \$(OBJS_SOLVER) \$(OBJS_SOLVER_CPP) \
\$(OBJS_DRT) \$(OBJS_DRT_LIB) \$(OBJS_PSS) \$(OBJS_INPUT) \$(OBJS_INPUT_CPP) \$(OBJS_PAR) \$(OBJS_MATH) \
\$(OBJS_OUTPUT) \$(OBJS_FORTRAN) \$(OBJS_VISUAL) \$(OBJS_STRUCTURE) \$(OBJS_IO) \$(OBJS_IO_LIB) \
\$(OBJS_DRT_FSI) \$(OBJS_DRT_FLUID) \$(OBJS_DRT_MAT_LIB) \$(OBJS_DRT_MAT) \$(OBJS_DRT_STRU_MULTI)"

FILTER_OBJECTS="\$(OBJS_POST_DRT_COMMON) \
\$(OBJS_IO_LIB) \$(OBJS_DRT_LIB) \$(OBJS_DRT_MAT_LIB) \$(OBJS_PSS) \$(OBJS_PAR) \
\$(OBJS_DRT_F2_LIB) \$(OBJS_DRT_F3_LIB)  \$(OBJS_DRT_XF3_LIB) \$(OBJS_DRT_S8_LIB) \
\$(OBJS_SOH8_LIB) \$(OBJS_DRT_W1_LIB) \$(OBJS_DRT_ALE_LIB)"

# ALE
if grep '^[[:blank:]]*D_ALE' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_ALE) \$(OBJS_DRT_ALE_LIB) \$(OBJS_DRT_ALE)"
fi

# AXISHELL
if grep '^[[:blank:]]*D_AXISHELL' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_AXISHELL)"
fi

# BEAM3
if grep '^[[:blank:]]*D_BEAM3' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_BEAM3)"
fi

# drt solid hex 8
if grep '^[[:blank:]]*D_SOH8' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SOH8) \$(OBJS_SOH8_LIB)"
fi

# BRICK1
if grep '^[[:blank:]]*D_BRICK1' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_BRICK1)"
fi

# FLUID2
if grep '^[[:blank:]]*D_FLUID2' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_FLUID2) \$(OBJS_FLUID) \$(OBJS_DRT_F2) \$(OBJS_DRT_F2_LIB)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID2) \$(OBJS_DRT_F2) \$(OBJS_DRT_F2_LIB)"
    fi
fi

# FLUID2_ML
if grep '^[[:blank:]]*FLUID2_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2ML) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2ML)"
    fi
fi

# FLUID2_TDS
if grep '^[[:blank:]]*D_FLUID2_TDS' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS  \$(OBJS_F2_TDS)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2_TDS)"
    fi
fi

# FLUID3
if grep '^[[:blank:]]*D_FLUID3' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_FLUID3) \$(OBJS_FLUID) \$(OBJS_DRT_F3_LIB) \$(OBJS_DRT_F3)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID3) \$(OBJS_DRT_F3_LIB) \$(OBJS_DRT_F3)"
    fi
fi

# FLUID3_FAST
if grep '^[[:blank:]]*D_FLUID3_F' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3F) \$(OBJS_F3F_F) \$(OBJS_FLUID)"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3F) \$(OBJS_F3F_F)"
    fi
fi

# FLUID3_ML
if grep '^[[:blank:]]*FLUID3_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3ML) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3ML)"
    fi
fi

# FLUID2_PRO
if grep '^[[:blank:]]*D_FLUID2_PRO' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2PRO) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2PRO)"
    fi
fi

# FLUID3_PRO
if grep '^[[:blank:]]*D_FLUID3_PRO' "$definefile" 3>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3PRO) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3PRO)"
    fi
fi

# FLUID2_IS
if grep '^[[:blank:]]*D_FLUID2_IS' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2IS) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2IS)"
    fi
fi

# FLUID3_IS
if grep '^[[:blank:]]*D_FLUID3_IS' "$definefile" 3>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3IS) \$(OBJS_FLUID)"
        DEFINES="$DEFINES -DD_FLUID"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3IS)"
    fi
fi

# SHELL8
if grep '^[[:blank:]]*D_SHELL8' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SHELL8) \$(OBJS_DRT_S8) \$(OBJS_DRT_S8_LIB)"
fi

# SHELL9
if grep '^[[:blank:]]*D_SHELL9' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SHELL9)"
fi

# WALL1
if grep '^[[:blank:]]*D_WALL1' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_WALL1) \$(OBJS_DRT_W1_LIB) \$(OBJS_DRT_W1)"
    OBJ_W1=1
fi

# INTERF
if grep '^[[:blank:]]*D_INTERF' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_IF)"
    OBJ_IF=1
fi

# WALLGE
if grep '^[[:blank:]]*D_WALLGE' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_WALLGE)"
    OBJ_WGE=1
fi

# MLSTRUCT
if grep '^[[:blank:]]*D_MLSTRUCT' "$definefile" 2>&1 > /dev/null ; then
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

# THERM2
if grep '^[[:blank:]]*D_THERM2' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_THERM2)"
fi

# THERM3
if grep '^[[:blank:]]*D_THERM3' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_THERM3)"
fi

# SOLID3
if grep '^[[:blank:]]*D_SOLID3' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SOLID3)"
fi

# FSI
if grep '^[[:blank:]]*D_FSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_FSI) \$(OBJS_FSI_CPP)"
fi

# SSI
if grep '^[[:blank:]]*D_SSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_SSI)"
fi

# OPTIM
if grep '^[[:blank:]]*D_OPTIM' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_OPT)"
fi

# TSI
if grep '^[[:blank:]]*D_TSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_TSI)"
fi

# MAT
if grep '^[[:blank:]]*D_MAT' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_MAT)"
fi

# NURBS
if grep '^[[:blank:]]*NURBS' "$definefile" 2>&1 > /dev/null ; then
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


