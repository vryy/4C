
# fluid objects:
OBJ_FLUID=0


# objects:
OBJECTS="\$(OBJS_MAIN) \$(OBJS_GLOBAL) \$(OBJS_SOLVER) \$(OBJS_PSS) \$(OBJS_INPUT) \
\$(OBJS_PAR) \$(OBJS_MATH) \$(OBJS_OUTPUT) \$(OBJS_FORTRAN) \$(OBJS_VISUAL) \$(OBJS_STRUCTURE)"

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
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID2)"
    fi
fi

# FLUID2_ML
if grep '^FLUID2_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2ML) \$(OBJS_FLUID)"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F2ML)"
    fi
fi

# FLUID3
if grep '^D_FLUID3' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_FLUID3) \$(OBJS_FLUID)"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_FLUID3)"
    fi
fi

# FLUID3_ML
if grep '^FLUID3_ML' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F3ML) \$(OBJS_FLUID)"
        OBJ_FLUID=1
    else
        OBJECTS="$OBJECTS \$(OBJS_F3ML)"
    fi
fi

# FLUID2_PR
if grep '^D_FLUID2_PR' "$definefile" 2>&1 > /dev/null ; then
    if [ "x$OBJ_FLUID" = "x0" ] ; then
        OBJECTS="$OBJECTS \$(OBJS_F2PRO) \$(OBJS_FLUID)"
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
fi

# INTERF
if grep '^D_INTERF' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_IF)"
fi

# WALLGE
if grep '^D_WALLGE' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_WALLGE)"
fi

# FSI
if grep '^D_FSI' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_FSI)"
fi

# OPTIM 
if grep '^D_OPTIM' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_OPT)"
fi

# MAT
if grep '^D_MAT' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_MAT)"
fi

# LS
if grep '^D_LS' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_LS)"
fi

# XFEM
if grep '^D_XFEM' "$definefile" 2>&1 > /dev/null ; then
    OBJECTS="$OBJECTS \$(OBJS_XFEM)"
fi
