#--------------------------------------------------------------------
#    Makefile CCARAT (parallel and sequentiel)
#--------------------------------------------------------------------
#
#--------------------------------------------------------------------
# list of possible defines:
#--------------------------------------------------------------------
#  
#  ELEMENTS:
#  =========
# -DD_ALE        defines element ale (2D+3D)
# -DD_AXISHELL   defines axishell element
# -DD_BEAM3      defines element beam3
# -DD_BRICK1     defines element brick1
# -DD_FLUID2     defines element fluid2
# -DD_FLUID2_PRO defines element fluid2 with projection method
# -DFLUID2_ML    defines element fluid2 with multilevel 
# -DD_FLUID3     defines element fluid3
# -DFLUID3_ML    defines element fluid3 with multilevel 
# -DD_SHELL8     defines element shell8
# -DD_SHELL9     defines element shell9
# -DD_WALL1      defines element wall1
# -DD_INTERF     defines interface element
# -DD_WALLGE     defines gradient enhanced wall element
# -DD_XFEM       defines for enriched fluid element
#
#  PROBLEMTYPES:
#  =============
# -DD_FLUID      defines general fluid problem
# -DD_FSI        defines general FSI problem
# -DD_OPTIM      defines optimization
# -DD_CONTACT    defines contact algorithms
# -DD_LS         defines for level set problem
#
#  ADDITIONAL FEATURES:
#  ====================
# -DD_MAT        defines materials
# -DS8CONTACT    defines contact algorithms for shell8 contact
# -DWALLCONTACT  defines contact algorithms for wall contact
# -DGEMM         defines Generalized Energy-Momentum Int. Scheme
#
#  TOOLS:
#  ======
# -DRESULTTEST   defines defines result testing before checking in a new version
# -DPERF         defines routines for exact time measurements 
#
#  PLATFORMS:
#  ==========
# -DHPUX10       defines for compilation on HP-UX 10.20
# -DHPUX11       defines for compilation on HP-UX 11.0
# -DHPUXITA      defines compilation for Itanium2 running HP-UX 11.23
# -DSUSE73       defines for compilation and linkage with gcc and g77 on Linux
# -DAZUSA        defines for compilation on azusa.hww.de
# -DSUN          defines for compilation and linkage on CUSS Sunfire machine
#
# -DDEBUG        defines debug mode: routine and array tracing 
# -DPARALLEL     defines parallel version of programm
# -DSIXTYFOUR    defines 64bit compilation and linkage
#
#  SOLVERS:
#  ========
# -DAZTEC_PACKAGE      defines interface to solver library AZTEC
# -DUMFPACK            defines interface to solver library Umfpack
# -DSPOOLES_PACKAGE    defines interface to solver library Spooles (only in parallel HPUX11)
# -DTRILINOS_PACKAGE   defines use of Sandia's Trilinos library (currently only in parallel under Mac OSX)
#
#  INACTIVE SOLVERS:
#  =================
# -DMLIB_PACKAGE       defines interface to solver library MLIB (HPUX11, sequentiell)
# -DHYPRE_PACKAGE      defines interface to solver library HYPRE (not yet impl.)
# -DPARSUPERLU_PACKAGE defines interface to solver library ParSuperLU
# -DSUPERLU_PACKAGE    defines interface to solver library SuperLU 
# -DMUMPS_PACKAGE      defines interface to solver library Mumps (only in parallel HPUX11)
#
#  VISUAL:
#  =======
# -DVISUAL2_PACKAGE    defines visualisation tool VISUAL2 (only in seriell)
# -DVISUAL3_PACKAGE    defines visualisation tool VISUAL3 (only in seriell)
#
#
#--------------------------------------------------------------------
#  use these defines for compilation
#--------------------------------------------------------------------
# parallel:
#serial:
# -------------------------------------------------------------------
#
#--------------------------------------------------------------------
#
# SERIAL PART 32bit
#
#--------------------------------------------------------------------
PFL1  = -DMACOSX -DD_BRICK1 -DD_MAT \
        -DD_SHELL8 -DAZTEC_PACKAGE -DTRILINOS_PACKAGE
FXX           = g77
CCX           = gcc
CXX           = g++
LDX           = g++
PFL           = $(PFL1)
SPOOLES       = 
TRILINOS_ARCH = LINUX_SERIAL
# --------------------------------------------------------------------
CC	 = $(CCX)
F77      = $(FXX)
LD       = $(LDX)
RM	 = /bin/rm -f
# -------------------------------------------------------------------
# this is debug
# -------------------------------------------------------------------
PROGRAM  = resulter.exe 
CFLAGS   = -v -g -ansi -Wall $(PFL)  -DDEBUG 
CXXFLAGS = -g -fno-inline $(PFL) $(TRILINOS_DEFINES) -DHAVE_CONFIG_H -DDEBUG 
FFLAGS   = -g -I/usr/include/
# -------------------------------------------------------------------
# this is fast 
# -------------------------------------------------------------------
#PROGRAM  = resulter.exe  
#CFLAGS   = -v -O3 -mcpu=pentium4 -march=pentium4 -ansi $(PFL) 
#CXXFLAGS = -v -o -ansi -Wall $(PFL) $(TRILINOS_DEFINES) 
#FFLAGS   = -O3 -I/usr/include
# -------------------------------------------------------------------
LDFLAGS	 =  

INCLUDES = -I/usr/include $(TRILINOS_INCLUDES) \
           -I/usr/local/sw/metis-4.0/Lib \
           -I/usr/local/sw/spooles
#
LDOPTS   = -a archive -v -W1 -M
#
#
LIBDIR = /usr/local/sw
#
#
LIBS	 =  -L$(LIBDIR)/spooles $(SPOOLES) $(TRILINOS_LIBRARIES) \
            -L/usr/lib \
            -L/usr/local/sw/lam/lib \
            -L$(LIBDIR)/metis-4.0 -lmetis \
	    -llapack -lblas \
            -L/usr/local/lib -L/usr/lib  \
            -lc -lm -lg2c 
#               
#
#
#----------------------- binaries -----------------------------------
OBJS_MAIN = main.o \
            read_mesh.o \
            printf_info.o \
            take_every.o \
            merge.o \
            pss_am.o
#--------------------------------------------------------------------
#
#
SRC=.
#
#
#
$(PROGRAM):	$(OBJS_MAIN) 
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) \
		$(OBJS_MAIN) \
		$(LIBS) $(INCLUDES) -o $(PROGRAM)
		@echo "done"



.SUFFIXES : .o .cc .cpp .c .f

.cc.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.f.o:
	$(F77) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf 	$(OBJS_MAIN) 
#
