#--------------------------------------------------------------------
# switch between sequentiell and parallel compilation by uncommenting 
# the approbiate compilers:
#--------------------------------------------------------------------
#
# sequentiell part
#--------------------------------------------------------------------
PROGRAM_NAME  = resulter.exe
FXX     = f77
CCX     = cc
LDX     = cc
PFL     = 
AZTEC   = -laztec_seriell21
UMFPACK = -lumfpack40_hpux10_32
##--------------------------------------------------------------------
#
# parallel part (not implemented!!!)
#--------------------------------------------------------------------
#PROGRAM_NAME  = resulter.exe
#FXX     = /opt/mpi/bin/mpif77
#CCX     = /opt/mpi/bin/mpicc
#LDX     = /opt/mpi/bin/mpicc
#PFL     = -DPARALLEL  -DD_FLUID2  -DD_SHELL8 -DD_FLUID3 -DD_ALE -DD_BRICK1 -DD_WALL1 \
#          -DAZTEC_PACKAGE -DUMFPACK
#AZTEC   = -laztec_par21
#UMFPACK = -lumfpack40_hpux10_32
#
#--------------------------------------------------------------------
#
CC	 = $(CCX)
LD       = $(CCX)
RM	 = /bin/rm -f

#--------------------------------------------------------------------
# switch between fast and debug compilation by uncommenting 
# the approbiate flags:
#------------------------------------------------ this is debug-32bit
CFLAGS        = -D__HP__ +w1 +e -Aa -Ae  -g ${PFL} 
FFLAGS        = -g -I/opt/mpi/include/
#------------------------------------------------- this is fast-32bit
#CFLAGS       = -D__HP__  +e -Ae +w1  -OAa +Onolimit ${PFL}
#FFLAGS        = +O2 +Onolimit -I/opt/mpi/include/
#--------------------------------------------------------------------
LDFLAGS	 = -L/usr/lib -L/opt/fortran/lib -lcl -lc -lM -lm 
INCLUDES = -I/opt/mpi/include/ \
           -I/bau/stat33/users/statik/lib/hypre-1.6.0/src/hypre/include/ \
           -I/bau/stat33/users/statik/lib/SULUDIST0210/SuperLU_DIST/SRC \
           -I/bau/stat33/users/statik/lib/METIS/ \
           -L/bau/stat33/users/statik/lib/METIS/ -lmetis \
           -L/bau/stat33/users/statik/lib/AZTEC21_MPI/ \
           -L/bau/stat33/users/statik/lib/AZTEC21/ $(AZTEC) \
           -L/bau/stat33/users/statik/lib \
           -L/bau/stat33/users/statik/lib/hypre-1.6.0/src/hypre/lib \
           -L/bau/stat33/users/statik/lib/SULUDIST0210/SuperLU_DIST \
           -L/bau/stat33/users/statik/lib/SCALAPACK \
           -L/bau/stat33/users/statik/lib/BLACS/LIB \
           -L/bau/stat33/users/statik/lib/MUMPS4.1.6/lib \
           -L/bau/stat33/users/statik/lib/UMFPACK/lib \
           -llapack_hppa \
           -lHYPRE_DistributedMatrix \
           -lHYPRE_DistributedMatrixPilutSolver \
           -lHYPRE_Euclid \
           -lHYPRE_FEI \
           -lHYPRE_IJ_mv \
           -lHYPRE_LSI \
           -lHYPRE_MatrixMatrix \
           -lHYPRE_ParaSails \
           -lHYPRE_blas \
           -lHYPRE_parcsr_ls \
           -lHYPRE_parcsr_mv \
           -lHYPRE_seq_mv \
           -lHYPRE_sstruct_ls \
           -lHYPRE_sstruct_mv \
           -lHYPRE_struct_ls \
           -lHYPRE_struct_mv \
           -lHYPRE_utilities \
           -lkrylov \
           -lsuperlu_hppa \
           -lblacs_MPI-HPPA-0 \
           -lblacsCinit_MPI-HPPA-0 \
           -lblacsF77init_MPI-HPPA-0 \
           -lscalapack_HPPA \
           -lmumps \
           -lblas \
           -lumfpack40_hpux10_32

# sparsekit is not used at the moment
#           -lskit 
# libHYPRE_superlu.a convlicts with libsuperlu_hppa.a
#           -lHYPRE_superlu

LDOPTS        = -a archive -v


LIBS	      =  
#
# binaries
#                main
OBJS1         = main.o \
                read_mesh.o \
                printf_info.o \
                take_every.o \
                merge.o \
                pss_am.o
#                fortran
#OBJS13	      = ${SRC}/fortran/colsol.o
#
#
# targets
#
PROGRAM       = $(PROGRAM_NAME)
#
# dependencies
#
$(PROGRAM):     $(OBJS1)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS1) $(LIBS) $(INCLUDES) -o $(PROGRAM)
		@echo "done"

#
$(OBJS1):     
		@echo "Making $@ ..."
		cd ${@D}; $(CC) $(CFLAGS) $(INCLUDES) -c  $<
		@echo "MainMake is back."
# fortran
#$(OBJS13):     
#		@echo "Making $@ ..."
#		cd ${@D}; $(FXX) $(FFLAGS)  -c  $<
#		@echo "MainMake is back."
##
.c.o:
	@ echo "$@ wird kompiliert...."
	$(CC) $(CFLAGS) $(INCLUDES) -c  $<

.f.o:
	@ echo "$@ wird kompiliert...."
	$(FXX) $(FFLAGS) -c  $<
clean:
	rm -rf  $(OBJS1)  
$
