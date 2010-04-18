c
c     Maintainer: Malte Neumann
c                 neumann@statik.uni-stuttgart.de
c                 http://www.uni-stuttgart.de/ibs/members/neumann/
c                 0711 - 685-6121
c
c     ---------------------------------------------------------------
!**********************************************************************
!
!
!   MUMPS VERSION 4.1.6
!   This Version generated on Thu Mar 16 15:39:34 PST 2000
!
!
! COPYRIGHT
! P. R. Amestoy( ENSEEIHT-IRIT/NERSC ), I. S. Duff( CERFACS/RAL ),
! J. Koster(RAL/Parallab), J.-Y. L'Excellent(CERFACS/ENSEEIHT-IRIT),
!             and M. Tuma( CERFACS )
!
!  CERFACS      , Toulouse    (France)  (http://www.cerfacs.fr)
!  ENSEEIHT-IRIT, Toulouse    (France)  (http://www.enseeiht.fr)
!  NERSC-LBL    , Berkeley CA (USA)     (http://www.nersc.gov)
!  PARALLAB     , Bergen      (Norway)  (http://www.parallab.uib.no)
!  RAL          , Oxfordshire (UK)      (http://www.clrc.ac.uk)
!
!
!      THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
!      EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
!
!  The user shall acknowledge the contribution
!  (using references [1] and [2]) of this
!  package in any publication of material dependent upon the use of
!  the package. The user shall use reasonable endeavours to notify
!  the authors of the package of this publication.
!
!  The user can modify this code but, at no time
!  shall the right or title to all or any part of this package pass
!  to the user. The user shall make available free of charge
!  to the authors for any purpose all information relating to any
!  alteration or addition made to this package for the purposes of
!  extending the capabilities or enhancing the performance of this
!  package.
!
!  This package and modified version of this package may
!  not be sold or redistributed.
!
!  One may make copies of the package or modify it
!  for their use provided that the copies, modified or otherwise,
!  are not sold or distributed, and are used under the same
!  terms and conditions.
!
!  [1] Amestoy, Duff and  L'Excellent (1998),
!      Multifrontal parallel distributed symmetric and unsymmetric solvers,
!      to appear in Comput. Methods in Appl. Mech. Eng. (1999).
!
!  [2] Amestoy, Duff, Koster and  L'Excellent (1999),
!      A fully asynchronous multifrontal solver using distributed
!      dynamic scheduling, Technical Report ENSEEIHT-IRIT, RT/APO/99/2.
!      Submitted to SIAM Journal of Matrix Analysis and Appl.
!
!
!  None of the comments from the Copyright notice up to and
!  including this one shall be removed or altered in any way.
!
!**********************************************************************
!
      INCLUDE 'mumps_root.h'
      TYPE MUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite,
!         SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
         INTEGER SYM, PAR
         INTEGER JOB
!
!    ------------------
!    Control parameters
!    ------------------
         INTEGER ICNTL(20)
         DOUBLE PRECISION CNTL(5)
!    --------------------
!    Order of Input matrix
!    --------------------
         INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
         INTEGER NZ
         DOUBLE PRECISION, DIMENSION(:), POINTER :: A
         INTEGER, DIMENSION(:), POINTER :: IRN, JCN
         DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA
!
!        ------------------------------------
!        Case of distributed assembled matrix
!        matrix on entry:
!        ------------------------------------
         INTEGER NZ_loc
         INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
         DOUBLE PRECISION, DIMENSION(:), POINTER :: A_loc
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
         INTEGER NELT
         INTEGER, DIMENSION(:), POINTER :: ELTPTR
         INTEGER, DIMENSION(:), POINTER :: ELTVAR
         DOUBLE PRECISION, DIMENSION(:), POINTER :: A_ELT
!
!    -----------------
!    MPI Communicator
!    -----------------
         INTEGER COMM
!
!    -------------------------------------
!    Ordering, if given by user (optional)
!    -------------------------------------
         INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data
! ******************
!    --------------------------------------------------------
!    RHS
!    ---
!       on input it holds the right hand side
!       on ouput : always hold the assembled solution
!    -------------------------------------------------------
         DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS
!
! **************************
! OUTPUT data and Statistics
! **************************
!
         INTEGER INFO(20)
!        Cost (flops) of subtrees on local process
         DOUBLE PRECISION COST_SUBTREES
         DOUBLE PRECISION RINFO(20)
!
!        Global information -- host only
         DOUBLE PRECISION RINFOG(20)
         INTEGER INFOG(20)
!
!    -------------------------------------
!    Case of distributed matrix on entry:
!    MUMPS potentially provides mapping
!    -------------------------------------
         INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -----------------------------------------
!    Deficiency and null space basis (optional
!    -----------------------------------------
         INTEGER Deficiency
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: NULL_SPACE
!    -----
!    Schur
!    -----
         INTEGER SIZE_SCHUR
         INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
         DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR

!
!
! **********************
! INTERNAL Working data
! **********************
!        For MPI
         INTEGER COMM_NODES, MYID_NODES
         INTEGER  MYID, NPROCS, NSLAVES
         INTEGER ASS_IRECV
         INTEGER, DIMENSION(:), POINTER :: POIDS
         INTEGER LBUFR
         INTEGER LBUFR_BYTES
         INTEGER, DIMENSION(:), POINTER ::  BUFR
!N        Instance number (for access to variables inside modules)
         INTEGER INST_Number
!        for anlaysis/facto/solve phases
         INTEGER MAXIS, MAXS
         INTEGER MAXIS1
         INTEGER KEEP(150)
!N        IS is used for the factors + workspace for contrib. blocks
         INTEGER, DIMENSION(:), POINTER :: IS
!N        is1 (maxis1) contains the arrays computed during analysis and
!N           used for factorization.
         INTEGER, DIMENSION(:), POINTER :: IS1
!N        The only array computed in facto and used by the solve
!N           (except the factors) is PTLUST.
         INTEGER, DIMENSION(:), POINTER :: PTLUST
!        main real working arrays for factorization/solve phases
         DOUBLE PRECISION, DIMENSION(:), POINTER :: S
!        Information on mapping
         INTEGER, DIMENSION(:), POINTER :: PROCNODE
         INTEGER nbsa
!N        Input matrix ready for numerical assembly
!N            -arrowhead format in case of assembled matrix
!N            -element format otherwise
         INTEGER, DIMENSION(:), POINTER :: INTARR
         DOUBLE PRECISION, DIMENSION(:), POINTER :: DBLARR
!N        Element entry: internal data
         INTEGER NELT_LOC, LELTVAR, NA_ELT
         INTEGER, DIMENSION(:), POINTER :: ELTPROC

!   ------------------------
!   Root structure(internal)
!   ------------------------
         TYPE (TYPE_ROOT_STRUC) :: root
      END TYPE MUMPS_STRUC

