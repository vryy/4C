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
      TYPE TYPE_ROOT_STRUC
        SEQUENCE
        INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER MYROW, MYCOL
        INTEGER ROOT_SIZE, TOT_ROOT_SIZE
        INTEGER :: CNTXT_BLACS
        INTEGER, DIMENSION(:), POINTER :: RG2L_ROW
        INTEGER, DIMENSION(:), POINTER :: RG2L_COL
!N used in call to scalapack
        INTEGER , DIMENSION(:), POINTER :: IPIV
        INTEGER, DIMENSION( 9 ) :: DESCRIPTOR, DESCB
        LOGICAL yes
        INTEGER LPIV
!
!      Data for nullspace/QR
!
        DOUBLE PRECISION, DIMENSION(:), POINTER :: QR_TAU
!
!      Givens rotations
!
        INTEGER MAXG, GIND
        DOUBLE PRECISION, DIMENSION(:),POINTER::GROW, GCOS, GSIN
!
!      RRRLU data
!
        INTEGER ELG_MAX,NULL_MAX
        INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
        INTEGER,DIMENSION(:),POINTER::PERM_ROW,PERM_COL
        INTEGER,DIMENSION(:),POINTER::ELROW, EUROW, PTREL, PTREU
        DOUBLE PRECISION, DIMENSION(:), POINTER :: ELELG, EUELG, DL
!
      END TYPE TYPE_ROOT_STRUC

