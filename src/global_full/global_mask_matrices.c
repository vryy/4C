#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |  calculate the storage mask of the global matrices    m.gee 5/01     |
 |  for various kinds of distributed sparsity patterns                  |
 |  the sparsity pattern implemented at the moment can be found in      |
 |  solution.h                                                          |
 *----------------------------------------------------------------------*/
void mask_global_matrices()
{
int i,j,k,l;               /* some counters */
int isaztec_msr  =0;       /* flag for a certain sparsity pattern */
int ishypre      =0;
int isucchb      =0;
int isdense      =0;
int ismlib_d_sp  =0;       /*mlib direct solver  - sparse */
int isrc_ptr     =0;
int iscolsol     =0;
int isspooles    =0;
FIELD      *actfield;      /* the active field */
PARTITION  *actpart;       /* my partition of the active field */
SOLVAR     *actsolv;       /* the active SOLVAR */
INTRA      *actintra;      /* the field's intra-communicator */
#ifdef DEBUG 
dstrc_enter("mask_global_matrices");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ loop over all fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   actsolv  = &(solv[i]);
   actpart  = &(partition[i]);
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
#else
   /* if we are not parallel here, we have to allocate a pseudo-intracommunicator */
   actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
   if (!actintra) dserror("Allocation of INTRA failed");
   actintra->intra_fieldtyp = actfield->fieldtyp;
   actintra->intra_rank   = 0;
   actintra->intra_nprocs   = 1;
#endif
/*-------------------------------------- not member of this field group */
   if (actintra->intra_fieldtyp==none) continue;
/*--------------------------------------------- first check some values */
   /*----------------------------- check solver and typ of partitioning */
   /*-------- column pointer, row index sparse matrix representation ---*/
   if (actsolv->solvertyp == mlib_d_sp)
   {
      ismlib_d_sp=1;
   }
   /*--------- matrix is distributed modified sparse row DMSR for Aztec */
   if (actsolv->solvertyp==aztec_msr)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with HYPRE"); 
      else isaztec_msr=1;
   }
   /*------------------------------------------- matrix is hypre_parcsr */
   if (
       actsolv->solvertyp==hypre_amg     ||
       actsolv->solvertyp==hypre_pcg     ||
       actsolv->solvertyp==hypre_gmres   ||
       actsolv->solvertyp==hypre_bicgstab 
      )
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with Aztec"); 
      else ishypre=1;
   }
   /*---- matrix is unsym. column compressed Harwell Boeing for superLU */
   if (actsolv->solvertyp==parsuperlu)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with HYPRE"); 
      else isucchb=1;
   }
   /*------------------------ matrix is (non)symmetric dense for Lapack */
   if (actsolv->solvertyp==lapack_nonsym ||
       actsolv->solvertyp==lapack_sym)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with LAPACK"); 
      else isdense=1;
   }
   /*-------------------- matrix is row-column pointer format for Mumps */
   if (actsolv->solvertyp==mumps_sym || actsolv->solvertyp==mumps_nonsym)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with MUMPS"); 
      else isrc_ptr=1;
   }
   /*------------------------------ matrix is skyline format for colsol */
   if (actsolv->solvertyp==colsol_solver)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with Colsol"); 
      else iscolsol=1;
   }
   /*-------------------- matrix is matrix objec for solver lib Spooles */
   if (actsolv->solvertyp==SPOOLES_sym || actsolv->solvertyp==SPOOLES_nonsym)
   {
      if (actsolv->parttyp != cut_elements)
      dserror("Partitioning has to be Cut_Elements for solution with SPOOLES"); 
      else isspooles=1;
   }
   /* allocate only one sparse matrix for each field. The sparsity
      pattern of the matrices for mass and damping and stiffness are  
      supposed to be the same, so they are calculated only once (expensive!) */
   /*-------------------------- for the lower triangle of the matrix ---*/
   if (ismlib_d_sp==1)
   {
      actsolv->nsysarray = 2;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = mds;
         actsolv->sysarray[i].mds = (ML_ARRAY_MDS*)CALLOC(1,sizeof(ML_ARRAY_MDS));
         if (actsolv->sysarray[i].mds==NULL) dserror("Allocation of ML_ARRAY_MDS failed");
      }
      strcpy(actsolv->sysarray[0].mds->arrayname,"gstif1");
      mask_mds(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].mds); 
      ismlib_d_sp=0;
   }
   /*------------------------- matrix is ditributed modified sparse row */
   if (isaztec_msr==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = msr;
         actsolv->sysarray[i].msr = (AZ_ARRAY_MSR*)CALLOC(1,sizeof(AZ_ARRAY_MSR));
         if (actsolv->sysarray[i].msr==NULL) dserror("Allocation of AZ_ARRAY_MSR failed");
         actsolv->sysarray[i].msr->bins=NULL;
      }
      mask_msr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].msr);
      isaztec_msr=0;
   }
   /*------------------------------------- matrix is Spooles's matrix  */
   if (isspooles==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = spoolmatrix;
         actsolv->sysarray[i].spo = (SPOOLMAT*)CALLOC(1,sizeof(SPOOLMAT));
         if (actsolv->sysarray[i].spo==NULL) dserror("Allocation of SPOOLMAT failed");
      }
      mask_spooles(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].spo);
      isspooles=0;
   }
   /*------------------------------------------- matrix is hypre_parcsr */
   if (ishypre==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = parcsr;
         actsolv->sysarray[i].parcsr = (H_PARCSR*)CALLOC(1,sizeof(H_PARCSR));
         if (actsolv->sysarray[i].parcsr==NULL) dserror("Allocation of H_PARCSR failed");
      }
      mask_parcsr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].parcsr);
      ishypre=0;
   }
   /*---------------------------------------------------- matrix is ucchb */
   if (isucchb==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = ucchb;
         actsolv->sysarray[i].ucchb = (UCCHB*)CALLOC(1,sizeof(UCCHB));
         if (actsolv->sysarray[i].ucchb==NULL) dserror("Allocation of UCCHB failed");
      }
      mask_ucchb(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].ucchb);
      isucchb=0;
   }
   /*---------------------------------------------------- matrix is dense */
   if (isdense==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = dense;
         actsolv->sysarray[i].dense = (DENSE*)CALLOC(1,sizeof(DENSE));
         if (actsolv->sysarray[i].dense==NULL) dserror("Allocation of DENSE failed");
      }
      mask_dense(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].dense);
      isdense=0;
   }
   /*----------------------------- matrix is row-column pointer format  */
   if (isrc_ptr==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = rc_ptr;
         actsolv->sysarray[i].rc_ptr = (RC_PTR*)CALLOC(1,sizeof(RC_PTR));
         if (actsolv->sysarray[i].rc_ptr==NULL) dserror("Allocation of RC_PTR failed");
      }
      mask_rc_ptr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].rc_ptr);
      isrc_ptr=0;
   }
   /*---------------------------------------- matrix is skyline format  */
   if (iscolsol==1)
   {
      actsolv->nsysarray = 1;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ || !actsolv->sysarray)
         dserror("Allocation of SPARSE_ARRAY failed");
      for (i=0; i<actsolv->nsysarray; i++)
      {
         actsolv->sysarray_typ[i] = skymatrix;
         actsolv->sysarray[i].sky = (SKYMATRIX*)CALLOC(1,sizeof(SKYMATRIX));
         if (actsolv->sysarray[i].sky==NULL) dserror("Allocation of SKY failed");
      }
      mask_skyline(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].sky);
      iscolsol=0;
   }
/*----------------------------------------------------------------------*/
} /* end of loop over numfld fields */
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_global_matrices */




