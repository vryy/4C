/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#ifdef TRILINOS_PACKAGE
#include "solver_trilinos_service.H"
#endif

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
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld submeshFIELDs, defined in global_control.c          |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
/*----------------------------------------------------------------------*
 | global variable *sm_solv, vector of lenght numfld of structures SOLVAR  |
 | defined in input_submesh.c                                          |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *sm_solv;
/*!----------------------------------------------------------------------
proc's info about his partition, global variable defined in input_submesh.c
*----------------------------------------------------------------------*/
extern struct _PARTITION  *sm_part;
#endif /* D_MLSTRUCT */

/*----------------------------------------------------------------------*
 |  calculate the storage mask of the global matrices    m.gee 5/01     |
 |  for various kinds of distributed sparsity patterns                  |
 |  the sparsity pattern implemented at the moment can be found in      |
 |  solution.h                                                          |
 |  for several discretisations:                                        |
 |    - at the moment only implemented for AZTEC_msr                    |
 |    - all discretisations have to use the same solver with the same   |
 |      parameter                                                       |
 |    - evertying is done seperate for each discritisations             |
 *----------------------------------------------------------------------*/
void mask_global_matrices()
{


  INT i,j;               /* some counters */

#ifdef AZTEC_PACKAGE
  INT isaztec_msr  =0;       /* flag for a certain sparsity pattern */
#endif

#ifdef HYPRE
  INT ishypre      =0;
#endif

#ifdef PARSUPERLU_PACKAGE
  INT isucchb      =0;
#endif

  INT isdense      =0;

#ifdef MLIB_PACKAGE
  INT ismlib_d_sp  =0;       /*mlib direct solver  - sparse */
#endif

#ifdef MUMPS_PACKAGE
  INT isrc_ptr     =0;
#endif

  INT iscolsol     =0;

#ifdef SPOOLES_PACKAGE
  INT isspooles    =0;
#endif

#ifdef UMFPACK
  INT isumfpack    =0;
#endif

#ifdef MLPCG
  INT ismlpcg      =0;
#endif

  /*
     INT nsysarray    =1;
     INT disnum       =0;
     */

  INT numeq;
  INT numeq_total;

  FIELD      *actfield;        /* the active field */
  PARTITION  *actpart;         /* my partition of the active field */
  SOLVAR     *actsolv;         /* the active SOLVAR */
  INTRA      *actintra = NULL; /* the field's intra-communicator */

  /* use trilinos for internal vectors/matrices */
  INT  usetrilinosalgebra = genprob.usetrilinosalgebra;

#ifdef DEBUG
  dstrc_enter("mask_global_matrices");
#endif


  /* loop over all fields */
  /*======================*/
  for (j=0; j<genprob.numfld; j++)
  {

    actfield = &(field[j]);
    actsolv  = &(solv[j]);
    actpart  = &(partition[j]);

#ifdef PARALLEL
    actintra = &(par.intra[j]);
#else
    /* if we are not parallel here, we have to allocate a pseudo-intracommunicator */
    actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
    actintra->intra_fieldtyp = actfield->fieldtyp;
    actintra->intra_rank   = 0;
    actintra->intra_nprocs   = 1;
#endif


    /* not member of this field group */
    if (actintra->intra_fieldtyp==none)
      continue;


    /* determine number of sysarrays in this field */
    actsolv->nsysarray = actfield->ndis;


    /* special case:  matrix typ is trilinos for ALL fields */
    /* (to be precise: Epetra_CrsMatrix)                    */
    /*------------------------------------------------------*/
    if (usetrilinosalgebra)
    { 
      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = trilinos;
        actsolv->sysarray[i].trilinos = (TRILINOSMATRIX*)CCACALLOC(1,sizeof(TRILINOSMATRIX));
#ifdef TRILINOS_PACKAGE
        init_trilinos_matrix(actfield,actpart,actsolv,actintra,actsolv->sysarray[i].trilinos,i);
#else
        dserror("TRILINOS_PACKAGE not defined");
#endif
      }  /* for (i=0; i<actfield->ndis; i++) */
      continue;
    }

    /* special case:  matrix typ is OLL */
    /*----------------------------------*/
    if (actsolv->matrixtyp == oll_matrix)
    {
      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actfield->ndis; i++)
      {
        actsolv->sysarray_typ[i] = trilinos;
        actsolv->sysarray[i].oll = (OLL*)CCACALLOC(1,sizeof(OLL));

        numeq_total = actfield->dis[i].numeq;
        oll_numeq(actfield, actpart, actintra, i, &numeq);

        oll_open(actsolv->sysarray[i].oll, numeq, numeq_total,
            actfield, actpart, actintra, i);

      }  /* for (i=0; i<actfield->ndis; i++) */

      continue;
    }


    /* first check some values */
    /*-------------------------*/


#ifdef MLIB_PACKAGE
    if (actsolv->solvertyp == mlib_d_sp)
    {
      ismlib_d_sp=1;
    }
#endif


    /* matrix is distributed modified sparse row DMSR for Aztec */
#ifdef AZTEC_PACKAGE
    if (actsolv->solvertyp==aztec_msr)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with Aztec");
      else isaztec_msr=1;
    }
#endif


    /* matrix is hypre_parcsr */
#ifdef HYPRE_PACKAGE
    if (
        actsolv->solvertyp==hypre_amg     ||
        actsolv->solvertyp==hypre_pcg     ||
        actsolv->solvertyp==hypre_gmres   ||
        actsolv->solvertyp==hypre_bicgstab
       )
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with Hypre");
      else ishypre=1;
    }
#endif


    /* matrix is unsym. column compressed Harwell Boeing for superLU */
#ifdef PARSUPERLU_PACKAGE
    if (actsolv->solvertyp==parsuperlu)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with Superlu");
      else isucchb=1;
    }
#endif


    /* matrix is (non)symmetric dense for Lapack */
    if (actsolv->solvertyp==lapack_nonsym ||
        actsolv->solvertyp==lapack_sym)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with LAPACK");
      else isdense=1;
    }


    /* matrix is row-column pointer format for Mumps */
#ifdef MUMPS_PACKAGE
    if (actsolv->solvertyp==mumps_sym || actsolv->solvertyp==mumps_nonsym)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with MUMPS");
      else isrc_ptr=1;
    }
#endif


    /* matrix is compressed column format for Umfpack */
#ifdef UMFPACK
    if (actsolv->solvertyp==umfpack)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with Umfpack");
      else isumfpack=1;
    }
#endif


    /* matrix is skyline format for colsol */
    if (actsolv->solvertyp==colsol_solver)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with Colsol");
      else iscolsol=1;
    }


    /* matrix is matrix objec for solver lib Spooles */
#ifdef SPOOLES_PACKAGE
    if (actsolv->solvertyp==SPOOLES_sym || actsolv->solvertyp==SPOOLES_nonsym)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with SPOOLES");
      else isspooles=1;
    }
#endif


    /* matrix is block distributed csr format for mlpcg */
#ifdef MLPCG
    if (actsolv->solvertyp==mlpcg)
    {
      if (actsolv->parttyp != cut_elements)
        dserror("Partitioning has to be Cut_Elements for solution with MLPCG");
      else ismlpcg=1;
    }
#endif


#ifdef MLIB_PACKAGE
    /* allocate only one sparse matrix for each field. The sparsity
       pattern of the matrices for mass and damping and stiffness are
       supposed to be the same, so they are calculated only once (expensive!)
       for the lower triangle of the matrix */

    if (ismlib_d_sp==1)
    {
      actsolv->nsysarray = 2;
      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));
      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = mds;
        actsolv->sysarray[i].mds = (ML_ARRAY_MDS*)CCACALLOC(1,sizeof(ML_ARRAY_MDS));
      }
      strcpy(actsolv->sysarray[0].mds->arrayname,"gstif1");
      mask_mds(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].mds);
      ismlib_d_sp=0;
    }
#endif



    /* matrix is distributed modified sparse row */
#ifdef AZTEC_PACKAGE
    if (isaztec_msr==1)
    {

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = msr;
        actsolv->sysarray[i].msr = (AZ_ARRAY_MSR*)CCACALLOC(1,sizeof(AZ_ARRAY_MSR));
        actsolv->sysarray[i].msr->bins=NULL;
        mask_msr(actfield,actpart,actsolv,actintra,actsolv->sysarray[i].msr,i);
      }
      isaztec_msr=0;
    }
#endif



    /* matrix is Spooles's matrix  */
#ifdef SPOOLES_PACKAGE
    if (isspooles==1)
    {

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = spoolmatrix;
        actsolv->sysarray[i].spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
        mask_spooles(actfield,actpart,actsolv,actintra,actsolv->sysarray[i].spo,i);
      }

      isspooles=0;
    }
#endif



    /* matrix is hypre_parcsr */
#ifdef HYPRE_PACKAGE
    if (ishypre==1)
    {
      if(actsolv->nsysarray>1)
        dserror("different discretisations not possible with SOLVER_TYP 'HYPRE'\n");

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = parcsr;
        actsolv->sysarray[i].parcsr = (H_PARCSR*)CCACALLOC(1,sizeof(H_PARCSR));
      }

      mask_parcsr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].parcsr);
      ishypre=0;
    }
#endif



    /* matrix is ucchb */
#ifdef PARSUPERLU_PACKAGE
    if (isucchb==1)
    {
      if(actsolv->nsysarray>1)
        dserror("different discretisations not possible with SOLVER_TYP 'superlu'\n");

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = ucchb;
        actsolv->sysarray[i].ucchb = (UCCHB*)CCACALLOC(1,sizeof(UCCHB));
      }

      mask_ucchb(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].ucchb);
      isucchb=0;
    }
#endif




    /* matrix is dense */
    if (isdense==1)
    {
      if(actsolv->nsysarray>1)
        dserror("different discretisations not possible with SOLVER_TYP 'Lapack'\n");

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = dense;
        actsolv->sysarray[i].dense = (DENSE*)CCACALLOC(1,sizeof(DENSE));
      }

      mask_dense(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].dense);
      isdense=0;
    }




    /* matrix is row-column pointer format  */
#ifdef MUMPS_PACKAGE
    if (isrc_ptr==1)
    {
      if(actsolv->nsysarray>1)
        dserror("different discretisations not possible with SOLVER_TYP 'Mumps'\n");

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = rc_ptr;
        actsolv->sysarray[i].rc_ptr = (RC_PTR*)CCACALLOC(1,sizeof(RC_PTR));
      }

      mask_rc_ptr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].rc_ptr);
      isrc_ptr=0;
    }
#endif




    /* matrix is row-column pointer format for umfpack solver */
#ifdef UMFPACK
    if (isumfpack==1)
    {

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = ccf;
        actsolv->sysarray[i].ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
        mask_ccf(actfield,actpart,actsolv,actintra,actsolv->sysarray[i].ccf,i);
      }

      isumfpack=0;
    }
#endif




    /* matrix is skyline format  */
    if (iscolsol==1)
    {

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = skymatrix;
        actsolv->sysarray[i].sky = (SKYMATRIX*)CCACALLOC(1,sizeof(SKYMATRIX));
        mask_skyline(actfield,actpart,actsolv,actintra,actsolv->sysarray[i].sky,i);
      }

      iscolsol=0;
    }




#ifdef MLPCG
    if (ismlpcg==1)
    {
      if(actsolv->nsysarray>1)
        dserror("different discretisations not possible with SOLVER_TYP 'MLPCG'\n");

      actsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_TYP));
      actsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsolv->nsysarray,sizeof(SPARSE_ARRAY));

      for (i=0; i<actsolv->nsysarray; i++)
      {
        actsolv->sysarray_typ[i] = bdcsr;
        actsolv->sysarray[i].bdcsr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
      }

      mask_bdcsr(actfield,actpart,actsolv,actintra,actsolv->sysarray[0].bdcsr);
      actsolv->mlpcgvars->fielddis = &(actfield->dis[0]);
      actsolv->mlpcgvars->partdis  = &(actpart->pdis[0]);
      ismlpcg=0;
    }
#endif




  }  /* for (j=0; j<genprob.numfld; j++) */


#ifndef PARALLEL
  CCAFREE(actintra);
#endif


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
} /* end of mask_global_matrices */








#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  calculate the storage mask of the submesh matrix k''    ah 03/04    |
 |  for (various kinds of distributed sparsity patterns) CCF->UMPFPACK  |
 |  and initialize the solver and the assemblation                      |
 *----------------------------------------------------------------------*/
void mask_submesh_matrices()
{
INT i;               /* some counters */
INT isumfpack    =0;

INT actsmsysarray=0;

INT numeq;
INT numeq_total;
INT init;                    /* flag for solver_control call */

FIELD      *actsmfield;      /* the active submeshfield */
PARTITION  *actsmpart;       /* my partition of the active submeshfield */
SOLVAR     *actsmsolv;       /* the active submeshSOLVAR */
INTRA      *actsmintra;      /* the field's submesh-intra-communicator */

#ifdef DEBUG
dstrc_enter("mask_submesh_matrices");
#endif
/*----------------------------------------------------------------------*/
/*              mask for submesh stiffness k prime prime                */
/*----------------------------------------------------------------------*/
actsmfield = &(sm_field[0]);
actsmsolv  = &(sm_solv[0]);
actsmpart  = &(sm_part[0]);
/* because we are not parallel here, we have to allocate a pseudo-intracommunicator */
actsmintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
actsmintra->intra_fieldtyp = structure;
actsmintra->intra_rank   = 0;
actsmintra->intra_nprocs   = 1;
/*-------------------- matrix is compressed column format for Umfpack */
if (actsmsolv->solvertyp==umfpack)
{
  if (actsmsolv->parttyp != cut_elements)
    dserror("Partitioning has to be Cut_Elements for solution with Umfpack");
  else isumfpack=1;
}
else dserror("only umpfpack for submesh solver ");
/*---------------- matrix is row-column pointer format for umfpack solver */
if (isumfpack==1)
{
  actsmsolv->nsysarray = 1;
  actsmsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(actsmsolv->nsysarray,sizeof(SPARSE_TYP));
  actsmsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(actsmsolv->nsysarray,sizeof(SPARSE_ARRAY));
  actsmsolv->sysarray_typ[0] = ccf;
  actsmsolv->sysarray[0].ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
  /*--------------------------------------------------------------------*/
  mask_ccf(actsmfield,actsmpart,actsmsolv,actsmintra,actsmsolv->sysarray[0].ccf,0);
  isumfpack=0;
}
/*----------------------------------------------------------------------*/
/*            init solver for sm-equation and assemblation of k ''      */
/*----------------------------------------------------------------------*/

numeq_total = actsmfield->dis[0].numeq;
numeq = numeq_total;

/*-------------------------- allocate 10 dist. load vectors for the RHS */
actsmsolv->nrhs = 10;
solserv_create_vec(&(actsmsolv->rhs),actsmsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsmsolv->nrhs; i++) solserv_zero_vec(&(actsmsolv->rhs[i]));
/*--------------------- allocate 10 dist. load vectors for the solution */
actsmsolv->nsol= 10;
solserv_create_vec(&(actsmsolv->sol),actsmsolv->nsol,numeq_total,numeq,"DV");
for (i=0; i<actsmsolv->nsol; i++) solserv_zero_vec(&(actsmsolv->sol[i]));
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsmsolv, actsmintra,
             &(actsmsolv->sysarray_typ[0]),
             &(actsmsolv->sysarray[0]),
             &(actsmsolv->sol[0]),
             &(actsmsolv->rhs[0]),
              init);
/*--------------------- init the Assemblation of submesh stiffness k'' */
/*----------------------------- init the assembly for ONE sparse matrix */
init_assembly(actsmpart,actsmsolv,actsmintra,actsmfield,actsmsysarray,0);


CCAFREE(actsmintra);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}
#endif /* D_MLSTRUCT */

