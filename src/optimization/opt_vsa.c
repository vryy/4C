/*!----------------------------------------------------------------------
\file
\brief contains the routine 'optvsa', 
       variational sensitivity analysis

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../headers/optimization.h"
/*! 
\addtogroup OPTIMIZATION 
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01   
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optvsa(DOUBLE *grdobj, DOUBLE *grdcon,INT init)
{
/*----------------------------------------------------------------------*/
  INT    i, iloc;               /* a counter */
  INT    kstep;
  INT numvar, numeq, numeq_total;
  DOUBLE act_loadfac, ddum;
  static  DOUBLE *svec;   /* vector with sensitivities on element level */
  static  DOUBLE *sveh;   /* necsessary for allreducing element values  */
/*----------------------------------------------------------------------*/
  static DIST_VECTOR  *grdisp;  /* vector for displacements derivatives */
/*----------------------------------------------------------------------*/
  INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */
  SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
  PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
  FIELD        *actfield;         /* pointer to the structural FIELD */
  INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
  CONTAINER     container;        /* contains variables defined in container.h */
  ELEMENT *actele;                /* active element                            */
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("optvsa");
  #endif
/*------------ the distributed system matrix, which is used for solving */
  actsysarray=0;
  numvar = 0; 
  if(opt->strategy==os_fsd) numvar = opt->strat.fsd->numvar;
/*--------------------------------------------------- set some pointers */
  actfield    = &(field[0]);
  actsolv     = &(solv[0]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  #ifdef PARALLEL 
  actintra    = &(par.intra[0]);
  #else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintra) dserror("Allocation of INTRA failed");
  actintra->intra_fieldtyp = structure;
  actintra->intra_rank   = 0;
  actintra->intra_nprocs   = 1;
  #endif
  container.fieldtyp  = actfield->fieldtyp;
  container.isdyn = 0;            /* static calculation */
  container.actndis      = 0;
/*--------------------------- initialize evaluation of sensitivities ---*/
  if(init==1)
  {
    svec  = (DOUBLE*)CCACALLOC(actfield->dis[0].numele,sizeof(DOUBLE));
    #ifdef PARALLEL 
    sveh  = (DOUBLE*)CCACALLOC(actfield->dis[0].numele,sizeof(DOUBLE));
    #endif
    /*---------------- create 1 vector grdisp for displacements deriv. */
    /*----------------------- get global and local number of equations */
    if (statvar->nonlinear==1) 
    {
      solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                    actsolv->sysarray_typ[actsysarray],
                    &numeq,
                    &numeq_total);
      solserv_create_vec(&(grdisp),1,numeq_total,numeq,"DV");
      solserv_zero_vec(&(grdisp[0]));
    }
    goto end;
  }
/*----------------------------------------------------------------------*/
/*-------------------------- evaluate sensitivities on element level ---*/
  /* derivative of strain energy     */
  if(opt->objective == oj_strain_energy)  *action = calc_struct_dee; 
  /* derivative of eigen frequencies */
  if(opt->objective == oj_frequency    )  *action = calc_struct_def; 

  container.dvec         = NULL;
  container.dirich       = NULL;
  container.global_numeq = 0;
  container.kstep        = 0;
  container.actndis      = 0;
  container.getvalue      = 0.;
  container.getvector     = svec;
  for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
#ifdef PARALLEL 
  for (i=0; i<actfield->dis[0].numele; i++) sveh[i]=0.;
#endif
  
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
    
#ifdef PARALLEL 
   /*---------------------------------------- allreduce objective value */
    MPI_Allreduce(svec,sveh,actfield->dis[0].numele,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
  
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    iloc       = actele->optdata[0];
#ifdef PARALLEL 
    grdobj[iloc-1] +=  sveh[i]; 
#else
    grdobj[iloc-1] +=  svec[i]; 
#endif
  }
  /*----------------------------------------- nonlinear static analysis */
  if (statvar->nonlinear==1) 
  {/*geolinear==0*/
    /*----------------------------- init the dist sparse matrix to zero */
    /*              NOTE: Has to be called after solver_control(init=1) */
    solserv_zero_mat(
                   actintra,
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sysarray_typ[actsysarray])
                    );
    /*------------------- calculate tangential stiffness in actsysarray */
    kstep      = statvar->nstep; 
    *action = calc_struct_nlnstiff;
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = kstep;
    container.actndis       = 0;
    container.isdyn         = 0;            /* static calculation */
    calelm(actfield,        /* active field                          */
           actsolv,         /* active solver typ                     */
           actpart,         /* my partition of this field            */
           actintra,        /* my intra-comunicators                 */
           actsysarray,     /* system-stiffness matrix               */
           -1,              /*system-mass matrix (there is no)       */
           &container,      /* contains variables defined in container.h */
           action);         /* what to do                            */
    /*- copy original load vector from [actsysarray+1] to [actsysarray] */
    solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(actsolv->rhs[actsysarray]));
    /* get final load factor */
    act_loadfac = opt->nlndat;
    /*----------------------- if load controlled; multiply with factor  */
    solserv_scalarprod_vec(&(actsolv->rhs[actsysarray]),act_loadfac);
/*------------------------------------------ solve for incremental load */
    init=0;
    solver_control(  
                 actsolv,
                 actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(grdisp[0]),
               &(actsolv->rhs[actsysarray]),
                 init
              );


   /*-------------------------- put actual total displacements to nodes */
    /* save old displacements for next step */

    solserv_result_total(
                     actfield,
                     actintra,
                     &(grdisp[0]),
                     1, /*---------- attentione ---------*//*!!!!!!!!!!!!!!!!!*/
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
    /*  evaluate objectives,constraints and derivatives  */
    /*     on element level for selfadjoint problems     */
    *action = calc_deriv_self_adj;
    container.getvector     = svec;
    container.actndis       = 0;
    container.isdyn         = 0;            /* static calculation */
    for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
#ifdef PARALLEL 
    for (i=0; i<actfield->dis[0].numele; i++) sveh[i]=0.;
#endif
  
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
    
#ifdef PARALLEL 
   /*---------------------------------------- allreduce objective value */
    MPI_Allreduce(svec,sveh,actfield->dis[0].numele,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
#ifdef PARALLEL 
    for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
    for (i=0; i<actfield->dis[0].numele; i++)
    {
      actele     = &(actfield->dis[0].element[i]);
      if(actele->optdata==NULL) continue; /* element does not take part in opt. */
      if(actele->optdata[0]==0) continue; /* position in variable vector        */
      iloc       = actele->optdata[0];
      svec[iloc-1] +=  sveh[i]; 
    }
#endif
  
    for (i=0; i<actfield->dis[0].numele; i++)
    {
      actele     = &(actfield->dis[0].element[i]);
      if(actele->optdata==NULL) continue; /* element does not take part in opt. */
      if(actele->optdata[0]==0) continue; /* position in variable vector        */
      iloc       = actele->optdata[0];
#ifdef PARALLEL 
      ddum            =  svec[iloc-1] - grdobj[iloc-1];
      grdobj[iloc-1]  =  ddum; 
#else
      ddum            =  svec[i] - grdobj[iloc-1];
      grdobj[iloc-1]  =  ddum; 
#endif
    }
  /*--------------------------------------------------------------------*/
  /* perform next structural analysis in one step only */
  statvar->nstep = 1;
  /* set first and final value of curve for incr. */
  curve[0].value.a.da[0][0] = opt->nlndat; 
  /*--------------------------------------------------------------------*/
  }/*geolinear==0*/
/*----------------------------------------------------------------------*/
  /*constraint*/
  if(opt->oeqc[0].oeqc_type == mass)
  {
    *action = calc_struct_dmc;       /* derivative of mass constraint */
  }
  for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
#ifdef PARALLEL 
  for (i=0; i<actfield->dis[0].numele; i++) sveh[i]=0.;
#endif
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
#ifdef PARALLEL 
   /*---------------------------------------- allreduce objective value */
    MPI_Allreduce(svec,sveh,actfield->dis[0].numele,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
  
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    iloc       = actele->optdata[0];
#ifdef PARALLEL 
    grdcon[iloc-1] +=  sveh[i]; /* add masses of elements */
#else
    grdcon[iloc-1] +=  svec[i]; /* add masses of elements */
#endif
  }
  for (i=0; i<numvar; i++)
  {
    grdcon[i] =  -grdcon[i]/opt->totmas; 
  }

/*----------------------------------------------------------------------*/
  end: ;
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
return; 
} /* end of optvsa */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
