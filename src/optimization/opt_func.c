/*!----------------------------------------------------------------------
\file
\brief contains the routine 'func',
       evaluate values of objective function and constraints

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../headers/optimization.h"
#include "opt_prototypes.h"
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
 |                                                         al 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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
/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01   
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | evaluate values of objective function and constraints                |
 *----------------------------------------------------------------------*/
void func(int *m, double *f, double *g)
{
/*----------------------------------------------------------------------*/
int indcon;
/*----------------------------------------------------------------------*/
double objective;
double constraint;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("func");
#endif
/*----------------------------------------------------------------------*/
(*m)=0;
/*---------------- update all arrays with actual values of variables ---*/
  optupd(0);
/*-------------------------------------------------- linear analysis ---*/
  opt_stalin(calsta_solve);
/*------------------------- evaluate objective and constraint values ---*/
  optobj(&objective);             

  f[0] = objective;
 
     
  indcon=1;/* number of global (independent) constraints ?*/
  if (opt->numeqc>0) opteqc(&constraint,0);                                    
  g[0] = constraint;

/*  if (opt->numiqc>0) optcon(indcon);  

  for (i=0;i<(*m);i++)
  {
    g[i] = opt->ov_constraint[i];
  }
*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of func */
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for evaluation of equality constraints               |
 *----------------------------------------------------------------------*/
void opteqc(double *constraint,int init)
{
/*----------------------------------------------------------------------*/
int ione=1;
static double refvolu; /* reference or initial volume */
static double refmass; /* reference or initial mass   */
double tmpobj;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */
SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
CONTAINER     container;        /* contains variables defined in container.h */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("opteqc");
#endif
/*----------------------------------------------------------------------*/
  if(init==1)
  {
    refvolu = 0.;
    refmass = 0.;
    goto end;
  }
/*------------ the distributed system matrix, which is used for solving */
  actsysarray=0;
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
  container.isdyn   = 0;            /* static calculation */
  container.actndis = 0;            /* only one discretisation */
/*--------------------------- volume and weight equality constraints ---*/
  if(opt->oeqc[0].oeqc_type == volume)
  {
    *action = calc_struct_stv;       /* calculate structural volume */
    
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = 0;
    
    container.getvalue     = 0.;
    
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
    
    constraint[0] = container.getvalue;
    /*----- store integral values (only for last element in group) ---*/
    /*--------------------------------- volume equality constraint ---*/
    /*
    if( opt->ov_feqv[0]<=0.0) opt->ov_feqv[0] = opt->totvol;
                              opt->ov_reqv[0] = opt->totvol;
    opt->ov_constraint[0] = 1.0-opt->ov_reqv[0]/opt->ov_feqv[0];
    /**/
  }
  if(opt->oeqc[0].oeqc_type == mass)
  {
    *action = calc_struct_stm;       /* calculate structural mass */
    
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = 0;
    container.actndis      = 0;            /* only one discretisation */
    
    container.getvalue     = 0.;
    
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
    
#ifdef PARALLEL 
   /*---------------------------------------- allreduce objective value */
    tmpobj = container.getvalue;
    MPI_Allreduce(&tmpobj,&constraint[0],ione,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
   /*------------------------------------------------------------------ */
    constraint[0] = container.getvalue;
   /*------------------------------------------------------------------ */
#endif

    if( refmass<=0.0)
    {
      refmass = constraint[0];
      opt->totmas = refmass; 
    }

    constraint[0] = 1.0-constraint[0]/refmass;
  }
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
/*----------------------------------------------------------------------*/
end: ;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of opteqc */
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for update of fem-arrays for act. variable           |
 *----------------------------------------------------------------------*/
void optupd(int init)
{
/*----------------------------------------------------------------------*/
  int i, j, cc;
  static  double *svec;   /* vector with values for element level */
/*----------------------------------------------------------------------*/
int           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
CONTAINER     container;        /* contains variables defined in container.h */
ELEMENT *actele;                /* active element                            */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("optupd");
#endif
/*------------ the distributed system matrix, which is used for solving */
  actsysarray=0;
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
  container.isdyn   = 0;            /* static calculation */
  container.actndis = 0;            /* only one discretisation */
  container.fieldtyp  = actfield->fieldtyp;
/*--------------------------- initialize update of element materials ---*/
  if(init==1)
  {
    if(opt->opttype == ot_topology_optimization)
    {
      svec  = (double*)CCACALLOC(actfield->dis[0].numele,sizeof(double));
    } 
    goto end;
  }
/*---------------------------------------------- update of variables ---*/
  updvar();
/*------------------------------------- update of design and fe-mesh ---*/
  if(opt->opttype == ot_shape_optimization)
  {
   /* deoupd();*/
/*--------- check for element variables and element loads for update ---*/
   /* s1_upd(actfield,0);*/
  } 
/*------------------------------------ update of materials - topoopt ---*/
  if(opt->opttype == ot_topology_optimization)
  {
    *action = update_struct_odens;
   /* calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,NULL,NULL,0,0,action);*/
    
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = 0;
    container.actndis      = 0;            /* only one discretisation */
    
    container.getvalue      = 0.;
    container.getvector     = svec;
    for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;

    cc=0;
    for (i=0; i<opt->numvar; i++)
    {
      if(opt->ovar[i].ovatt==eleofmat)
      {
        for (j=0; j<actfield->dis[0].numele; j++)
        {
          actele = &(actfield->dis[0].element[j]);
          if(actele->mat==opt->ovar[i].objId)
          {
            if(actele->optdata==NULL) continue; /* element does not take part in opt. */
            if(actele->optdata[0]==0) continue; /* position in variable vector        */
            svec[actele->optdata[0]-1] = opt->strat.fsd->var[cc];
            cc++;
          }
        }
      }
    }
    /*********/
    
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);

  } 
/*----------------------------------------- new rhs for updated mesh ---*/
 /* wird in stalin automatisch gemacht ! */
/*----------------------------------------------------------------------*/
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
} /* end of optupd */
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | direct update of optimization variables                              |
 *----------------------------------------------------------------------*/
void updvar(void)
{
/*----------------------------------------------------------------------*/
OSNLP *nlp;
/*----------------------------------------------------------------------*/
int iscvar;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("updvar");
#endif
/*----------------------------------------------------------------------*/
  if(opt->strategy == os_nlp)
  {
    nlp = opt->strat.nlp;
    iscvar = nlp->iscvar;
  }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of updvar */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                      a.lipka 5/01    |
 | control program for calculation of objective functions               |
 *----------------------------------------------------------------------*/
void optobj(double *objective)
{
/*----------------------------------------------------------------------*/
int ione=1;
double tmpobj;
/*----------------------------------------------------------------------*/
int           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
CONTAINER     container;        /* contains variables defined in container.h */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("optobj");
#endif
/*------------ the distributed system matrix, which is used for solving */
  actsysarray=0;
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

  container.isdyn   = 0;            /* static calculation */
  container.actndis = 0;            /* only one discretisation */
/*------------------------- set values of objective functions = zero ---*/
  objective[0] = 0.0; /* only one value, right now */
  container.fieldtyp  = actfield->fieldtyp;
/*----------------------------------------------------------------------*/
  if(opt->objective == oj_strain_energy)
  {
    *action = calc_struct_ste;
    
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = 0;
    container.actndis      = 0;            /* only one discretisation */
    
    container.getvalue     = 0.;
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
#ifdef PARALLEL 
   /*---------------------------------------- allreduce objective value */
    tmpobj = container.getvalue;
    MPI_Allreduce(&tmpobj,&objective[0],ione,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
   /*------------------------------------------------------------------ */
    objective[0] = container.getvalue;
   /*------------------------------------------------------------------ */
#endif
  }
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of optobj */

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
