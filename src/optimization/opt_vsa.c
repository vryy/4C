/*!----------------------------------------------------------------------
\file
\brief contains the routine 'optvsa', 
       variational sensitivity analysis

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../headers/optimization.h"
/*! 
\addtogroup OPTIMIZATION 
*//*! @{ (documentation module open)*/



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
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optvsa(double *grdobj, double *grdcon,int init)
{
/*----------------------------------------------------------------------*/
  int    i, j, iloc, cc;               /* a counter */
  static  double *svec;   /* vector with sensitivities on element level */
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
  dstrc_enter("optvsa");
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
  container.fieldtyp  = actfield->fieldtyp;
  container.isdyn = 0;            /* static calculation */
/*--------------------------- initialize evaluation of sensitivities ---*/
  if(init==1)
  {
    svec  = (double*)CCACALLOC(actfield->dis[0].numele,sizeof(double));
    goto end;
  }
/*----------------------------------------------------------------------*/
/*-------------------------- evaluate sensitivities on element level ---*/
  /*objective*/
  if(opt->objective == oj_strain_energy)
  {
    *action = calc_struct_dee;       /* derivative of strain energy */
  }
  container.dvec         = NULL;
  container.dirich       = NULL;
  container.global_numeq = 0;
  container.kstep        = 0;
  
  container.getvalue      = 0.;
  container.getvector     = svec;
  for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
  
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
  
  cc=0;
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    iloc       = actele->optdata[0];
    grdobj[cc] =  container.getvector[iloc-1]; 
    cc++;
  }
  /*constraint*/
  if(opt->oeqc[0].oeqc_type == mass)
  {
    *action = calc_struct_dmc;       /* derivative of mass constraint */
  }
  for (i=0; i<actfield->dis[0].numele; i++) svec[i]=0.;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
  
  cc=0;
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    iloc       = actele->optdata[0];
    grdcon[cc] =  -container.getvector[iloc-1]/opt->totmas; 
    cc++;
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
