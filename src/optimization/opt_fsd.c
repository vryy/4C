/*!----------------------------------------------------------------------
\file
\brief contains the routine 'optfsd',
       main control routine for oc-methods in fsd

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
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
/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01   
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 | main control routine for oc-methods in fsd               al 05/01    |
 *----------------------------------------------------------------------*/
void optfsd()
{
/*----------------------------------------------------------------------*/
  INT  i;            /* a counter */
  INT  fsdstep, max_fsdstep;
  DOUBLE objective, constraint;
/*----------------------------------------------------------------------*/
  INT numvar, iprint, indcon;
  DOUBLE beta, accit, accitfsd, delta, etha, iota;
  DOUBLE diff, obj_old;
  DOUBLE *grdobj;
  DOUBLE *grdcon;
  DOUBLE *var   ;
  DOUBLE *resu  ;
  DOUBLE *resl  ;
  
  DOUBLE *etai  ;
  DOUBLE *varup ;
  DOUBLE *varlo ;
/*----------------------------------------------------------------------*/
  INT    gcounter, graph_out;  /* counter for graphical output AS*/
  INT    numvarlin;
  DOUBLE *grdobj_lin;              /* AS */
  DOUBLE *grdcon_lin;              /* AS */
  DOUBLE *var_lin   ;              /* AS */
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("optfsd");
  #endif
/*----------------------------------------------------------------------*/
  indcon = 1;/*number of independent constraints */
/*----------------------------------------------------------------------*/
  max_fsdstep   = opt->strat.fsd->numiter;
  numvar        = opt->strat.fsd->numvar;
  beta          = opt->strat.fsd->beta;
  accit         = opt->strat.fsd->acc;
  accitfsd      = opt->strat.fsd->alpha;
  delta         = opt->strat.fsd->delta;
  numvarlin     = opt->strat.fsd->numvar_lin; /* AS */
  graph_out     = opt->graph_out;             /* AS */
  
  etha     = 0.000001;         
  iota     = 1.0;              
  iprint   = 0;         
/*----------------------------------------------------------------------*/
  grdobj = opt->strat.fsd->grdobj;
  grdcon = opt->strat.fsd->grdcon;
  var    = opt->strat.fsd->var   ; 
  resu   = opt->strat.fsd->resu  ;
  resl   = opt->strat.fsd->resl  ; 
  grdobj_lin = opt->strat.fsd->grdobj_lin;  /* AS */
  grdcon_lin = opt->strat.fsd->grdcon_lin;  /* AS */
  var_lin    = opt->strat.fsd->var_lin   ;  /* AS */
/*----------------------------------------------------------------------*/
  etai  = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
  varup = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
  varlo = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
/*----------------------------------------------------------------------*/
  
/*------------------------------------------------ optimization loop ---*/
gcounter=0;
printf("graphical output requested every %d steps \n",graph_out);
for (fsdstep=0; fsdstep<max_fsdstep; fsdstep++)
{
/*---------------------------------------- initialize optimization step */
  for (i=0; i<numvar; i++) grdobj[i] =0.0;   
  for (i=0; i<numvar; i++) grdcon[i] =0.0;
  for (i=0; i<numvar; i++)   etai[i] =0.0;   
  for (i=0; i<numvar; i++)   varup[i]=0.0;
  for (i=0; i<numvar; i++)   varlo[i]=0.0;
  if(opt->numlin != 0)                         /* AS start */
  {
    for (i=0; i<numvarlin; i++) 
       { grdobj_lin[i] =0.0;  grdcon_lin[i] =0.0;}
  }                                            /* AS end   */
/*------------------------------ determine objective and constraints ---*/
  func(&indcon, &objective, &constraint);
/*------------------------- write discretization to graphical output ---*/
  if(fsdstep==0) opt_g_out(gr_mesh);
/*------------------------------------ write output of updated structure*/
  opt_g_out(gr_dens); 
  opt_g_out(gr_disp); 
/*------------------------------------ write output of updated structure*/
  gcounter++;                                 /* AS start */
  if (gcounter==graph_out)
  {
    opt_g_out(gr_dens); 
    opt_g_out(gr_disp); 
    gcounter = gcounter - graph_out;
    printf("\ngraphical output at optimisation step No. %d\n",fsdstep+1);
  }                                            /* AS end   */
/*---------------------------------------------------- print iterate ---*/
  if (par.myrank==0)
  {
    printf("----------------------------------------------------------------------------------------------\n");
    printf("Optimization step No. %d\n",fsdstep+1);
    printf("----------------------------------------------------------------------------------------------\n");
    printf("objective  function value: %12.3#E \n", objective );
    printf("constraint function value: %12.3#E \n", constraint);
  }
/*------------------------------------------- and check for convergence */
  if(fsdstep>0)
  {
    if(fabs(obj_old)<1E-10) obj_old = 1E-10;
    diff = fabs(objective-obj_old)/fabs(obj_old);
    if (diff<accit) break;
  }
  obj_old = objective;
/*--------------- determine derivatives of objective and constraints ---*/       
  optvsa(grdobj, grdcon,0);
/*------------------------------- perform smoothing of sensitivities ---*/
  if(opt->optsmooth == sm_on) optsmo(grdobj,0); 
/*------------------------------------------------ call oc-algorithm ---*/
  if(opt->numlin != 0)                         /* AS start */
  {
  linking(1, grdobj, grdcon, var, grdobj_lin, grdcon_lin, var_lin);   /* resize sens. vectors */
  fsdoc(var_lin,grdobj_lin,grdcon_lin,
        etai,&etha,&iota,resu,resl,varup,varlo,
        &numvarlin,&beta,&accitfsd,&delta,&iprint);
  linking(2, grdobj, grdcon, var, grdobj_lin, grdcon_lin, var_lin);  /* resize density vector */        
  }
  else
  {
  fsdoc(var,grdobj,grdcon,
        etai,&etha,&iota,resu,resl,varup,varlo,
        &numvar,&beta,&accitfsd,&delta,&iprint);
  }	
/*----------------------------------------------------- update structure*/
  /* will be done in func->optupd */
/*----------------------------------------------------------------------*/
}
/*----------------------------------------- end of optimization loop ---*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
return; 
} /* end of optfsd */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
