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
  int  i;            /* a counter */
  int  fsdstep, max_fsdstep;
  double objective, constraint;
/*----------------------------------------------------------------------*/
  int numvar, iprint, indcon;
  double beta, accit, accitfsd, delta, etha, iota;
  double diff, obj_old;
  double *grdobj;
  double *grdcon;
  double *var   ;
  double *resu  ;
  double *resl  ;
  
  double *etai  ;
  double *varup ;
  double *varlo ;
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
  
  etha     = 0.000001;         
  iota     = 1.0;              
  iprint   = 0;         
/*----------------------------------------------------------------------*/
  grdobj = opt->strat.fsd->grdobj;
  grdcon = opt->strat.fsd->grdcon;
  var    = opt->strat.fsd->var   ; 
  resu   = opt->strat.fsd->resu  ;
  resl   = opt->strat.fsd->resl  ; 
/*----------------------------------------------------------------------*/
  etai  = (double*)CCACALLOC(numvar,sizeof(double));
  varup = (double*)CCACALLOC(numvar,sizeof(double));
  varlo = (double*)CCACALLOC(numvar,sizeof(double));
  for (i=0; i<numvar; i++) etai[i] =0.0;   
  for (i=0; i<numvar; i++) varup[i]=0.0;
  for (i=0; i<numvar; i++) varlo[i]=0.0;
/*----------------------------------------------------------------------*/
  
/*------------------------------------------------ optimization loop ---*/
for (fsdstep=0; fsdstep<max_fsdstep; fsdstep++)
{
/*---------------------------------------- initialize optimization step */
  for (i=0; i<numvar; i++) grdobj[i] =0.0;   
  for (i=0; i<numvar; i++) grdcon[i] =0.0;
/*------------------------------ determine objective and constraints ---*/
  func(&indcon, &objective, &constraint);
/*------------------------- write discretization to graphical output ---*/
  if(fsdstep==0) opt_g_out(gr_mesh);
/*------------------------------------ write output of updated structure*/
  opt_g_out(gr_dens); 
  opt_g_out(gr_disp); 
/*---------------------------------------------------- print iterate ---*/
  printf("----------------------------------------------------------------------------------------------\n");
  printf("Optimization step No. %d\n",fsdstep+1);
  printf("----------------------------------------------------------------------------------------------\n");
  printf("objective  function value: %12.3#E \n", objective );
  printf("constraint function value: %12.3#E \n", constraint);
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
  fsdoc(var,grdobj,grdcon,
        etai,&etha,&iota,
        resu,resl,varup,varlo,
        &numvar,&beta,&accitfsd,&delta,&iprint);
/*----------------------------------------------------- update structure*/
  /* will be done in func->optupd */
/*------------------------------------ write output of updated structure*/
/*  opt_g_out(gr_dens); */
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
