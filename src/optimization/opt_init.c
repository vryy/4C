/*!----------------------------------------------------------------------
\file
\brief contains the routine 'opcini',
       initialize execution stage of optimization 

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
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
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
 |                                                         al 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 | initialize execution stage of optimization           a.lipka 5/01    |
 *----------------------------------------------------------------------*/
void opcini()
{
/*----------------------------------------------------------------------*/
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
/*----------------------------------------------------------------------*/
CONTAINER     container;        /* contains variables defined in container.h */
MATERIAL    *actmat;
/*----------------------------------------------------------------------*/
  int i, j, k;
  int numvar;                       /* number of optimization variables */
  double dens;
  ELEMENT *actele;                  /* active element                   */
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("opcini");
  #endif
/*--------------------------------------------------- set some pointers */
  actfield    = &(field[0]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  container.fieldtyp  = actfield->fieldtyp;
/*-------------- init the element integration routines for optimization */
  *action = calc_struct_opt_init;
  calinit(actfield,actpart,action,&container);
/*------------------------------------------------ initialize solver ---*/
  if (statvar->geolinear==1) 
  {
    opt_calsta(calsta_init);
  }
/*--------------- initialize element arrays for sensitivity analysis ---*/        
/*----------------------------------------------- number of opt.var. ---*/
  numvar=0;
  for (i=0; i<opt->numvar; i++)
  {
    if(opt->ovar[i].ovatt==eleofmat)
    {
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if(actele->mat==opt->ovar[i].objId)
         {
           numvar++;
           /* initialize element struct for opti. */
           actele->optdata = (int*)CCACALLOC(2,sizeof(int));
         }
      }
    
    }
  }
/*------------------------------------------------- initialize nlpql ---*/
  if (opt->strategy==os_nlp)
  {
  }
/*--------------------------------------------------- initialize fsd ---*/
  if (opt->strategy==os_fsd)
  {
    opt->strat.fsd->numvar  = numvar;
    
    opt->strat.fsd->grdobj  = (double*)CCACALLOC(numvar,sizeof(double));
    opt->strat.fsd->grdcon  = (double*)CCACALLOC(numvar,sizeof(double));
    opt->strat.fsd->var     = (double*)CCACALLOC(numvar,sizeof(double)); 
    opt->strat.fsd->resu    = (double*)CCACALLOC(numvar,sizeof(double));
    opt->strat.fsd->resl    = (double*)CCACALLOC(numvar,sizeof(double)); 
    for (i=0; i<numvar; i++)
    {
      opt->strat.fsd->grdobj[i]  =  0.;
      opt->strat.fsd->grdcon[i]  =  0.;
      opt->strat.fsd->var[i]     =  0.; 
      opt->strat.fsd->resu[i]    =  1.0E20;
      opt->strat.fsd->resl[i]    = -1.0E20; 
    }
  }

/*---------------------------------------------- initialize opt.var. ---*/
  numvar=0;
  for (i=0; i<opt->numvar; i++)
  {
    if(opt->ovar[i].ovatt==eleofmat)
    {
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if(actele->mat==opt->ovar[i].objId)
         {
           actmat = &(mat[actele->mat-1]);

           switch(actmat->mattyp)
           {
           case m_stvenant:/* ST.VENANT-KIRCHHOFF-MATERIAL */
              dens = actmat->m.stvenant->density;
           break;
           case m_neohooke:/* kompressible neo-hooke */
              dens = actmat->m.neohooke->density;
           break;
           case m_stvenpor:/* porous linear elastic ---*/
              dens = actmat->m.stvenpor->density;
           break;
           default:
              dserror("Ilegal typ of material");
           break;
           }
           /* element gets position in variable vector */
           actele->optdata[0] = numvar+1;
           /* and material Id of porous material */
           actele->optdata[1] = opt->ovar[i].objId;
           /* fill variable vector with initial values */
           opt->strat.fsd->var[numvar]  = dens;
           opt->strat.fsd->resu[numvar] = opt->ovar[i].bupper;
           opt->strat.fsd->resl[numvar] = opt->ovar[i].blower;

           numvar++;
           
         }
      }
    
    }
  }
/*-------------------- initialize upodate of optimization variables  ---*/
  optupd(1); 
/*-------------------- initialize evaluation of equality constraints ---*/
  opteqc(NULL,1); 
/*--------------------------- initialize evaluation of sensitivities ---*/
  optvsa(NULL,NULL,1); 
/*---------------------------- initialize smoothing of sensitivities ---*/
  if(opt->optsmooth == sm_on) optsmo(NULL,1); 
/*------------------ initialize reference values of mass, volume ... ---*/
  opt->totmas = 0.; 
  opt->totvol = 0.; 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of opcini */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
