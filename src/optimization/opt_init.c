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
 |                                                          al 08/02    |
 | pointer to allocate eigensolution variables                          |
 | dedfined in global_control.c                                         |
 | struct _ALLEIG       *alleig;                                        |
 *----------------------------------------------------------------------*/
extern ALLEIG              *alleig;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


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
  INT i, j;
  INT numvar;                       /* number of optimization variables */
  INT dId, desmat;
  INT *dsurfopt;
  INT *dvolopt;
  DOUBLE dens;
  ELEMENT *actele;                  /* active element                   */
/*----------------------------------------------------------------------*/
  INT c, c_max, mymat, nemat;
  INT numvarlin;                    /* AS */
  INT  *adr;
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
  if(opt->objective==oj_frequency)
  {
    calfrq(1);
    opt->oeig = (OEIG*)CCACALLOC(1,sizeof(OEIG));
    opt->oeig->numeigv = alleig->nroot;  /* number of eigenvalues in opt.process */
    opt->oeig->rhoks   = 5.0;            /* KREISSELMEIER-STEINHAUSER            */
  }
  else
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
           if(opt->numlin != 0)                              
	   {
	   actele->optdata = (INT*)CCACALLOC(3,sizeof(INT));
	   }                                                 
	   else actele->optdata = (INT*)CCACALLOC(2,sizeof(INT));
         }
      }
    
    }
  }
  /*-----------------------------------------------*/
  desmat = 0;
  for (i=0; i<opt->numvar; i++) if(opt->ovar[i].ovatt==eleofdesofmat) desmat = 1;
  /*-----------------------------------------------*/
  if(desmat > 0)
  {/*desmat*/
  /*-----------------------------------------------*/
  if(design->ndsurf> 0 )
  {
    dsurfopt = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));
    for (i=0; i<design->ndsurf; i++) dsurfopt[i] = 0;
  }
  if(design->ndvol > 0 )
  {
    dvolopt  = (INT*)CCACALLOC(design->ndvol ,sizeof(INT));
    for (i=0; i<design->ndvol ; i++)  dvolopt[i] = 0;
  }
  /*-----------------------------------------------*/
  for (i=0; i<opt->numvar; i++)
  {
    if(opt->ovar[i].ovatt==eleofdesofmat)
    {
      /*---------*/
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if(actele->mat==opt->ovar[i].objId)
         {
           /* initialize element struct for opti. */
           actele->optdata = (INT*)CCACALLOC(2,sizeof(INT));
           
           if(actele->eltyp==el_wall1)
           {
             dId = actele->g.gsurf[0].dsurf->Id;
             dsurfopt[dId] = 1;
           }
           if(actele->eltyp==el_brick1)
           {
             dId = actele->g.gvol[0].dvol->Id;
             dvolopt[dId] = 1;
           }
         
         }
      }
      /*---------*/
      /*
      for (j=0; j<design->ndsurf; j++) if(dsurfopt[dId]==1) numvar++;
      for (j=0; j<design->ndvol ; j++) if( dvolopt[dId]==1) numvar++;
      /*---------*/
      /* position in variable vector */
      for (j=0; j<design->ndsurf; j++)
      {
       if(dsurfopt[j]==1)
       {
         numvar++;
         dsurfopt[j]=numvar;
       }
      }
      for (j=0; j<design->ndvol; j++)
      {
       if(dvolopt[j]==1)
       {
         numvar++;
         dvolopt[j]=numvar;
       }
      }
      
      /*---------*/
    }
  }
  /*-----------------------------------------------*/
  /*
  if(design->ndsurf> 0 ) CCAFREE(dsurfopt);
  if(design->ndvol > 0 ) CCAFREE(dvolopt );
  /*-----------------------------------------------*/
  }/*desmat*/
/*------------------------------------------------- initialize nlpql ---*/
  if (opt->strategy==os_nlp)
  {
  }
/*--------------------------------------------------- initialize fsd ---*/
  if (opt->strategy==os_fsd)
  {
    opt->strat.fsd->numvar  = numvar;
    
    opt->strat.fsd->grdobj  = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
    opt->strat.fsd->grdcon  = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
    opt->strat.fsd->var     = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE)); 
    opt->strat.fsd->resu    = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE));
    opt->strat.fsd->resl    = (DOUBLE*)CCACALLOC(numvar,sizeof(DOUBLE)); 
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
           case m_mfoc:/* Porous open cell ST.VENANT-KIRCHHOFF-material */
              dens = actmat->m.mfoc->dens;
           break;
           case m_mfcc:/* Porous closed cell ST.VENANT-KIRCHHOFF-material */
              dens = actmat->m.mfcc->dens;
           break;
           case m_nhmfcc:/* foam, closed cell, based on modified Neo Hook */
              dens = actmat->m.nhmfcc->dens;
           break;
           default:
              dserror("Ilegal typ of material");
           break;
           }
           /* element gets position in variable vector */
           actele->optdata[0] = numvar+1;
	   if(opt->numlin != 0) actele->optdata[2] = numvar+1;
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
  /*--------------------------------------------------------------------*/
  if(desmat > 0)
  {/*desmat*/
  /*-----------------------------------------------*/
  numvar=0;
  for (i=0; i<opt->numvar; i++)
  {
    if(opt->ovar[i].ovatt==eleofdesofmat)
    {
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if(actele->mat==opt->ovar[i].objId)
         {
           actmat = &(mat[actele->mat-1]);

           switch(actmat->mattyp)
           {
           case m_stvenant:/* ST.VENANT-KIRCHHOFF-material */
              dens = actmat->m.stvenant->density;
           break;
           case m_mfoc:/* Porous open cell ST.VENANT-KIRCHHOFF-material */
              dens = actmat->m.mfoc->dens;
           break;
           case m_mfcc:/* Porous closed cell ST.VENANT-KIRCHHOFF-material */
              dens = actmat->m.mfcc->dens;
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
           /*----------------------------------------*/
           if(actele->eltyp==el_wall1)
           {
             /* element gets position in variable vector */
             dId = actele->g.gsurf[0].dsurf->Id;
             actele->optdata[0] = dsurfopt[dId];
             /* fill variable vector with initial values */
             opt->strat.fsd->var[ dsurfopt[dId]-1] = dens;
             opt->strat.fsd->resu[dsurfopt[dId]-1] = opt->ovar[i].bupper;
             opt->strat.fsd->resl[dsurfopt[dId]-1] = opt->ovar[i].blower;
           }
           if(actele->eltyp==el_brick1)
           {
             /* element gets position in variable vector */
             dId = actele->g.gvol[0].dvol->Id;
             actele->optdata[0] = dvolopt[dId];
             /* fill variable vector with initial values */
             opt->strat.fsd->var[ dvolopt[dId]-1] = dens;
             opt->strat.fsd->resu[dvolopt[dId]-1] = opt->ovar[i].bupper;
             opt->strat.fsd->resl[dvolopt[dId]-1] = opt->ovar[i].blower;
           }
           /* and material Id of porous material */
           actele->optdata[1] = opt->ovar[i].objId;
         }
      }
    
    }
  }
  
  /*-----------------------------------------------*/
  if(design->ndsurf> 0 ) CCAFREE(dsurfopt);
  if(design->ndvol > 0 ) CCAFREE(dvolopt );
  /*-----------------------------------------------*/
  }/*desmat*/
  /*--------------------------------------------------------------------*/
/*--- AS ---------------------------------------------------------------*/
/*---------------------------------------------- initialize opt.var. ---*/
if(opt->numlin != 0)
{
  numvarlin = numvar;
  adr = (int*)CCACALLOC(actfield->dis[0].numele,sizeof(int));
  for (i=0; i<opt->numlin; i++)  /* for all linking rules*/
  {
    mymat = opt->olin[i].objIds[0];
    nemat = opt->olin[i].objIds[1];
    c=0;
    for (j=0; j<actfield->dis[0].numele; j++)
    {
       actele = &(actfield->dis[0].element[j]);
       if(actele->mat == nemat)
       {
       actele->mylinweight = opt->olin[i].neweight;
       adr[c] = actele->optdata[2];
       c++;
       }
    }
    c_max = c;
    c=0;
    for (j=0; j<actfield->dis[0].numele; j++)
    {
       actele = &(actfield->dis[0].element[j]);
       if(actele->mat == mymat)
       {
       actele->mylinweight = opt->olin[i].myweight;
       actele->optdata[2] = adr[c];
       c++;
       }  
       if(c > c_max) dserror("Linking problem in opt_init.c, sorry!");
    }
   if(c != c_max) dserror("problem in opt_init.c: different no. of elements for linking materials!");
   numvarlin = numvarlin-c;      //angepasste laenge von numvar 
   }
   opt->strat.fsd->numvar_lin = numvarlin;
   CCAFREE(adr);  

/* -------------- allocate memory for resized vectors for linking  --- */
  opt->strat.fsd->grdobj_lin = (double*)CCACALLOC(numvarlin,sizeof(double));
  opt->strat.fsd->grdcon_lin = (double*)CCACALLOC(numvarlin,sizeof(double));
  opt->strat.fsd->var_lin    = (double*)CCACALLOC(numvarlin,sizeof(double));

  for (i=0; i<numvarlin; i++)
  {
     opt->strat.fsd->grdobj_lin[i]  =  0.;
     opt->strat.fsd->grdcon_lin[i]  =  0.;
     opt->strat.fsd->var_lin[i]     =  0.; 
  }
}  
/*--- AS ENDE ----------------------------------------------------------*/
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
