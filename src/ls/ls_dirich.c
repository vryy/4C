#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "ls_prototypes.h"



extern struct _PAR        par; 
extern struct _MATERIAL  *mat;
extern INT                numcurve;
extern struct _CURVE     *curve;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls_initdirich(
  FIELD          *actfield, 
  LS_DYNAMIC     *lsdyn
  )
{
  INT        i,j;
  INT        actcurve;
  DOUBLE     timefac[MAXTIMECURVE];
  DOUBLE     T=0.0;
  DOUBLE     acttimefac;
  DOUBLE     initval;
  GNODE     *actgnode;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_initdirich");
#endif  
/*----------------------------------------------------------------------*/
  
  /* check dirichlet conditions */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode  = &(actfield->dis[0].node[i]); 
    actgnode = actnode->gnode; 
    if (actgnode->dirich==NULL)
      continue;
    if (actgnode->dirich->dirich_type==dirich_none)
    { 
      for (j=0;j<actnode->numdf;j++)
      {
        if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
          continue;
        actcurve = actgnode->dirich->curve.a.iv[j];
        if(actcurve>numcurve)
          dserror("Load curve: actual curve > number defined curves\n");   
      }
    }
  }
  /*  set dirichlet conditions at time (0) for zero intial field */
  if (lsdyn->init==0)
  {
    /* values from time curve */
    for (actcurve=0;actcurve<numcurve;actcurve++)
    {
      dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
    }
    /* loop */
    for (i=0;i<actfield->dis[0].numnp;i++)
    {
      actnode  = &(actfield->dis[0].node[i]); 
      actgnode = actnode->gnode;      
      if (actgnode->dirich==NULL)
        continue;
      switch(actgnode->dirich->dirich_type)
      {
          case dirich_none:
            for (j=0;j<actnode->numdf;j++)
            {
              if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
                continue;
              actcurve = actgnode->dirich->curve.a.iv[j]-1;
              if (actcurve<0)
                acttimefac = ONE;
              else
                acttimefac = timefac[actcurve];
              initval  = actgnode->dirich->dirich_val.a.dv[j];               
              actnode->sol_increment.a.da[0][j] = initval*acttimefac;
              actnode->sol.a.da[0][j] = initval*acttimefac;
            }
            break;
          default:
            dserror("dirch_type unknown!\n");
      }
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_initdirich */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls_setdirich(
  FIELD           *actfield, 
  LS_DYNAMIC      *lsdyn,
  INT              pos
  )
{
  INT        i,j;
  INT        actcurve;
  DOUBLE     timefac[MAXTIMECURVE];
  DOUBLE     T;
  DOUBLE     acttimefac;
  DOUBLE     initval;
  GNODE     *actgnode;
  NODE      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls_setdirich");
#endif 
/*----------------------------------------------------------------------*/
  
  /* set time */
  T = lsdyn->time;
  /* values from time curve */
  for (actcurve=0;actcurve<numcurve;actcurve++)
  {
    dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
  }
  /* loop */
  for (i=0; i<actfield->dis[0].numnp; i++) 
  {
    actnode  = &(actfield->dis[0].node[i]); 
    actgnode = actnode->gnode;      
    if (actgnode->dirich==NULL)
      continue;
    switch(actgnode->dirich->dirich_type)
    {
        case dirich_none:
          for (j=0; j<actnode->numdf; j++)
          {
            if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
              continue;
            actcurve = actgnode->dirich->curve.a.iv[j]-1;
            if (actcurve < 0)
              acttimefac = ONE;
            else
              acttimefac = timefac[actcurve];
            initval  = actgnode->dirich->dirich_val.a.dv[j];               
            actnode->sol_increment.a.da[pos][j] = initval*acttimefac;	 
          }
          break;
        default:
          dserror("dirch_type unknown!\n");
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_settdirich */
#endif
