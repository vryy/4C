#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "xfem_prototypes.h"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL      *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC *fdyn;


static INT     PREDOF=2;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_f2_calset( 
  ELEMENT            *ele,     
  DOUBLE            **xyze,
  DOUBLE            **eveln,    
  DOUBLE            **evelng,   
  DOUBLE             *epren,
  DOUBLE             *edeadn,
  DOUBLE             *edeadng,
  INT                *hasext
  )
{
  INT      i;
  INT      iel;
  NODE    *actnode;
  INT      actmat;
  INT      actcurve;
  DOUBLE   acttimefac;
  DOUBLE   dens;
  GSURF   *actgsurf;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_calset");
#endif
/*---------------------------------------------------------------------*/
  
  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  iel = ele->numnp;

  /* set element coordinates */
  for(i=0; i<iel; i++)
  {
    xyze[0][i] = ele->node[i]->x[0];
    xyze[1][i] = ele->node[i]->x[1];
  }
  
/*---------------------------------------------------------------------*
  | position of the different solutions:                               |
  | node->sol_incement: solution history used for calculations         |
  |	 sol_increment[1][i]: solution at (n)                          |
  |	 sol_increment[3][i]: solution at (n+1)                        |
  *--------------------------------------------------------------------*/
  for(i=0; i<iel; i++) /*-------------- loop nodes of element */
  {
    actnode=ele->node[i];
    /* set element velocities (n+gamma) */      
    /* => standart */
    evelng[0][i]=actnode->sol_increment.a.da[3][0];
    evelng[1][i]=actnode->sol_increment.a.da[3][1];
    /* => enriched */
    evelng[0][i+iel]=actnode->sol_increment.a.da[3][3];
    evelng[1][i+iel]=actnode->sol_increment.a.da[3][4];
    /* set supported pressures (n+1) */   
  }

  /*
   * NOTE =>
   * the time forces at a certain time step
   * is to be computed only once
   */
  if(fdyn->nif!=0) 
  {
    for(i=0; i<iel; i++)
    {
      actnode=ele->node[i];
      /* set element velocities at (n) */
      /* => standart */
      eveln[0][i]=actnode->sol_increment.a.da[1][0];
      eveln[1][i]=actnode->sol_increment.a.da[1][1];
      /* => enriched */
      eveln[0][i+iel]=actnode->sol_increment.a.da[1][3];
      eveln[1][i+iel]=actnode->sol_increment.a.da[1][4];
      /* set pressures (n) */   
      epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
    }
  }

  /*  check for dead load */
  actgsurf = ele->g.gsurf;
  if (actgsurf->neum!=NULL)
  {
    actmat = ele->mat-1;
    dens = mat[actmat].m.fluid->density;
    actcurve = actgsurf->neum->curve-1;
    if (actcurve<0)
      acttimefac = ONE;
    else
      dyn_facfromcurve(actcurve,fdyn->acttime,&acttimefac);
    for (i=0; i<2; i++)
    {
      if (actgsurf->neum->neum_onoff.a.iv[i]==0)
      {
        edeadn[i] = ZERO;
        edeadng[i] = ZERO;
      }
      if (actgsurf->neum->neum_type==neum_dead &&
          actgsurf->neum->neum_onoff.a.iv[i]!=0)
      {
        edeadn[i] = actgsurf->neum->neum_val.a.dv[i] * acttimefac;
        edeadng[i] = actgsurf->neum->neum_val.a.dv[i] * acttimefac;
        (*hasext)++;
      }
    }
  }
  
/*---------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return; 
} /* end of xfem_f2_calset */
#endif

