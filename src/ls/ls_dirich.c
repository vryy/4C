/*!----------------------------------------------------------------------
\file
\brief ls_dirich.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR        par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT                numcurve;
extern struct _CURVE     *curve;

extern struct _LS_DYNAMIC     *lsdyn;



/*!----------------------------------------------------------------------
\brief subroutine to set dirichlet conditions at the beginning of the
process

<pre>                                                            irhan 05/04
subroutine to set dirichlet conditions at the beginning of the
process
</pre>

*----------------------------------------------------------------------*/
void ls_initdirich(
  FIELD          *actfield
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



/*!----------------------------------------------------------------------
\brief subroutine to set dirichlet conditions at an intermediate time step
in the process

<pre>                                                            irhan 05/04
subroutine to set dirichlet conditions at an intermediate time step
in the process
</pre>

*----------------------------------------------------------------------*/
void ls_setdirich(
  FIELD           *actfield,
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
  T = lsdyn->acttime;
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
/*! @} (documentation module close)*/
#endif
