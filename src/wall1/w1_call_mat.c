#include "../headers/standardtypes.h"
#include "wall1.h"
/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 WALL_TYPE wtype,
                 double **bop,
                 int ip,       
                 double *stress,
                 double **d,
                 int istore)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_call_mat");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_lin_el:/*------------------------------- linear elastic ---*/
    w1_mat_linel(mat->m.lin_el->youngs,
                 mat->m.lin_el->possionratio,
                 wtype,
                 d);
  break;
  case m_pl_mises:/*--------------------- von mises material law ---*/
    w1_mat_plast_mises(mat->m.pl_mises->youngs,
                       mat->m.pl_mises->possionratio,
                       mat->m.pl_mises->ALFAT,
                       mat->m.pl_mises->Sigy,
                       mat->m.pl_mises->Hard,
                       mat->m.pl_mises->GF,
                       ele,
                       wtype,
                       bop,
                       ip,
                       stress,
                       d,
                       istore);
  break;
  case m_pl_dp:/*------------------------ von mises material law ---*/
    w1_mat_plast_dp(   mat->m.pl_dp->youngs,
                       mat->m.pl_dp->possionratio,
                       mat->m.pl_dp->ALFAT,
                       mat->m.pl_dp->Sigy,
                       mat->m.pl_dp->Hard,
                       mat->m.pl_dp->PHI,
                       ele,
                       wtype,
                       bop,
                       ip,
                       stress,
                       d,
                       istore);
  break;
  default:
    dserror(" unknown type of material law");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_call_mat */
/*----------------------------------------------------------------------*/
