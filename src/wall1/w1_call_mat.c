#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 WALL_TYPE wtype,
                 double **bop,
                 double **xjm,
                 int ip,       
                 double *stress,
                 double **d,
                 int istore,/* controls storing of new stresses to wa */
                 int newval)/* controls evaluation of new stresses    */
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_call_mat");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*------------------------------- linear elastic ---*/
    w1_mat_linel(mat->m.stvenant->youngs,
                 mat->m.stvenant->possionratio,
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
                       istore,
                       newval);
  break;
  case m_pl_dp:/*------------------- drucker prager material law ---*/
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
                       istore,
                       newval);
  break;
  case m_pl_epc:/*---------- elastoplastic concrete material law ---*/
    w1_mat_plast_epc(  mat->m.pl_epc->dens        ,
                       mat->m.pl_epc->youngs      ,
                       mat->m.pl_epc->possionratio,
                       mat->m.pl_epc->alfat       ,
                       mat->m.pl_epc->xsi         ,
                       mat->m.pl_epc->sigy        ,
                       mat->m.pl_epc->ftm         ,
                       mat->m.pl_epc->fcm         ,
                       mat->m.pl_epc->gt          ,
                       mat->m.pl_epc->gc          ,
                       mat->m.pl_epc->gamma1      ,
                       mat->m.pl_epc->gamma2      ,
                       mat->m.pl_epc->nstiff      ,
                       mat->m.pl_epc->maxreb      ,
                       mat->m.pl_epc->rebar       ,
                       mat->m.pl_epc->reb_area    ,
                       mat->m.pl_epc->reb_ang     ,
                       mat->m.pl_epc->reb_so      ,
                       mat->m.pl_epc->reb_ds      ,
                       mat->m.pl_epc->reb_rgamma  ,
                       mat->m.pl_epc->reb_dens    ,
                       mat->m.pl_epc->reb_alfat   ,
                       mat->m.pl_epc->reb_emod    ,
                       mat->m.pl_epc->reb_rebnue  ,
                       mat->m.pl_epc->reb_sigy    ,
                       mat->m.pl_epc->reb_hard    ,
                       ele,                      
                       wtype,                    
                       bop,
                       xjm,                      
                       ip,                       
                       stress,                   
                       d,                        
                       istore,                   
                       newval);                  
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
