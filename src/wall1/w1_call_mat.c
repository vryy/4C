#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*-----------------------------------------------------------------------|
|      topic: kondense 3D conditions                              sh 7/02|
|             to plane stress/strain conditions                          |
|-----------------------------------------------------------------------*/
void w1mat_trans_down (double **d,/*current material matrix d14-d44     */
                     ELEMENT   *ele,                                        
                     WALL_TYPE wtype,
                     double  *stress,  /*actuel stress [4]              */
                     double  *strain,  /*actual strain [4]              */
                     double  *qn,
                     double  **bop);
/*-----------------------------------------------------------------------|
|      topic: blowing up plane stress/strain conditions           sh 7/02|
|             to 3D --> 3D-Material Law                                  |
|-----------------------------------------------------------------------*/
void w1mat_trans_up (double **d,/*current material matrix d14-d44       */
                     ELEMENT   *ele,                                        
                     WALL_TYPE wtype,
                     double  *stress,  /*actuel stress [4]              */
                     double  *strain,  /*actual strain [4]              */
                     double  *qn,
                     double  **bop);
/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 WALL_TYPE wtype,
                 double **bop,
                 double  *gop,
                 double  *alpha,
                 double **xjm,
                 int ip,       
                 double *stress,
                 double **d,
                 int istore,/* controls storing of new stresses to wa */
                 int newval)/* controls evaluation of new stresses    */
{
int i;
int j;
double disd2;
double disd[4];
/*----------------------------------------------------------------------*/
double strain[6];
double *qn;
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
/*--------------------------------------------------fh 03/02----------*/
    if (newval=1)
    {
       w1_disd(ele,bop,gop,alpha,wtype,disd);
       switch(wtype)
       {
       case plane_stress:
	  disd2=(disd[2]+disd[3]);
    	  stress[0]=d[0][0]*disd[0]+d[0][1]*disd[1]+d[0][2]*disd2;
    	  stress[1]=d[1][0]*disd[0]+d[1][1]*disd[1]+d[1][2]*disd2;
    	  stress[2]=d[2][0]*disd[0]+d[2][1]*disd[1]+d[2][2]*disd2;
	  stress[3]=0;
       break;
       
       case rotat_symmet:
	  disd2=(disd[2]+disd[3]);        	 
	  stress[0]=d[0][0]*disd[0]+d[0][1]*disd[1]+d[0][2]*disd2+d[0][3]*disd[4];
	  stress[1]=d[1][0]*disd[0]+d[1][1]*disd[1]+d[1][2]*disd2+d[1][3]*disd[4];
	  stress[2]=d[2][0]*disd[0]+d[2][1]*disd[1]+d[2][2]*disd2+d[2][3]*disd[4];
	  stress[3]=d[3][0]*disd[0]+d[3][1]*disd[1]+d[3][2]*disd2+d[3][3]*disd[4];
       }
    }    
/*------------------------------------------------------------------*/	
  
  
  break;
  case m_stvenpor:/*------------------------ porous linear elastic ---*/
    w1_mat_stvpor(mat, ele->e.w1->elewa->matdata, wtype, d);
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
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
  case m_pl_mises_3D:/*--------------------- von mises material law -> 3D---*/
    w1mat_trans_up(d,
                   ele,
                   wtype,
                   stress,
                   strain,
                   qn,
                   bop);
    /*mat_plast_mises_3D(mat->m.pl_mises->youngs,
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
                       newval);*/
    w1mat_trans_down(d,
                   ele,
                   wtype,
                   stress,
                   strain,
                   qn,
                   bop);
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
                       gop,
                       alpha,
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
                       gop,
                       alpha,
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

#if 1
/*----------------------------------------------------------------------*/
/* get density out of material law                        ah 06/02      */
/*----------------------------------------------------------------------*/
void w1_getdensity(MATERIAL   *mat, double *density)
{
#ifdef DEBUG 
dstrc_enter("w1_getdensity");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*-------------------------------------- linear elastic */
   *density = mat->m.stvenant->density;
break;
case m_pl_mises:/*--------------------------- von mises material law ---*/
  dserror("Ilegal typ of material for this element");
break;
case m_pl_dp:/*------------------------- drucker prager material law ---*/
   dserror("Ilegal typ of material for this element");
break;
default:
   dserror("Ilegal typ of material for this element");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_getdensity */
#endif
#endif
