/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1inp' which reads the wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | read wall element                                         al 9/01    |
 *----------------------------------------------------------------------*/
void w1inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("winp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.w1 = (WALL1*)CCACALLOC(1,sizeof(WALL1));
if (ele->e.w1==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("QUAD4",&ierr);
if (ierr==1) 
{
   ele->distyp = quad4;
   ele->numnp=4;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD8",&ierr);
if (ierr==1) 
{
   ele->distyp = quad8;
   ele->numnp=8;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD9",&ierr);
if (ierr==1) 
{
   ele->distyp = quad9;
   ele->numnp=9;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*---------------------------------------------- read element thickness */
frdouble("THICK",&(ele->e.w1->thick),&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*-------------------------- read gaussian points for wall elements */
frint_n("GP",&(ele->e.w1->nGP[0]),2,&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*---------------------------------------------- read 2D problem type */
ele->e.w1->wtype = plane_stress; /* default */
frchk("PLANE_STRESS",&ierr);
if (ierr==1) ele->e.w1->wtype = plane_stress;
frchk("PLANE_STRAIN",&ierr);
if (ierr==1) ele->e.w1->wtype = plane_strain;
frchk("ROTATIONAL_SYMMETRY",&ierr);
if (ierr==1) ele->e.w1->wtype = rotat_symmet;
/*------------------------------------ read kinematic type (ah 06/02)*/
ele->e.w1->kintype = geo_lin;   /* default */
frchk("W_GeoLin",&ierr);
if (ierr==1) ele->e.w1->kintype = geo_lin;
frchk("W_TotalLagr",&ierr);
if (ierr==1) ele->e.w1->kintype = total_lagr;
frchk("W_UpdatedLagr",&ierr);
if (ierr==1) 
{
 ele->e.w1->kintype = updated_lagr;
 dserror("updated lagrange for WALL1 not implemented");
}
/*--------------------------------------------- read model for QUAD4 */
ele->e.w1->modeltype=displ_model;      /* default */
frchk("Displ_Model",&ierr);
if (ierr==1) ele->e.w1->modeltype=displ_model;
frchk("Incomp_Mode",&ierr);
if (ierr==1) ele->e.w1->modeltype=incomp_mode;
/*--------------------------------------- read local or global stresses */
frchar("STRESSES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"XY",2)==0)       ele->e.w1->stresstyp = w1_xy;
   if (strncmp(buffer,"RS",2)==0)       ele->e.w1->stresstyp = w1_rs;
}
if (ierr!=1) ele->e.w1->stresstyp = w1_xy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;

} /* end of w1inp */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
