/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_inp' which reads the
       gradient enhanced wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifdef D_WALLGE
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief reads an gradient enhanced wall element from the input file

<pre>                                                        ah 05/03
This routine reads an gradient enhanced wall element from the input file
</pre>

\param *ele  ELEMENT  (o)   the element

\warning buffer[50] is not needed locally
\return void
\sa caling:    ---;
    called by: inp_struct_field()

*----------------------------------------------------------------------*/
void wge_inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
char buffer[50];
#ifdef DEBUG
dstrc_enter("wge_inp");
#endif
/*------------------------------------------------ allocate the element */
ele->e.wallge = (WALLGE*)CCACALLOC(1,sizeof(WALLGE));
if (ele->e.wallge==NULL) dserror("Allocation of element failed");
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
if (ierr!=1) dserror("Reading of WALLGE material failed");
/*---------------------------------------------- read element thickness */
frdouble("THICK",&(ele->e.wallge->thick),&ierr);
if (ierr!=1) dserror("Reading of WALLGE thickness failed");
/*-------------------------- read gaussian points for wall elements */
frint_n("GP",&(ele->e.wallge->nGP[0]),2,&ierr);
ele->e.wallge->nGP[1]=ele->e.wallge->nGP[0];
if (ierr!=1) dserror("Reading of WALLGE GP element failed");
/*---------------------------------------------- read 2D problem type */
ele->e.wallge->wgetype = pl_stress; /* default */
frchk("PLANE_STRESS",&ierr);
if (ierr==1) ele->e.wallge->wgetype = pl_stress;
frchk("PLANE_STRAIN",&ierr);
if (ierr==1) ele->e.wallge->wgetype = pl_strain;
frchk("ROTATIONAL_SYMMETRY",&ierr);
if (ierr==1) ele->e.wallge->wgetype = rot_symmet;
/*-------------------------- read if local or global stresses are asked */
frchar("STRESSES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"XY",2)==0)       ele->e.wallge->stresstyp = wge_xy;
   if (strncmp(buffer,"RS",2)==0)       ele->e.wallge->stresstyp = wge_rs;
}
if (ierr!=1) ele->e.wallge->stresstyp = wge_xy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of wge_inp */
/*----------------------------------------------------------------------*/
#endif /*D_WALLGE*/
/*! @} (documentation module close)*/
