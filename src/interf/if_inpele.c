/*!----------------------------------------------------------------------
\file
\brief contains the routine 'interf_inp' which reads the 
       1D interface element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h" 

/*! 
\addtogroup INTERF
*/
/*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief reads an 1D interface element from the input file

<pre>                                                              ah 05/03 
This routine reads an 1D interface element from the input file

</pre>
\param *ele  ELEMENT  (I)   the element

\return void                                               
*----------------------------------------------------------------------*/
void interf_inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("interf_inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.interf = (INTERF*)CCACALLOC(1,sizeof(INTERF));
if (ele->e.interf==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("QUAD4",&ierr);
if (ierr==1) 
{
   ele->distyp = quad4;
   ele->numnp=4;
   ele->e.interf->nGP = 2;
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
   ele->e.interf->nGP = 3;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of INTERFACE material failed");
/*---------------------------------------------- read element thickness */
frdouble("THICK",&(ele->e.interf->thick),&ierr);
if (ierr!=1) dserror("Reading of INTERFACE thickness failed");
/*-------------------------------------- read number of gaussian points */
frint_n("GP",&(ele->e.interf->nGP),1,&ierr);
if (ierr!=1) dserror("Reading of INTERFACE gausspoints failed");
/*---------------------------------------------- read 2D problem type */
ele->e.interf->iftype = plane; /* default:plane stress or plane strain */
frchk("PLANE_STRESS",&ierr);
if (ierr==1) ele->e.interf->iftype = plane;
frchk("PLANE_STRAIN",&ierr);
if (ierr==1) ele->e.interf->iftype = plane;
frchk("ROTATIONAL_SYMMETRY",&ierr);
if (ierr==1) ele->e.interf->iftype = rotat_sym;
/*-------------------------- read if local or global stresses are asked */
frchar("STRESSES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"XY",2)==0)       ele->e.interf->stresstyp = if_xy;
   if (strncmp(buffer,"TN",2)==0)       ele->e.interf->stresstyp = if_tn;
}
if (ierr!=1) ele->e.interf->stresstyp = if_xy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of interf_inp */
/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
