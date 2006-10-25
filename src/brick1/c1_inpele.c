/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1inp' which reads data for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief input of data for a 3D hex element

<pre>                                                              al 06/02
This routine reads data for an 3D-hex-element.

</pre>
\param *ele     ELEMENT    (o)   the element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: inp_struct_field()

*----------------------------------------------------------------------*/
void c1inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
char buffer[50];
#ifdef DEBUG
dstrc_enter("c1inp");
#endif
/*------------------------------------------------ allocate the element */
ele->e.c1 = (BRICK1*)CCACALLOC(1,sizeof(BRICK1));
if (ele->e.c1==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("HEX8",&ierr);
if (ierr==1)
{
   ele->distyp = hex8;
   ele->numnp=8;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX20",&ierr);
if (ierr==1)
{
   ele->distyp = hex20;
   ele->numnp=20;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX27",&ierr);
if (ierr==1)
{
   ele->numnp=27;
   ele->distyp = hex27;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX27",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET4",&ierr);
if (ierr==1)
{
   ele->numnp=4;
   ele->distyp = tet4;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET10",&ierr);
if (ierr==1)
{
   ele->numnp=10;
   ele->distyp = tet10;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.c1->nGP[0]),3,&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
/*-------------------------- read eas-information */
frint("HYB",&(ele->e.c1->nhyb),&ierr);
/*------------------------ read element formulation: linear, T.L., U.L. */
frint("FORM",&(ele->e.c1->form),&ierr);
/*--------------------------------------- read local or global stresses */
ele->e.c1->stresstyp = c1_nostr; /* default */
frchar("STRESSES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"GPXYZ",3)==0)       ele->e.c1->stresstyp = c1_gpxyz;
   if (strncmp(buffer,"GPRST",3)==0)       ele->e.c1->stresstyp = c1_gprst;
   if (strncmp(buffer,"GP123",3)==0)       ele->e.c1->stresstyp = c1_gp123;
   if (strncmp(buffer,"NDXYZ",3)==0)       ele->e.c1->stresstyp = c1_npxyz;
   if (strncmp(buffer,"NDRST",3)==0)       ele->e.c1->stresstyp = c1_nprst;
   if (strncmp(buffer,"ND123",3)==0)       ele->e.c1->stresstyp = c1_np123;
   if (strncmp(buffer,"NDEQS",3)==0)       ele->e.c1->stresstyp = c1_npeqs;
}
/*--------------------------------------------------- read TSI coupling */
#ifdef D_TSI
ele->e.c1->tsi_couptyp = tsi_coup_none;  /* default */
frchar("TSI_COUPTYP",buffer,&ierr);
if (ierr)
{
  if (strncmp(buffer,"None",4)==0)
  {
    ele->e.c1->tsi_couptyp = tsi_coup_none;
  }
  if (strncmp(buffer,"Thermconf",9)==0)
  {
    ele->e.c1->tsi_couptyp = tsi_coup_thermconf;
  }
  if (strncmp(buffer,"Thermcreate",11)==0)
  {
    ele->e.c1->tsi_couptyp = tsi_coup_thermcreate;
  }
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of c1inp */
#endif
/*! @} (documentation module close)*/

