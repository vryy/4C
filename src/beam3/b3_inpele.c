/*!----------------------------------------------------------------------
\file
\brief contains the routines 'b3inp' which reads the beam element from
input file

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief reads all datas of the beam element from input file

<pre>                                                              fh 09/02
This routine reads all datas of the actual beam element from input file

</pre>
\param *ele      ELEMENT  (i/o)  actual element


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: inp_struct_field()

*----------------------------------------------------------------------*/
void b3inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
#ifdef DEBUG
dstrc_enter("b3inp");
#endif
/*------------------------------------------------ allocate the element */
ele->e.b3 = (BEAM3*)CCACALLOC(1,sizeof(BEAM3));
/*-------------- read elements topology --------------------------------*/
frchk("LIN2",&ierr);
if (ierr==1)
{
   ele->distyp = line2;
   ele->numnp=2;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("LIN2",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("LIN3",&ierr);
if (ierr==1)
{
   ele->distyp = line3;
   ele->numnp=3;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("LIN3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*------------- read coordinates of reference node for local z-axis ----*/
frdouble_n("REF",&(ele->e.b3->nref[0]),3,&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
/*------------- read the material number -------------------------------*/
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
/*------------- read cross section data --------------------------------*/
frdouble("AREA",&(ele->e.b3->area),&ierr);
ele->e.b3->width=ele->e.b3->area;
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("GS",&(ele->e.b3->gs),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("HEIGHT",&(ele->e.b3->height),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("IYY",&(ele->e.b3->iyy),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("IZZ",&(ele->e.b3->izz),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("IYZ",&(ele->e.b3->iyz),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frdouble("IT",&(ele->e.b3->it),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
frint("IKE",&(ele->e.b3->ike),&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
/*------------- read number of gaussian points for beam elements -------*/
frint_n("GP",&(ele->e.b3->nGP[0]),1,&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
/*------------- read hinge code (modified later on if necessary)--------*/
frint_n("CODE",&(ele->e.b3->hc[1]),12,&ierr);
if (ierr!=1) dserror("Reading of BEAM3 element failed");
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif

return;

} /* end of b3inp */
#endif
/*! @} (documentation module close)*/
