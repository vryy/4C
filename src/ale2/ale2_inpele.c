/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2inp' which reads a 2d ale element

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "ale2.h"

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief reads a 2d ale element from the input file

<pre>                                                              mn 06/02
This routine reads a 2d ale element from the input file.

</pre>
\param *ele  ELEMENT  (o)   the element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: inp_ale_field()

*----------------------------------------------------------------------*/
void ale2inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
#ifdef DEBUG
dstrc_enter("ale2inp");
#endif
/*------------------------------------------------ allocate the element */
ele->e.ale2 = (ALE2*)CCACALLOC(1,sizeof(ALE2));
/*----------------------------------- read stuff needed for ALE element */
/*---------------------------------------------- read the element nodes */
frchk("QUAD4",&ierr);
if (ierr==1)
{
   ele->distyp = quad4;
   ele->numnp=4;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD8",&ierr);
if (ierr==1)
{
   ele->distyp = quad8;
   ele->numnp=8;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD9",&ierr);
if (ierr==1)
{
   ele->distyp = quad9;
   ele->numnp=9;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI3",&ierr);
if (ierr==1)
{
   ele->distyp = tri3;
   ele->numnp=3;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI6",&ierr);
if (ierr==1)
{
   ele->distyp = tri6;
   ele->numnp=6;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}

/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of ALE2 element failed");
/*-------------------------------------------- read the gaussian points */
if (ele->numnp==4 || ele->numnp==8 || ele->numnp==9)
{
   frint_n("GP",&(ele->e.ale2->nGP[0]),2,&ierr);
   if (ierr!=1) dserror("Reading of ALE2 element failed: integration\n");
}
/*-------------------------- read gaussian points for triangle elements */
if (ele->numnp==3 || ele->numnp==6)
{
   frint("GP_TRI",&(ele->e.ale2->nGP[0]),&ierr);
   if (ierr!=1) dserror("Reading of ALE2 element failed: integration\n");
}
/*-------------------------- read gaussian points for triangle elements */
frint("JAC",&(ele->e.ale2->jacobi),&ierr);
if (ierr!=1) dserror("Reading of ALE element failed");

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2inp */
#endif
/*! @} (documentation module close)*/
