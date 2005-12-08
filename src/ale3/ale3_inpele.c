/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains the routine 'ale3inp' which reads a 3d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/


/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/


#ifdef D_ALE


#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2.h"
#include "ale3.h"
#include "../ale2/ale2.h"


/*-----------------------------------------------------------------------*/
/*!
  \brief reads a 3d ale element from the input file

  This routine reads a 3d ale element from the input file

  \param *ele   ELEMENT  (o)   the element

  \return void

  \sa caling: ---; called by: inp_ale_field()

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3inp(
    ELEMENT       *ele
    )
{

  INT  i;
  INT  ierr=0;
  INT  lmtmp;


#ifdef DEBUG
  dstrc_enter("ale3inp");
#endif


  /* allocate the element */
  ele->e.ale3 = (ALE3*)CCACALLOC(1,sizeof(ALE3));
  if (ele->e.ale3==NULL) dserror("Allocation of element ALE3 failed");


  /* read stuff needed for ALE element */
  /* read the element nodes */
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
  frchk("TET4",&ierr);
  if (ierr==1)
  {
    ele->distyp = tet4;
    ele->numnp=4;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
    /* rearrange element node numbers for tet4 */
    lmtmp=ele->lm[0];
    ele->lm[0]=ele->lm[1];
    ele->lm[1]=lmtmp;
  }
  frchk("TET10",&ierr);
  if (ierr==1)
  {
    ele->distyp = tet10;
    ele->numnp=10;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
  }

  /* reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;

  /* read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1) dserror("Reading of ALE element failed");

  /* read the gaussian points */
  if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
  {
    frint_n("GP",&(ele->e.ale3->nGP[0]),3,&ierr);
    if (ierr!=1) dserror("Reading of ALE3 element failed\n");
  }

  /* read gaussian points for tetrahedral elements */
  if (ele->numnp==4 || ele->numnp==10)
  {
    frint("GP_TET",&(ele->e.ale3->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of ALE3 element failed\n");
  }

  /* read gaussian points for triangle elements */
  frint("JAC",&(ele->e.ale3->jacobi),&ierr);
  if (ierr!=1) dserror("Reading of ALE element failed");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3inp */



/*! @} (documentation module close)*/

#endif

