/*======================================================================*/
/*!
\file
\brief Select proper material law

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>

\author mf
\date 10/06
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Select proper material law

\param *ele       ELEMENT   (i)   pointer to current element
\param *mat       MATERIAL  (i)   pointer to current material
\param **bop      DOUBLE    (i)   B-operator
\param ip         INT       (i)   current Gauss point index
\param *stress    DOUBLE    (o)   stress
\param **cmat     DOUBLE    (o)   constitutive matrix
\return void

\author mf
\date 10/06
*/
void so3_mat_sel(ELEMENT *ele,
                 MATERIAL *mat,
                 DOUBLE **bop,
                 INT ip,
                 DOUBLE *stress,
                 DOUBLE **cmat)
{
  INT imat,jmat,istrn,istss,inode;  /* counters */
  DOUBLE E,nu,mfac,strainsum,stresssum;
  DOUBLE strain[NUMSTRN_SOLID3];  /* strain vector */
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* the material law (it's a material world!) */
  /* soon here we want to connect to global material library */
  switch (mat->mattyp)
  {
    case struct_stvenant:
         /* !!! here we work hardwired for testing !!!*/
	 E   = mat->youngs;
	 nu  = mat->possionratio;
	 mfac=E/(1+nu)/(1-2*nu);
         /* constitutive matrix */
         /* set the whole thing to zero */
         for (imat=0; imat<6; imat++)
         { 
    	   for (jmat=0; jmat<6; jmat++)
	   {
	     cmat[imat][jmat] = 0.0;
	   }
         }
	 cmat[0][0] = mfac*(1-nu);
	 cmat[0][1] = mfac*nu;
	 cmat[0][2] = mfac*nu;
	 cmat[1][0] = mfac*nu;
	 cmat[1][1] = mfac*(1-nu);
	 cmat[1][2] = mfac*nu;
	 cmat[2][0] = mfac*nu;
	 cmat[2][1] = mfac*nu;
	 cmat[2][2] = mfac*(1-nu);
	 cmat[3][3] = mfac*0.5*(1-2*nu);
	 cmat[4][4] = mfac*0.5*(1-2*nu);
	 cmat[5][5] = mfac*0.5*(1-2*nu);
	 /* compute strains and stresses */
	 for (istrn=0; istrn<NUMSTRN_SOLID3; istrn++)
	 {
	   strainsum = 0.0;
	   for (inode=0; inode<ele->numnp; inode++)
	   {
	     strainsum += bop[istrn][inode] * ele->node[inode]->sol.a.da[0][0];
	   }
	   strain[istrn] = strainsum;
	 }
	 for (istss=0; istss<NUMSTSS_SOLID3; istss++)
	 {
	   stresssum = 0.0;
	   for (istrn=0; istrn<NUMSTRN_SOLID3; istrn++)
	   {
	     stresssum += cmat[istss][istrn] * strain[istrn];
	   }
	   stress[istress] = stresssum;
	 }
      break;
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_mat_sel */


/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
