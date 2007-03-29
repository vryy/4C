/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_final' which finalises the WALL1 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*/
/*!
\brief finalise element, basically deallocate stuff
\author bborn
\date 03/07
*/
void w1_final(PARTITION *actpart,
              MATERIAL *mat)
{
  INT          i,k;
  INT          size_j;
  ELEMENT     *actele;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_final");
#endif

  /*--------------------------------------------------------------------*/
  /* loop all elements */
  for (i=0; i<actpart->pdis[0].numele; i++)
  {/*matplast00*/
    actele = actpart->pdis[0].element[i];
    if (actele->eltyp == el_wall1)
    {
      /*----------------------------------------------------------------*/
      /* nothing to deallocate in integration points */
      /*----------------------------------------------------------------*/
      /* deallocate the space for stresses */
      am4del(&(actele->e.w1->stress_GP));
      am4del(&(actele->e.w1->stress_ND));
      /*----------------------------------------------------------------*/
      /* finalise history for GEMM scheme*/
#ifdef GEMM
      am4del(&(actele->e.w1->b_bar_history));
      am4del(&(actele->e.w1->PK_history));
#endif
      /*----------------------------------------------------------------*/
      /* init info for multiscale */  
#ifdef D_MLSTRUCT
      /* nothing to be done */
#endif
      /*----------------------------------------------------------------*/
      /* finalise working array */
      /* porous St.Venant-Kirchhoff */
      if(mat[actele->mat-1].mattyp == m_stvenpor)
      {
        if (actele->e.w1->elewa->optdata != NULL)
        {
          CCAFREE(actele->e.w1->elewa->optdata);
        }
        if (actele->e.w1->elewa->matdata != NULL)
        {
          CCAFREE(actele->e.w1->elewa->matdata);
        }
      }  /* end if */
      /* damage */
      if(mat[actele->mat-1].mattyp == m_dam_mp )
      {
        if (actele->e.w1->elewa[0].ipwa != NULL)
        {
          CCAFREE(actele->e.w1->elewa[0].ipwa);
        }
      }  /* end if */
      /* for plasticity and 3D-damage */
      if(mat[actele->mat-1].mattyp == m_pl_mises ||
         mat[actele->mat-1].mattyp == m_pl_mises_3D ||  /*Stefan's mises 3D*/
         mat[actele->mat-1].mattyp == m_pl_dp ||
         mat[actele->mat-1].mattyp == m_pl_epc ||
         mat[actele->mat-1].mattyp == m_pl_epc3D ||
         mat[actele->mat-1].mattyp == m_damage )
      {/*matplast01*/
        size_j = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
        for (k=0; k<size_j; k++)
        {
          if (actele->e.w1->elewa[0].ipwa[k].qn != NULL)
          {
            CCAFREE(actele->e.w1->elewa[0].ipwa[k].qn);
          }
          /* additional values needed for condensation       sh 08/02*/
          if(mat[actele->mat-1].mattyp == m_pl_mises_3D)
          {
            if (actele->e.w1->elewa[0].ipwa[k].sigc != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].sigc);
            }
            if (actele->e.w1->elewa[0].ipwa[k].sigi != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].sigi);
            }
            if ( actele->e.w1->elewa[0].ipwa[k].epsi != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].epsi);
            }
            if (actele->e.w1->elewa[0].ipwa[k].di != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].di);
            }
          }
          if(mat[actele->mat-1].mattyp == m_pl_epc ||
             mat[actele->mat-1].mattyp == m_pl_epc3D)
          {
            if (actele->e.w1->elewa[0].ipwa[k].sigc != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].sigc);
            }
            if (actele->e.w1->elewa[0].ipwa[k].grad != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].grad);
            }
            if ( actele->e.w1->elewa[0].ipwa[k].dlam != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].dlam);
            }
            if (actele->e.w1->elewa[0].ipwa[k].sigi != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].sigi);
            }
            if (actele->e.w1->elewa[0].ipwa[k].epsi != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].epsi);
            }
            if (actele->e.w1->elewa[0].ipwa[k].di != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].di);
            }
            if (actele->e.w1->elewa[0].ipwa[k].rsig != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].rsig);
            }
            if (actele->e.w1->elewa[0].ipwa[k].reps != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].reps);
            }
            if (actele->e.w1->elewa[0].ipwa[k].repstn != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].repstn);
            }
            if (actele->e.w1->elewa[0].ipwa[k].ryip != NULL)
            {
              CCAFREE(actele->e.w1->elewa[0].ipwa[k].ryip);
            }
          }
        }
        if (actele->e.w1->elewa[0].ipwa != NULL)
        {
          CCAFREE(actele->e.w1->elewa[0].ipwa);
        }
      }/*matplast01*/
      /* incomaptible modes */
      if(actele->e.w1->modeltype == incomp_mode)
      {
        if (actele->e.w1->elewa[0].imodewa != NULL)
        {
          CCAFREE(actele->e.w1->elewa[0].imodewa);
        }
      }
      /*-------------------------------------------------------------------*/
      if (actele->e.w1->elewa != NULL)
      {
        CCAFREE(actele->e.w1->elewa);
      }
    }  /* end if */
  }/*matplast00*/

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of w1_final */

#endif /*D_WALL1*/
/*! @} (documentation module close)*/
