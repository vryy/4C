/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_cal_stress' which gets converged stresses
    of GP-IPWA to have them for the output
<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  the routine 'wge_cal_stress' which gets converged stresses
        of GP-IPWA to have them for the output

\param *ele             ELEMENT     (I)   my element
*----------------------------------------------------------------------*/
void wge_cal_stress(ELEMENT   *ele)
{
#ifdef D_WALLGE

INT                 ip,i;
INT                 lr,ls;           /* loopers over GP */
INT                 nir,nis;         /* number of GP */

#ifdef DEBUG
dstrc_enter("wge_cal_stress");
#endif
/*----------------------------------------------------------------------*/
nir     = ele->e.wallge->nGP[0];
nis     = ele->e.wallge->nGP[1];
/*=================================================== GP- loops ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
   for (ls=0; ls<nis; ls++)
   {
      ip++;
      for (i=0; i<4; i++)
      {
         ele->e.wallge->stress_GP.a.d3[0][i][ip] =
                        ele->e.wallge->elwa[0].iptwa[ip].sig[i];
      }
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_cal_stress */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
#endif
