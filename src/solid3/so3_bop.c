/*======================================================================*/
/*!
\file
\brief Calculates the (linear) B-operator matrix for SOLID3 element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"


/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Calculate B-operator matrix at point (r,s,t) (e.g. Gauss point)

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   **bop      DOUBLE  (o)  B-operator contribution at (r,s,t)

\return void

\author mf
\date 10/06
*/
void so3_lin_bop(const INT        enod,
                 const DOUBLE   **deriv,
                 const DOUBLE   **xji,
                       DOUBLE   bop[NUMSTRN_SOLID3][MAXNOD_SOLID3*NUMDOF_SOLID3])
{
  /*--------------------------------------------------------------------*/
  INT inod;  /* current node */
  INT nodeindex;  /* node-column in bop */
  DOUBLE N_X,N_Y,N_Z; /* derivative w.r. to global coords */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_lin_bop");
#endif

  /*--------------------------------------------------------------------*/
  /* construct operator B with global derivatives */
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    nodeindex=inod*NUMDOF_SOLID3;  /* adress node-column in bop */
    /* transform derivatives local to global */
    N_X=xji[0][0]*deriv[inod][0]+xji[0][1]*deriv[inod][1]+xji[0][2]*deriv[inod][2];
    N_Y=xji[1][0]*deriv[inod][0]+xji[1][1]*deriv[inod][1]+xji[1][2]*deriv[inod][2];
    N_Z=xji[2][0]*deriv[inod][0]+xji[2][1]*deriv[inod][1]+xji[2][2]*deriv[inod][2];

    /* disperse global derivatives to bop-lines */
    /* bop is arranged as usual (refer to script FE or elsewhere):
     * [N1,X  0  0  |N2,X  0  0  |... |Ni,X  0  0 ]
     * [0  N1,Y  0  |0  N2,Y  0  |... |0  Ni,Y  0 ]
     * [0  0  N1,Z  |0  0  N2,Z  |... |0  0  Ni,Z ]
     * [N1,Y N1,X 0 |N2,Y N2,X 0 |... |Ni,Y Ni,X 0]
     * [0 N1,Z N1,Y |0 N2,Z N2,Y |... |0 Ni,Z Ni,Y]
     * [N1,Z 0 N1,X |N2,Z 0 N2,X |... |Ni,Z 0 Ni,X]*/
    bop[0][nodeindex+0] = N_X;
    bop[0][nodeindex+1] = 0.0;
    bop[0][nodeindex+2] = 0.0;
    bop[1][nodeindex+0] = 0.0;
    bop[1][nodeindex+1] = N_Y;
    bop[1][nodeindex+2] = 0.0;
    bop[2][nodeindex+0] = 0.0;
    bop[2][nodeindex+1] = 0.0;
    bop[2][nodeindex+2] = N_Z;
    bop[3][nodeindex+0] = N_Y;
    bop[3][nodeindex+1] = N_X;
    bop[3][nodeindex+2] = 0.0;
    bop[4][nodeindex+0] = 0.0;
    bop[4][nodeindex+1] = N_Z;
    bop[4][nodeindex+2] = N_Y;
    bop[5][nodeindex+0] = N_Z;
    bop[5][nodeindex+1] = 0.0;
    bop[5][nodeindex+2] = N_X;
    
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_lin_bop */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
