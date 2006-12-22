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
\brief Calculate linear B-operator matrix at point (r,s,t) 
       (e.g. Gauss point)

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   **bop      DOUBLE  (o)  B-operator contribution at (r,s,t)

\return void

\author mf
\date 10/06
*/
void so3_bop_lin(INT enod,
                 DOUBLE **deriv,
                 DOUBLE xji[NUMDOF_SOLID3][NUMDOF_SOLID3],
                 DOUBLE boplin[NUMSTR_SOLID3][MAXNOD_SOLID3*NUMDOF_SOLID3])
{
  /*--------------------------------------------------------------------*/
  INT inod;  /* current node */
  INT nodeindex;  /* node-column in bop */
  DOUBLE N_X,N_Y,N_Z; /* derivative w.r. to global coords */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_bop_lin");
#endif

  /*--------------------------------------------------------------------*/
  /* construct operator B with global derivatives */
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    nodeindex=inod*NUMDOF_SOLID3;  /* adress node-column in bop */
    /* transform derivatives local to global */
    N_X = xji[0][0]*deriv[inod][0]
        + xji[0][1]*deriv[inod][1]
        + xji[0][2]*deriv[inod][2];
    N_Y = xji[1][0]*deriv[inod][0]
        + xji[1][1]*deriv[inod][1]
        + xji[1][2]*deriv[inod][2];
    N_Z = xji[2][0]*deriv[inod][0]
        + xji[2][1]*deriv[inod][1]
        + xji[2][2]*deriv[inod][2];

    /* disperse global derivatives to bop-lines */
    /* bop is arranged as usual (refer to script FE or elsewhere):
     * [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
     * [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
     * [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
     * [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
     * [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
     * [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ] */
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

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_bop_lin */


/*======================================================================*/
/*!
\brief Calculate non-linear B-operator matrix at point (r,s,t) 
       (e.g. Gauss point)

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   *egrdis    DOUBLE  (i)  
\param   **bopnl    DOUBLE  (o)  B-operator contribution at (r,s,t)

\return void

\author mf
\date 10/06
*/
void so3_bop(ELEMENT *ele,
             INT enod,
             DOUBLE **deriv,
             DOUBLE xji[NUMDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE defgrd[NUMDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE bopn[MAXDOF_SOLID3][NUMDOF_SOLID3],
             DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3])
{
  /*--------------------------------------------------------------------*/
  INT inod;  /* current node */
  INT idof;  /* element DOF index */
  DOUBLE N_X, N_Y, N_Z;  /* derivative w.r. to global coords */
  DOUBLE grdshp[NUMDFGR_SOLID3][MAXDOF_SOLID3];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_bop_nl");
#endif

  /*--------------------------------------------------------------------*/
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    /*------------------------------------------------------------------*/
    /* referential derivatives of shape functions 
     * by isoparametric concept */
    N_X = xji[0][0]*deriv[inod][0]
        + xji[0][1]*deriv[inod][1]
        + xji[0][2]*deriv[inod][2];
    N_Y = xji[1][0]*deriv[inod][0]
        + xji[1][1]*deriv[inod][1]
        + xji[1][2]*deriv[inod][2];
    N_Z = xji[2][0]*deriv[inod][0]
        + xji[2][1]*deriv[inod][1]
        + xji[2][2]*deriv[inod][2];
    /*------------------------------------------------------------------*/
    /* nodally stored B-operator */
    /*           [ ... | N_{,1}^k | ... ]
     *    Bn^T = [ ... | N_{,2}^k | ... ]
     *           [ ... | N_{,3}^k | ... ]
     */
    bopn[inod][0] = N_X;
    bopn[inod][1] = N_Y;
    bopn[inod][2] = N_Z;
    /* address node-column in operator matrix */
    idof_X = inod * NUMDOF_SOLID3;
    idof_Y = idof_X + 1;
    idof_Z = idof_X + 2;
    /*------------------------------------------------------------------*/
    /* construct material gradient of shape function matrix  B_L */ 
    /*        [ d_X           ]   [ N^1           | ... | N^k           ]
     *        [      d_Y      ]   [      N^1      | ... |      N^k      ]
     *        [           d_Z ]   [           N^1 | ... |           N^k ]
     *        [ ~~~  ~~~  ~~~ ]
     *        [ d_Y           ]
     *        [      d_X      ]
     *        [ ~~~  ~~~  ~~~ ] . 
     *        [      d_Z      ]
     *        [           d_Y ]
     *        [ ~~~  ~~~  ~~~ ]
     *        [           d_X ]
     *        [ d_Z           ]
     *        
     *  Bl   =         L         .                   N
     *    
     */
    bopl[0][idof_X] = N_X;
    bopl[0][idof_Y] = 0.0;
    bopl[0][idof_Z] = 0.0;
    bopl[1][idof_X] = 0.0;
    bopl[1][idof_Y] = N_Y;
    bopl[1][idof_Z] = 0.0;
    bopl[2][idof_X] = 0.0;
    bopl[2][idof_Y] = 0.0;
    bopl[2][idof_Z] = N_Z;
    /* ~~~ */
    bopl[3][idof_X] = N_Y;
    bopl[3][idof_Y] = 0.0;
    bopl[3][idof_Z] = 0.0;
    bopl[4][idof_X] = 0.0;
    bopl[4][idof_Y] = N_X;
    bopl[4][idof_Z] = 0.0;
    /* ~~~ */
    bopl[5][idof_X] = 0.0;
    bopl[5][idof_Y] = N_Z;
    bopl[5][idof_Z] = 0.0;
    bopl[6][idof_X] = 0.0;
    bopl[6][idof_Y] = 0.0;
    bopl[6][idof_Z] = N_Y;
    /* ~~~ */
    bopl[7][idof_X] = N_Z;
    bopl[7][idof_Y] = 0.0;
    bopl[7][idof_Z] = 0.0;
    bopl[8][idof_X] = 0.0;
    bopl[8][idof_Y] = 0.0;
    bopl[8][idof_Z] = N_X;
    /*------------------------------------------------------------------*/
    /* B-operator */
    /*------------------------------------------------------------------*/
    /* linear B-operator */
    if (ele->kintype == so3_geo_lin)
    {
      /*----------------------------------------------------------------*/
      /* linear B-operator */
      /* 
       *     [ ... | N_{,X}^k                       | ... ]
       *     [ ... |            N_{,Y}^k            | ... ]
       *     [ ... |                       N_{,Z}^k | ... ]
       * B = [ ~~~   ~~~~~~~    ~~~~~~~~   ~~~~~~~~   ~~~ ]
       *     [ ... | N_{,Y}^k   N_{,X}^k            | ... ]
       *     [ ... |            N_{,Z}^k   N_{,Y}^k | ... ]
       *     [ ... | N_{,Z}^k              N_{,X}^k | ....]
       */
      bop[0][idof_X] = N_X;
      bop[0][idof_Y] = 0.0;
      bop[0][idof_Z] = 0.0;
      bop[1][idof_X] = 0.0;
      bop[1][idof_Y] = N_Y;
      bop[1][idof_Z] = 0.0;
      bop[2][idof_X] = 0.0;
      bop[2][idof_Y] = 0.0;
      bop[2][idof_Z] = N_Z;
      /* ~~~ */
      bop[3][idof_X] = N_Y;
      bop[3][idof_Y] = N_X;
      bop[3][idof_Z] = 0.0;
      bop[4][idof_X] = 0.0;
      bop[4][idof_Y] = N_Z;
      bop[4][idof_Z] = N_Y;
      bop[5][idof_X] = N_Z;
      bop[5][idof_Y] = 0.0;
      bop[5][idof_Z] = N_X;   
    }
    else if (ele->kintype == so3_total_lagr)
    {
      /*----------------------------------------------------------------*/
      /* non-linear B-operator (maybe so-called, meaning
       * of B-operator is not so sharp in the non-linear realm) */
      /*
       * B = F . Bl
       *
       *      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
       *      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
       *      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
       * B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
       *      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
       *      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
       *      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
       *      [                                                         ]
       *      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
       *      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
       *      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
       *      [                                                         ]
       *      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
       *      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
       *      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
       */
      bop[0][idof_X] = defgrd[0][0]*N_X;
      bop[0][idof_Y] = defgrd[1][0]*N_X;
      bop[0][idof_Z] = defgrd[2][0]*N_X;
      bop[1][idof_X] = defgrd[0][1]*N_Y;
      bop[1][idof_Y] = defgrd[1][1]*N_Y;
      bop[1][idof_Z] = defgrd[2][1]*N_Y;
      bop[2][idof_X] = defgrd[0][2]*N_Z; 
      bop[2][idof_Y] = defgrd[1][2]*N_Z; 
      bop[2][idof_Z] = defgrd[2][2]*N_Z; 
      /* ~~~ */
      bop[3][idof_X] = defgrd[0][0]*N_Y + defgrd[0][1]*N_X;
      bop[3][idof_Y] = defgrd[1][0]*N_Y + defgrd[1][1]*N_X;
      bop[3][idof_Z] = defgrd[2][0]*N_Y + defgrd[2][1]*N_X;
      bop[4][idof_X] = defgrd[0][1]*N_Z + defgrd[0][2]*N_Y;
      bop[4][idof_Y] = defgrd[1][1]*N_Z + defgrd[1][2]*N_Y;
      bop[4][idof_Z] = defgrd[2][1]*N_Z + defgrd[2][2]*N_Y;
      bop[5][idof_X] = defgrd[0][2]*N_X + defgrd[0][0]*N_Z;
      bop[5][idof_Y] = defgrd[1][2]*N_X + defgrd[1][0]*N_Z;
      bop[5][idof_Z] = defgrd[2][2]*N_X + defgrd[2][0]*N_Z;
    }
    else
    {
      dserror("Cannot compute B-operator for chosen kinematics.\n");
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_bop_nl */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
