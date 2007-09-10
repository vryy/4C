/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wgeintg' which calculates the needed values
       for the integration points

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALLGE
#include "../headers/standardtypes.h"
#include "wallge.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | integration points                                        ah 6/03    |
 -----------------------------------------------------------------------|
 |     ngp[0] --> number of integration points r-direction              |
 |     ngp[1] --> number of integration points s-direction              |
 |     xgrr   --> gauss sampling points        r-direction              |
 |     xgss   --> gauss sampling points        s-direction              |
 |     wgtr   --> weighting factors            r-direction              |
 |     wgts   --> weighting factors            s-direction              |
 *----------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief  calculates the needed values for the integration points

\param *ele          ELEMENT     (I)   the element
\param *data         WALLGE_DATA (I/O) the element's integration data
\param  option       INT         (I)   write the data or get them?

\warning There is nothing special to this routine.
\return void

*----------------------------------------------------------------------*/
void wgeintg(ELEMENT       *ele,
             WALLGE_DATA   *data,
             INT            option)
{
INT i, k;

DOUBLE one   = 1.0;
DOUBLE two   = 2.0;

static DOUBLE xg[6][6],wgt[6][6];
#ifdef DEBUG
dstrc_enter("w1intg");
#endif
/*----------------------------------------------------------------------*/
if (option==0)
{                                                  /* initialize arrays */
  for (i=0; i<6; i++)
  {
    for (k=0; k<6; k++)
    {
       xg[i][k] = 0.;
       wgt[i][k] = 0.;
    }
  }
/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    R E C T A N G U L A R   ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/

      xg[0][1]  =  -0.5773502691896;
      xg[1][1]  =  -xg[0][1]       ;
      xg[0][2]  =  -0.7745966692415;
      xg[2][2]  =  -xg[0][2]       ;
      xg[0][3]  =  -0.8611363115941;
      xg[1][3]  =  -0.3399810435849;
      xg[2][3]  =  -xg[1][3]       ;
      xg[3][3]  =  -xg[0][3]       ;
      xg[0][4]  =  -0.9061798459387;
      xg[1][4]  =  -0.5384693101057;
      xg[3][4]  =  -xg[1][4]       ;
      xg[4][4]  =  -xg[0][4]       ;
      xg[0][5]  =  -0.9324695142032;
      xg[1][5]  =  -0.6612093864663;
      xg[2][5]  =  -0.2386191860832;
      xg[3][5]  =  -xg[2][5]       ;
      xg[4][5]  =  -xg[1][5]       ;
      xg[5][5]  =  -xg[0][5]       ;

      wgt[0][0] =  two             ;
      wgt[0][1] =  one             ;
      wgt[1][1] =  one             ;
      wgt[0][2] =  0.5555555555556 ;
      wgt[1][2] =  0.8888888888889 ;
      wgt[2][2] =  wgt[0][2]       ;
      wgt[0][3] =  0.3478548451375 ;
      wgt[1][3] =  0.6521451548625 ;
      wgt[2][3] =  wgt[1][3]       ;
      wgt[3][3] =  wgt[0][3]       ;
      wgt[0][4] =  0.2369268850562 ;
      wgt[1][4] =  0.4786286704994 ;
      wgt[2][4] =  0.5688888888889 ;
      wgt[3][4] =  wgt[1][4]       ;
      wgt[4][4] =  wgt[0][4]       ;
      wgt[0][5] =  0.1713244923792 ;
      wgt[1][5] =  0.3607615730481 ;
      wgt[2][5] =  0.4679139345727 ;
      wgt[3][5] =  wgt[2][5]       ;
      wgt[4][5] =  wgt[1][5]       ;
      wgt[5][5] =  wgt[0][5]       ;

/*
      xg[0][1]  =  -1.0/sqrt(3.0);
      xg[1][1]  =  -xg[0][1]       ;
      xg[0][2]  =  -sqrt(3.0/5.0);
      xg[2][2]  =  -xg[0][2]       ;
      xg[0][3]  =  -0.8611363115941;
      xg[1][3]  =  -0.3399810435849;
      xg[2][3]  =  -xg[1][3]       ;
      xg[3][3]  =  -xg[0][3]       ;
      xg[0][4]  =  -0.9061798459387;
      xg[1][4]  =  -0.5384693101057;
      xg[3][4]  =  -xg[1][4]       ;
      xg[4][4]  =  -xg[0][4]       ;
      xg[0][5]  =  -0.9324695142032;
      xg[1][5]  =  -0.6612093864663;
      xg[2][5]  =  -0.2386191860832;
      xg[3][5]  =  -xg[2][5]       ;
      xg[4][5]  =  -xg[1][5]       ;
      xg[5][5]  =  -xg[0][5]       ;

      wgt[0][0] =  two             ;
      wgt[0][1] =  one             ;
      wgt[1][1] =  one             ;
      wgt[0][2] =  5.0/9.0;
      wgt[1][2] =  8.0/9.0;
      wgt[2][2] =  wgt[0][2]       ;
      wgt[0][3] =  0.3478548451375 ;
      wgt[1][3] =  0.6521451548625 ;
      wgt[2][3] =  wgt[1][3]       ;
      wgt[3][3] =  wgt[0][3]       ;
      wgt[0][4] =  0.2369268850562 ;
      wgt[1][4] =  0.4786286704994 ;
      wgt[2][4] =  0.5688888888889 ;
      wgt[3][4] =  wgt[1][4]       ;
      wgt[4][4] =  wgt[0][4]       ;
      wgt[0][5] =  0.1713244923792 ;
      wgt[1][5] =  0.3607615730481 ;
      wgt[2][5] =  0.4679139345727 ;
      wgt[3][5] =  wgt[2][5]       ;
      wgt[4][5] =  wgt[1][5]       ;
      wgt[5][5] =  wgt[0][5]       ;
*/
/*----------------------------------------------------------------------*/
}
else
{
/*-------------------------------------- rectangualar element values ---*/
      for(i=0;i<ele->e.wallge->nGP[0] ;i++)
      {
        data->xgrr[i] = xg[ i][ele->e.wallge->nGP[0]-1];
        data->wgtr[i] = wgt[i][ele->e.wallge->nGP[0]-1];
      }

      for(i=0;i<ele->e.wallge->nGP[1] ;i++)
      {
        data->xgss[i] = xg[ i][ele->e.wallge->nGP[1]-1];
        data->wgts[i] = wgt[i][ele->e.wallge->nGP[1]-1];
      }
/*----------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of wgeintg */
/*----------------------------------------------------------------------*/
#endif /*D_WALLGE*/
/*! @} (documentation module close)*/
#endif
