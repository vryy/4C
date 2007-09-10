/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_intg' which sets all integration
parameters for the actual beam element

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0711 - 685-6574
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief sets all integration parameters for actual beam element

<pre>                                                              fh 09/02
This routine sets the gauss-point coordinates and the corresponding
weighting factors for the actual Timoshenko-beam element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *data     B3_DATA   (o)  data for integration parameters
\param option    INT       (i)  initialization (0) or calculation (1)


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: b3_cal_ele(), b3_init()

*-----------------------------------------------------------------------*/
void b3_intg(ELEMENT   *ele,
            B3_DATA   *data,
            INT        option)
{
INT i, k;

DOUBLE zero  = 0.0;
DOUBLE one   = 1.0;
DOUBLE two   = 2.0;
DOUBLE three = 3.0;

static DOUBLE xg[6][6],wgt[6][6];
static DOUBLE xl[6][6],wlt[6][6];
static DOUBLE xr[3];
#ifdef DEBUG
dstrc_enter("b3_intg");
#endif
/*----------------------------------------------------------------------*/
if (option==0)
{                                                  /* initialize arrays */
  for (i=0; i<6; i++)
  {
    for (k=0; k<6; k++)
    {
       xg[i][k]  = 0.;
       wgt[i][k] = 0.;
       xl[i][k]  = 0.;
       wlt[i][k] = 0.;
    }
  }
/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    L I N E A R ELEMENTS               |
 |     GAUSS SAMPLING POINTS  AT     R COORDINATES                      |
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

/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    L I N E A R ELEMENTS               |
 |     LOBATTO SAMPLING POINTS  AT   S,T COORDINATES                    |
 |                              AND  CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |    to get the values of strains/stresses at the cross section        |
 |    boundaries, an alternative integration rule is used, which        |
 |    includes the boundary coordinates for integration in s and t      |
 |    direction: LOBATTO quadrature rule (Bronstein p.611)              |
 *----------------------------------------------------------------------*/
      xl[0][1]  = -1.              ;
      xl[1][1]  = 0.               ;
      xl[2][1]  = -xl[0][1]        ;
      xl[0][2]  = xl[0][1]         ;
      xl[1][2]  = -0.4472135955    ;
      xl[2][2]  = -xl[1][2]        ;
      xl[3][2]  = -xl[0][2]        ;

      wlt[0][1] = 0.33333333333333 ;
      wlt[1][1] = 1.33333333333333 ;
      wlt[2][1] = wlt[0][1]        ;
      wlt[0][2] = 0.16666666666667 ;
      wlt[1][2] = 0.83333333333333 ;
      wlt[2][2] = wlt[1][2]        ;
      wlt[3][2] = wlt[0][2]        ;

/* local r-coordinate of the element nodes (for 2- and 3-noded elements)*/
      xr[0]	= -1.		   ;
      xr[1]	=  1.		   ;
      xr[2]	=  0.		   ;

/*----------------------------------------------------------------------*/
}                                                  /* initialize arrays */
else
{

/*-------------------------------------- beam element values -----------*/
      for(i=0;i<ele->e.b3->nGP[0] ;i++)
      {
        data->xgrr[i] = xg[ i][ele->e.b3->nGP[0]-1];
        data->wgtr[i] = wgt[i][ele->e.b3->nGP[0]-1];
      }
      for(i=0;i<3;i++)
      {
        data->xlst[i] = xl[ i][1];
	data->wlst[i] = wlt[i][1];
	data->xr[i]   = xr[i];
      }
/*----------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_intg */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
#endif
