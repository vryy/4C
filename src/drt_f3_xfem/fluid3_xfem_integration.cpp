/*!----------------------------------------------------------------------
\file fluid3_xfem_integration.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem_integration.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  evaluate the element integration points (private)        g.bau 03/07|
 *----------------------------------------------------------------------*/
INTEGRATION_POINTS_3D integration_points_3d(const  GaussRule3D gaussrule)
{
  const double Q12  = 1.0/2.0;
  const double Q14  = 1.0/4.0;
  const double Q16  = 1.0/6.0;
  const double Q124 = 1.0/6.0/4.0;
  const double Q430 = 4.0/5.0/6.0;
  const double Q9120= 9.0/4.0/5.0/6.0;

  const double xi2 = 0.5773502691896;
  const double xi3 = 0.7745966692415;

  const double w1 = 0.5555555555556;
  const double w2 = 0.8888888888889;
  const double w3 = w1;

  const double palpha = (5.0 + 3.0*sqrt(5.0))/20.0;
  const double pbeta  = (5.0 - sqrt(5.0))/20.0;

  INTEGRATION_POINTS_3D  intpoints;

  switch(gaussrule)
  {
    case hex_1point:
      intpoints.nquad = 1;
      intpoints.qxg[0][0] = 0.0;
      intpoints.qxg[0][1] = 0.0;
      intpoints.qxg[0][2] = 0.0;
      intpoints.qwgt[0] = 8.0;
      break;
    case hex_8point:
      intpoints.nquad = 8;
      intpoints.qxg[0][0] = -xi2; intpoints.qxg[0][1] = -xi2; intpoints.qxg[0][2] = -xi2;
      intpoints.qxg[1][0] =  xi2; intpoints.qxg[1][1] = -xi2; intpoints.qxg[1][2] = -xi2;
      intpoints.qxg[2][0] =  xi2; intpoints.qxg[2][1] =  xi2; intpoints.qxg[2][2] = -xi2;
      intpoints.qxg[3][0] = -xi2; intpoints.qxg[3][1] =  xi2; intpoints.qxg[3][2] = -xi2;
      intpoints.qxg[4][0] = -xi2; intpoints.qxg[4][1] = -xi2; intpoints.qxg[4][2] =  xi2;
      intpoints.qxg[5][0] =  xi2; intpoints.qxg[5][1] = -xi2; intpoints.qxg[5][2] =  xi2;
      intpoints.qxg[6][0] =  xi2; intpoints.qxg[6][1] =  xi2; intpoints.qxg[6][2] =  xi2;
      intpoints.qxg[7][0] = -xi2; intpoints.qxg[7][1] =  xi2; intpoints.qxg[7][2] =  xi2;
      intpoints.qwgt[0] = 1.0;
      intpoints.qwgt[1] = 1.0;
      intpoints.qwgt[2] = 1.0;
      intpoints.qwgt[3] = 1.0;
      intpoints.qwgt[4] = 1.0;
      intpoints.qwgt[5] = 1.0;
      intpoints.qwgt[6] = 1.0;
      intpoints.qwgt[7] = 1.0;
      break;
    case hex_27point:
      intpoints.nquad = 27;
      intpoints.qxg[0][0]  = -xi3; intpoints.qxg[0][1]  = -xi3; intpoints.qxg[0][2]  = -xi3;
      intpoints.qxg[1][0]  =  0.0; intpoints.qxg[1][1]  = -xi3; intpoints.qxg[1][2]  = -xi3;
      intpoints.qxg[2][0]  =  xi3; intpoints.qxg[2][1]  = -xi3; intpoints.qxg[2][2]  = -xi3;
      intpoints.qxg[3][0]  = -xi3; intpoints.qxg[3][1]  =  0.0; intpoints.qxg[3][2]  = -xi3;
      intpoints.qxg[4][0]  =  0.0; intpoints.qxg[4][1]  =  0.0; intpoints.qxg[4][2]  = -xi3;
      intpoints.qxg[5][0]  =  xi3; intpoints.qxg[5][1]  =  0.0; intpoints.qxg[5][2]  = -xi3;
      intpoints.qxg[6][0]  = -xi3; intpoints.qxg[6][1]  =  xi3; intpoints.qxg[6][2]  = -xi3;
      intpoints.qxg[7][0]  =  0.0; intpoints.qxg[7][1]  =  xi3; intpoints.qxg[7][2]  = -xi3;
      intpoints.qxg[8][0]  =  xi3; intpoints.qxg[8][1]  =  xi3; intpoints.qxg[8][2]  = -xi3;
      intpoints.qxg[9][0]  = -xi3; intpoints.qxg[9][1]  = -xi3; intpoints.qxg[9][2]  =  0.0;
      intpoints.qxg[10][0] =  0.0; intpoints.qxg[10][1] = -xi3; intpoints.qxg[10][2] =  0.0;
      intpoints.qxg[11][0] =  xi3; intpoints.qxg[11][1] = -xi3; intpoints.qxg[11][2] =  0.0;
      intpoints.qxg[12][0] = -xi3; intpoints.qxg[12][1] =  0.0; intpoints.qxg[12][2] =  0.0;
      intpoints.qxg[13][0] =  0.0; intpoints.qxg[13][1] =  0.0; intpoints.qxg[13][2] =  0.0;
      intpoints.qxg[14][0] =  xi3; intpoints.qxg[14][1] =  0.0; intpoints.qxg[14][2] =  0.0;
      intpoints.qxg[15][0] = -xi3; intpoints.qxg[15][1] =  xi3; intpoints.qxg[15][2] =  0.0;
      intpoints.qxg[16][0] =  0.0; intpoints.qxg[16][1] =  xi3; intpoints.qxg[16][2] =  0.0;
      intpoints.qxg[17][0] =  xi3; intpoints.qxg[17][1] =  xi3; intpoints.qxg[17][2] =  0.0;
      intpoints.qxg[18][0] = -xi3; intpoints.qxg[18][1] = -xi3; intpoints.qxg[18][2] =  xi3;
      intpoints.qxg[19][0] =  0.0; intpoints.qxg[19][1] = -xi3; intpoints.qxg[19][2] =  xi3;
      intpoints.qxg[20][0] =  xi3; intpoints.qxg[20][1] = -xi3; intpoints.qxg[20][2] =  xi3;
      intpoints.qxg[21][0] = -xi3; intpoints.qxg[21][1] =  0.0; intpoints.qxg[21][2] =  xi3;
      intpoints.qxg[22][0] =  0.0; intpoints.qxg[22][1] =  0.0; intpoints.qxg[22][2] =  xi3;
      intpoints.qxg[23][0] =  xi3; intpoints.qxg[23][1] =  0.0; intpoints.qxg[23][2] =  xi3;
      intpoints.qxg[24][0] = -xi3; intpoints.qxg[24][1] =  xi3; intpoints.qxg[24][2] =  xi3;
      intpoints.qxg[25][0] =  0.0; intpoints.qxg[25][1] =  xi3; intpoints.qxg[25][2] =  xi3;
      intpoints.qxg[26][0] =  xi3; intpoints.qxg[26][1] =  xi3; intpoints.qxg[26][2] =  xi3;
      intpoints.qwgt[0]  = w1*w1*w1;
      intpoints.qwgt[1]  = w2*w1*w1;
      intpoints.qwgt[2]  = w3*w1*w1;
      intpoints.qwgt[3]  = w1*w2*w1;
      intpoints.qwgt[4]  = w2*w2*w1;
      intpoints.qwgt[5]  = w3*w2*w1;
      intpoints.qwgt[6]  = w1*w3*w1;
      intpoints.qwgt[7]  = w2*w3*w1;
      intpoints.qwgt[8]  = w3*w3*w1;
      intpoints.qwgt[9]  = w1*w1*w2;
      intpoints.qwgt[10] = w2*w1*w2;
      intpoints.qwgt[11] = w3*w1*w2;
      intpoints.qwgt[12] = w1*w2*w2;
      intpoints.qwgt[13] = w2*w2*w2;
      intpoints.qwgt[14] = w3*w2*w2;
      intpoints.qwgt[15] = w1*w3*w2;
      intpoints.qwgt[16] = w2*w3*w2;
      intpoints.qwgt[17] = w3*w3*w2;
      intpoints.qwgt[18] = w1*w1*w3;
      intpoints.qwgt[19] = w2*w1*w3;
      intpoints.qwgt[20] = w3*w1*w3;
      intpoints.qwgt[21] = w1*w2*w3;
      intpoints.qwgt[22] = w2*w2*w3;
      intpoints.qwgt[23] = w3*w2*w3;
      intpoints.qwgt[24] = w1*w3*w3;
      intpoints.qwgt[25] = w2*w3*w3;
      intpoints.qwgt[26] = w3*w3*w3;
      break;
    case tet_1point:
      // GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1
      intpoints.nquad = 1;
      intpoints.qxg[0][0] =  Q14 ;
      intpoints.qxg[0][1] =  Q14 ;
      intpoints.qxg[0][2] =  Q14 ;
      intpoints.qwgt[0]   =  Q16 ;
    case tet_4point:
      // GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 2
      intpoints.nquad = 4;
      intpoints.qxg[0][0]    =    pbeta ;
      intpoints.qxg[1][0]    =    palpha;
      intpoints.qxg[2][0]    =    pbeta ;
      intpoints.qxg[3][0]    =    pbeta ;
      intpoints.qxg[0][1]    =    pbeta ;
      intpoints.qxg[1][1]    =    pbeta ;
      intpoints.qxg[2][1]    =    palpha;
      intpoints.qxg[3][1]    =    pbeta ;
      intpoints.qxg[0][2]    =    pbeta ;
      intpoints.qxg[1][2]    =    pbeta ;
      intpoints.qxg[2][2]    =    pbeta ;
      intpoints.qxg[3][2]    =    palpha;
      intpoints.qwgt[0]   =    Q124  ;
      intpoints.qwgt[1]   =    Q124  ;
      intpoints.qwgt[2]   =    Q124  ;
      intpoints.qwgt[3]   =    Q124  ;
    case tet_4point_alternative:
      // ALT.GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 1
      intpoints.qxg[0][0] = 0.0;
      intpoints.qxg[1][0] = 1.0;
      intpoints.qxg[2][0] = 0.0;
      intpoints.qxg[3][0] = 0.0;
      intpoints.qxg[0][1] = 0.0;
      intpoints.qxg[1][1] = 0.0;
      intpoints.qxg[2][1] = 1.0;
      intpoints.qxg[3][1] = 0.0;
      intpoints.qxg[0][2] = 0.0;
      intpoints.qxg[1][2] = 0.0;
      intpoints.qxg[2][2] = 0.0;
      intpoints.qxg[3][2] = 1.0;
      intpoints.qwgt[0]   = Q124;
      intpoints.qwgt[1]   = Q124;
      intpoints.qwgt[2]   = Q124;
      intpoints.qwgt[3]   = Q124;
      break;
    case tet_10point:
      // GAUSS INTEGRATION        5 SAMPLING POINTS, DEG.OF PRECISION 3
      intpoints.qxg[0][0] =     Q14  ;
      intpoints.qxg[1][0] =     Q12  ;
      intpoints.qxg[2][0] =     Q16  ;
      intpoints.qxg[3][0] =     Q16  ;
      intpoints.qxg[4][0] =     Q16  ;
      intpoints.qxg[0][1] =     Q14  ;
      intpoints.qxg[1][1] =     Q16  ;
      intpoints.qxg[2][1] =     Q16  ;
      intpoints.qxg[3][1] =     Q16  ;
      intpoints.qxg[4][1] =     Q12  ;
      intpoints.qxg[0][2] =     Q14  ;
      intpoints.qxg[1][2] =     Q16  ;
      intpoints.qxg[2][2] =     Q16  ;
      intpoints.qxg[3][2] =     Q12  ;
      intpoints.qxg[4][2] =     Q16  ;
      intpoints.qwgt[0]   =    -Q430 ;
      intpoints.qwgt[1]   =     Q9120;
      intpoints.qwgt[2]   =     Q9120;
      intpoints.qwgt[3]   =     Q9120;
      intpoints.qwgt[4]   =     Q9120;
      break;
    default:
      dserror("unknown integration rule");
  }

  return intpoints;
}


/*----------------------------------------------------------------------*
 |  evaluate the element integration points                             |
 *----------------------------------------------------------------------*/
void integration_points_2d(struct _INTEGRATION_POINTS_2D& intpoints,
                           const   GaussRule2D gaussrule)
{
    switch(gaussrule)
    {
    case quad_4point :
        intpoints.nquad = 4;
        intpoints.qwgt[0]  =  1.0;
        intpoints.qwgt[1]  =  1.0;
        intpoints.qwgt[2]  =  1.0;
        intpoints.qwgt[3]  =  1.0;
        
        intpoints.qxg[0][0] = -0.5773502691896;
        intpoints.qxg[0][1] = -0.5773502691896;
        intpoints.qxg[1][0] =  0.5773502691896;
        intpoints.qxg[1][1] = -0.5773502691896;
        intpoints.qxg[2][0] = -0.5773502691896;
        intpoints.qxg[2][1] =  0.5773502691896;
        intpoints.qxg[3][0] =  0.5773502691896;
        intpoints.qxg[3][1] =  0.5773502691896;
        break;
        
    case quad_9point:
        intpoints.nquad = 9; 
        intpoints.qwgt[0]  =  0.5555555555556*0.5555555555556;
        intpoints.qwgt[1]  =  0.8888888888889*0.5555555555556;
        intpoints.qwgt[2]  =  0.5555555555556*0.5555555555556;
        intpoints.qwgt[3]  =  0.5555555555556*0.8888888888889;
        intpoints.qwgt[4]  =  0.8888888888889*0.8888888888889;
        intpoints.qwgt[5]  =  0.5555555555556*0.8888888888889;
        intpoints.qwgt[6]  =  0.5555555555556*0.5555555555556;
        intpoints.qwgt[7]  =  0.8888888888889*0.5555555555556;
        intpoints.qwgt[8]  =  0.5555555555556*0.5555555555556;
        
        intpoints.qxg[0][0] = -0.7745966692415;
        intpoints.qxg[0][1] = -0.7745966692415;
        intpoints.qxg[1][0] =  0.0;
        intpoints.qxg[1][1] = -0.7745966692415;
        intpoints.qxg[2][0] =  0.7745966692415;
        intpoints.qxg[2][1] = -0.7745966692415;
        intpoints.qxg[3][0] = -0.7745966692415;
        intpoints.qxg[3][1] =  0.0; 
        intpoints.qxg[4][0] =  0.0; 
        intpoints.qxg[4][1] =  0.0;
        intpoints.qxg[5][0] =  0.7745966692415; 
        intpoints.qxg[5][1] =  0.0;
        intpoints.qxg[6][0] = -0.7745966692415; 
        intpoints.qxg[6][1] =  0.7745966692415; 
        intpoints.qxg[7][0] =  0.0;  
        intpoints.qxg[7][1] =  0.7745966692415; 
        intpoints.qxg[8][0] =  0.7745966692415;
        intpoints.qxg[8][1] =  0.7745966692415; 
        break;
        
    case tri_3point :
        intpoints.nquad = 3;                
        intpoints.qwgt[0]  = 1.0/6.0 ;
        intpoints.qwgt[1]  = 1.0/6.0 ;
        intpoints.qwgt[2]  = 1.0/6.0 ;

        intpoints.qxg[0][0] = 0.5;
        intpoints.qxg[0][1] = 0.0;
        intpoints.qxg[1][0] = 0.5;
        intpoints.qxg[1][1] = 0.5;
        intpoints.qxg[2][0] = 0.0;
        intpoints.qxg[2][1] = 0.5;
        break;
        
    case tri_6point:
        intpoints.nquad = 6;
        intpoints.qwgt[0]  = 0.0549758718277;
        intpoints.qwgt[1]  = 0.0549758718277;
        intpoints.qwgt[2]  = 0.0549758718277;
        intpoints.qwgt[3]  = 0.1116907948390;
        intpoints.qwgt[4]  = 0.1116907948390;
        intpoints.qwgt[5]  = 0.1116907948390;

        intpoints.qxg[0][0] = 0.0915762135098;
        intpoints.qxg[0][1] = 0.0915762135098;
        intpoints.qxg[1][0] = 0.8168475729805;
        intpoints.qxg[1][1] = 0.0915762135098;
        intpoints.qxg[2][0] = 0.0915762135098;
        intpoints.qxg[2][1] = 0.8168475729805;
        intpoints.qxg[3][0] = 0.4459484909160;
        intpoints.qxg[3][1] = 0.1081030181681;
        intpoints.qxg[4][0] = 0.4459484909160;
        intpoints.qxg[4][1] = 0.4459484909160;
        intpoints.qxg[5][0] = 0.1081030181681; 
        intpoints.qxg[5][1] = 0.4459484909160; 
        break;
    default:
        dserror("unknown integration rule");
    }

    return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
