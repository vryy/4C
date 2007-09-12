/*!----------------------------------------------------------------------
\file drt_utils_integration.cpp

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_utils_integration.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"


DRT::Utils::IntegrationPoints3D::IntegrationPoints3D(const GaussRule3D gaussrule)
{
  const double Q12  = 1.0/2.0;
  const double Q14  = 1.0/4.0;
  const double Q16  = 1.0/6.0;
  const double Q124 = 1.0/6.0/4.0;

  switch(gaussrule)
  {
  case intrule_hex_1point:
  {
    nquad = 1;
    qxg[0][0] = 0.0;
    qxg[0][1] = 0.0;
    qxg[0][2] = 0.0;
    qwgt[0] = 8.0;
    break;
  }
  case intrule_hex_8point:
  {
    nquad = 8;
    const double xi2 = 0.5773502691896;
    qxg[0][0] = -xi2; qxg[0][1] = -xi2; qxg[0][2] = -xi2;
    qxg[1][0] =  xi2; qxg[1][1] = -xi2; qxg[1][2] = -xi2;
    qxg[2][0] =  xi2; qxg[2][1] =  xi2; qxg[2][2] = -xi2;
    qxg[3][0] = -xi2; qxg[3][1] =  xi2; qxg[3][2] = -xi2;
    qxg[4][0] = -xi2; qxg[4][1] = -xi2; qxg[4][2] =  xi2;
    qxg[5][0] =  xi2; qxg[5][1] = -xi2; qxg[5][2] =  xi2;
    qxg[6][0] =  xi2; qxg[6][1] =  xi2; qxg[6][2] =  xi2;
    qxg[7][0] = -xi2; qxg[7][1] =  xi2; qxg[7][2] =  xi2;
    qwgt[0] = 1.0;
    qwgt[1] = 1.0;
    qwgt[2] = 1.0;
    qwgt[3] = 1.0;
    qwgt[4] = 1.0;
    qwgt[5] = 1.0;
    qwgt[6] = 1.0;
    qwgt[7] = 1.0;
    break;
  }
  case intrule_hex_27point:
  {
    nquad = 27;
    const double xi3 = 0.7745966692415;
    qxg[0][0]  = -xi3; qxg[0][1]  = -xi3; qxg[0][2]  = -xi3;
    qxg[1][0]  =  0.0; qxg[1][1]  = -xi3; qxg[1][2]  = -xi3;
    qxg[2][0]  =  xi3; qxg[2][1]  = -xi3; qxg[2][2]  = -xi3;
    qxg[3][0]  = -xi3; qxg[3][1]  =  0.0; qxg[3][2]  = -xi3;
    qxg[4][0]  =  0.0; qxg[4][1]  =  0.0; qxg[4][2]  = -xi3;
    qxg[5][0]  =  xi3; qxg[5][1]  =  0.0; qxg[5][2]  = -xi3;
    qxg[6][0]  = -xi3; qxg[6][1]  =  xi3; qxg[6][2]  = -xi3;
    qxg[7][0]  =  0.0; qxg[7][1]  =  xi3; qxg[7][2]  = -xi3;
    qxg[8][0]  =  xi3; qxg[8][1]  =  xi3; qxg[8][2]  = -xi3;
    qxg[9][0]  = -xi3; qxg[9][1]  = -xi3; qxg[9][2]  =  0.0;
    qxg[10][0] =  0.0; qxg[10][1] = -xi3; qxg[10][2] =  0.0;
    qxg[11][0] =  xi3; qxg[11][1] = -xi3; qxg[11][2] =  0.0;
    qxg[12][0] = -xi3; qxg[12][1] =  0.0; qxg[12][2] =  0.0;
    qxg[13][0] =  0.0; qxg[13][1] =  0.0; qxg[13][2] =  0.0;
    qxg[14][0] =  xi3; qxg[14][1] =  0.0; qxg[14][2] =  0.0;
    qxg[15][0] = -xi3; qxg[15][1] =  xi3; qxg[15][2] =  0.0;
    qxg[16][0] =  0.0; qxg[16][1] =  xi3; qxg[16][2] =  0.0;
    qxg[17][0] =  xi3; qxg[17][1] =  xi3; qxg[17][2] =  0.0;
    qxg[18][0] = -xi3; qxg[18][1] = -xi3; qxg[18][2] =  xi3;
    qxg[19][0] =  0.0; qxg[19][1] = -xi3; qxg[19][2] =  xi3;
    qxg[20][0] =  xi3; qxg[20][1] = -xi3; qxg[20][2] =  xi3;
    qxg[21][0] = -xi3; qxg[21][1] =  0.0; qxg[21][2] =  xi3;
    qxg[22][0] =  0.0; qxg[22][1] =  0.0; qxg[22][2] =  xi3;
    qxg[23][0] =  xi3; qxg[23][1] =  0.0; qxg[23][2] =  xi3;
    qxg[24][0] = -xi3; qxg[24][1] =  xi3; qxg[24][2] =  xi3;
    qxg[25][0] =  0.0; qxg[25][1] =  xi3; qxg[25][2] =  xi3;
    qxg[26][0] =  xi3; qxg[26][1] =  xi3; qxg[26][2] =  xi3;
    const double w1 = 0.5555555555556;
    const double w2 = 0.8888888888889;
    const double w3 = w1;
    qwgt[0]  = w1*w1*w1;
    qwgt[1]  = w2*w1*w1;
    qwgt[2]  = w3*w1*w1;
    qwgt[3]  = w1*w2*w1;
    qwgt[4]  = w2*w2*w1;
    qwgt[5]  = w3*w2*w1;
    qwgt[6]  = w1*w3*w1;
    qwgt[7]  = w2*w3*w1;
    qwgt[8]  = w3*w3*w1;
    qwgt[9]  = w1*w1*w2;
    qwgt[10] = w2*w1*w2;
    qwgt[11] = w3*w1*w2;
    qwgt[12] = w1*w2*w2;
    qwgt[13] = w2*w2*w2;
    qwgt[14] = w3*w2*w2;
    qwgt[15] = w1*w3*w2;
    qwgt[16] = w2*w3*w2;
    qwgt[17] = w3*w3*w2;
    qwgt[18] = w1*w1*w3;
    qwgt[19] = w2*w1*w3;
    qwgt[20] = w3*w1*w3;
    qwgt[21] = w1*w2*w3;
    qwgt[22] = w2*w2*w3;
    qwgt[23] = w3*w2*w3;
    qwgt[24] = w1*w3*w3;
    qwgt[25] = w2*w3*w3;
    qwgt[26] = w3*w3*w3;
    break;
  }
  case intrule_tet_1point:
  {
    // GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1
    nquad = 1;
    qxg[0][0] =  Q14 ;
    qxg[0][1] =  Q14 ;
    qxg[0][2] =  Q14 ;
    qwgt[0]   =  Q16 ;
    break;
  }
  case intrule_tet_4point:
  {
    // GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 2
    nquad = 4;
    const double palpha = (5.0 + 3.0*sqrt(5.0))/20.0;
    const double pbeta  = (5.0 - sqrt(5.0))/20.0;
    qxg[0][0]    =    pbeta ;
    qxg[1][0]    =    palpha;
    qxg[2][0]    =    pbeta ;
    qxg[3][0]    =    pbeta ;
    qxg[0][1]    =    pbeta ;
    qxg[1][1]    =    pbeta ;
    qxg[2][1]    =    palpha;
    qxg[3][1]    =    pbeta ;
    qxg[0][2]    =    pbeta ;
    qxg[1][2]    =    pbeta ;
    qxg[2][2]    =    pbeta ;
    qxg[3][2]    =    palpha;
    qwgt[0]   =    Q124  ;
    qwgt[1]   =    Q124  ;
    qwgt[2]   =    Q124  ;
    qwgt[3]   =    Q124  ;
    break;
  }
  case intrule_tet_4point_alternative:
  {
    // ALT.GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 1
    nquad = 4;
    qxg[0][0] = 0.0;
    qxg[1][0] = 1.0;
    qxg[2][0] = 0.0;
    qxg[3][0] = 0.0;
    qxg[0][1] = 0.0;
    qxg[1][1] = 0.0;
    qxg[2][1] = 1.0;
    qxg[3][1] = 0.0;
    qxg[0][2] = 0.0;
    qxg[1][2] = 0.0;
    qxg[2][2] = 0.0;
    qxg[3][2] = 1.0;
    qwgt[0]   = Q124;
    qwgt[1]   = Q124;
    qwgt[2]   = Q124;
    qwgt[3]   = Q124;
    break;
  }
  case intrule_tet_10point:
  {
    // GAUSS INTEGRATION        5 SAMPLING POINTS, DEG.OF PRECISION 3
    nquad = 10;
    qxg[0][0] =     Q14  ;
    qxg[1][0] =     Q12  ;
    qxg[2][0] =     Q16  ;
    qxg[3][0] =     Q16  ;
    qxg[4][0] =     Q16  ;
    qxg[0][1] =     Q14  ;
    qxg[1][1] =     Q16  ;
    qxg[2][1] =     Q16  ;
    qxg[3][1] =     Q16  ;
    qxg[4][1] =     Q12  ;
    qxg[0][2] =     Q14  ;
    qxg[1][2] =     Q16  ;
    qxg[2][2] =     Q16  ;
    qxg[3][2] =     Q12  ;
    qxg[4][2] =     Q16  ;
    const double Q430 = 4.0/5.0/6.0;
    const double Q9120= 9.0/4.0/5.0/6.0;
    qwgt[0]   =    -Q430 ;
    qwgt[1]   =     Q9120;
    qwgt[2]   =     Q9120;
    qwgt[3]   =     Q9120;
    qwgt[4]   =     Q9120;
    break;
  }

  case intrule_wedge_1point:
    {
      const double Q13 = 1.0/3.0;
      nquad = 1;
      qxg[0][0] = Q13;
      qxg[0][1] = Q13;
      qxg[0][2] = 0.0;
      qwgt[0] = 1;
      break;
    }

  case intrule_wedge_6point:
  {
    // GAUSS INTEGRATION         6 SAMPLING POINTS
    const double xi3 = 1.0/sqrt(3.0);
    const double Q23 = 2.0/3.0;
    const double Q16 = 1.0/6.0;

    nquad = 6;
    qxg[0][0] = Q23;
    qxg[1][0] = Q16;
    qxg[2][0] = Q16;
    qxg[3][0] = Q23;
    qxg[4][0] = Q16;
    qxg[5][0] = Q16;
    qxg[0][1] = Q16;
    qxg[1][1] = Q23;
    qxg[2][1] = Q16;
    qxg[3][1] = Q16;
    qxg[4][1] = Q23;
    qxg[5][1] = Q16;
    qxg[0][2] = xi3;
    qxg[1][2] = xi3;
    qxg[2][2] = xi3;
    qxg[3][2] = -xi3;
    qxg[4][2] = -xi3;
    qxg[5][2] = -xi3;

    qwgt[0] = Q16;
    qwgt[1] = Q16;
    qwgt[2] = Q16;
    qwgt[3] = Q16;
    qwgt[4] = Q16;
    qwgt[5] = Q16;
    break;
  }

   case intrule_wedge_9point:
    {
     // GAUSS INTEGRATION         9 SAMPLING POINTS
     const double xi3 = 0.77459666924;
    const double Q23 = 2.0/3.0;
    const double Q5913 = 5.0/9.0*1.0/3.0;

    nquad = 9;
    qxg[0][0] = Q23;
    qxg[1][0] = Q16;
    qxg[2][0] = Q16;
    qxg[3][0] = Q23;
    qxg[4][0] = Q16;
    qxg[5][0] = Q16;
    qxg[6][0] = Q23;
    qxg[7][0] = Q16;
    qxg[8][0] = Q16;
    qxg[0][1] = Q16;
    qxg[1][1] = Q23;
    qxg[2][1] = Q16;
    qxg[3][1] = Q16;
    qxg[4][1] = Q23;
    qxg[5][1] = Q16;
    qxg[6][1] = Q16;
    qxg[7][1] = Q23;
    qxg[8][1] = Q16;
    qxg[0][2] = xi3;
    qxg[1][2] = xi3;
    qxg[2][2] = xi3;
    qxg[3][2] = 0;
    qxg[4][2] = 0;
    qxg[5][2] = 0;
    qxg[6][2] = -xi3;
    qxg[7][2] = -xi3;
    qxg[8][2] = -xi3;

    qwgt[0] = Q5913;
    qwgt[1] = Q5913;
    qwgt[2] = Q5913;
    qwgt[3] = 0;
    qwgt[4] = 0;
    qwgt[5] = 0;
    qwgt[6] = -Q5913;
    qwgt[7] = -Q5913;
    qwgt[8] = -Q5913;
   break;
  }

  case intrule_pyramid_1point:
  {
    nquad = 1;
    qxg[0][0] = 0;
    qxg[1][0] = 0;
    qxg[2][0] = Q14;
    qwgt[0] = 4.0/3.0;
    break;
  }

  case intrule_pyramid_8point:
  {
    nquad = 8;
    qxg[0][0] = -0.26318405556971;
    qxg[1][0] = -0.50661630334979;
    qxg[2][0] = -0.26318405556971;
    qxg[3][0] = -0.50661630334979;
    qxg[4][0] = 0.26318405556971;
    qxg[5][0] = 0.50661630334979;
    qxg[6][0] = 0.26318405556971;
    qxg[7][0] = 0.50661630334979;
    qxg[0][1] = -0.26318405556971;
    qxg[1][1] = -0.50661630334979;
    qxg[2][1] = 0.26318405556971;
    qxg[3][1] = 0.50661630334979;
    qxg[4][1] = -0.26318405556971;
    qxg[5][1] = -0.50661630334979;
    qxg[6][1] = 0.26318405556971;
    qxg[7][1] = 0.50661630334979;
    qxg[0][2] = 0.54415184401122;
    qxg[1][2] = 0.12251482265544;
    qxg[2][2] = 0.54415184401122;
    qxg[3][2] = 0.12251482265544;
    qxg[4][2] = 0.54415184401122;
    qxg[5][2] = 0.12251482265544;
    qxg[6][2] = 0.54415184401122;
    qxg[7][2] = 0.12251482265544;


    qwgt[0] = 0.10078588207983;
    qwgt[1] = 0.23254745125351;
    qwgt[2] = 0.10078588207983;
    qwgt[3] = 0.23254745125351;
    qwgt[4] = 0.10078588207983;
    qwgt[5] = 0.23254745125351;
    qwgt[6] = 0.10078588207983;
    qwgt[7] = 0.23254745125351;
    break;
  }

  default:
    dserror("unknown 3D integration rule");
  }
}


DRT::Utils::IntegrationPoints2D::IntegrationPoints2D(const GaussRule2D gaussrule)
{
  switch(gaussrule)
  {
  case intrule_quad_1point :
  {
    nquad = 1;
    qwgt[0]  =  4.0;
    qxg[0][0] =  0.0;
    qxg[0][1] =  0.0;
    break;
  }
  case intrule_quad_4point :
  {
    nquad = 4;
    qwgt[0]  =  1.0;
    qwgt[1]  =  1.0;
    qwgt[2]  =  1.0;
    qwgt[3]  =  1.0;

    qxg[0][0] = -0.5773502691896;
    qxg[0][1] = -0.5773502691896;
    qxg[1][0] =  0.5773502691896;
    qxg[1][1] = -0.5773502691896;
    qxg[2][0] = -0.5773502691896;
    qxg[2][1] =  0.5773502691896;
    qxg[3][0] =  0.5773502691896;
    qxg[3][1] =  0.5773502691896;
    break;
  }
  case intrule_quad_9point:
  {
    nquad = 9;
    qwgt[0]  =  0.5555555555556*0.5555555555556;
    qwgt[1]  =  0.8888888888889*0.5555555555556;
    qwgt[2]  =  0.5555555555556*0.5555555555556;
    qwgt[3]  =  0.5555555555556*0.8888888888889;
    qwgt[4]  =  0.8888888888889*0.8888888888889;
    qwgt[5]  =  0.5555555555556*0.8888888888889;
    qwgt[6]  =  0.5555555555556*0.5555555555556;
    qwgt[7]  =  0.8888888888889*0.5555555555556;
    qwgt[8]  =  0.5555555555556*0.5555555555556;

    qxg[0][0] = -0.7745966692415;
    qxg[0][1] = -0.7745966692415;
    qxg[1][0] =  0.0;
    qxg[1][1] = -0.7745966692415;
    qxg[2][0] =  0.7745966692415;
    qxg[2][1] = -0.7745966692415;
    qxg[3][0] = -0.7745966692415;
    qxg[3][1] =  0.0;
    qxg[4][0] =  0.0;
    qxg[4][1] =  0.0;
    qxg[5][0] =  0.7745966692415;
    qxg[5][1] =  0.0;
    qxg[6][0] = -0.7745966692415;
    qxg[6][1] =  0.7745966692415;
    qxg[7][0] =  0.0;
    qxg[7][1] =  0.7745966692415;
    qxg[8][0] =  0.7745966692415;
    qxg[8][1] =  0.7745966692415;
    break;
  }
  case intrule_tri_1point:
  {
    nquad = 1;
    qwgt[0]  = 0.5;

    qxg[0][0] = 1.0/3.0;
    qxg[0][1] = 1.0/3.0;
    break;
  }
  case intrule_tri_3point:
  {
    // gp on triangular lines/edges
    nquad = 3;
    qwgt[0]  = 1.0/6.0 ;
    qwgt[1]  = 1.0/6.0 ;
    qwgt[2]  = 1.0/6.0 ;

    qxg[0][0] = 0.5;
    qxg[0][1] = 0.0;
    qxg[1][0] = 0.5;
    qxg[1][1] = 0.5;
    qxg[2][0] = 0.0;
    qxg[2][1] = 0.5;
    break;
  }
  case intrule_tri_3point_on_corners:
  {
    nquad = 3;
    qwgt[0]  = 1.0/6.0 ;
    qwgt[1]  = 1.0/6.0 ;
    qwgt[2]  = 1.0/6.0 ;

    qxg[0][0] = 1.0/6.0;
    qxg[0][1] = 1.0/6.0;
    qxg[1][0] = 2.0/3.0;
    qxg[1][1] = 1.0/6.0;
    qxg[2][0] = 1.0/6.0;
    qxg[2][1] = 2.0/3.0;
    break;
  }
  case intrule_tri_6point:
  {
    nquad = 6;
    qwgt[0]  = 0.0549758718277;
    qwgt[1]  = 0.0549758718277;
    qwgt[2]  = 0.0549758718277;
    qwgt[3]  = 0.1116907948390;
    qwgt[4]  = 0.1116907948390;
    qwgt[5]  = 0.1116907948390;

    qxg[0][0] = 0.0915762135098;
    qxg[0][1] = 0.0915762135098;
    qxg[1][0] = 0.8168475729805;
    qxg[1][1] = 0.0915762135098;
    qxg[2][0] = 0.0915762135098;
    qxg[2][1] = 0.8168475729805;
    qxg[3][0] = 0.4459484909160;
    qxg[3][1] = 0.1081030181681;
    qxg[4][0] = 0.4459484909160;
    qxg[4][1] = 0.4459484909160;
    qxg[5][0] = 0.1081030181681;
    qxg[5][1] = 0.4459484909160;
    break;
  }
  default:
    dserror("unknown 2D integration rule");
  }
}


DRT::Utils::IntegrationPoints1D::IntegrationPoints1D(const GaussRule1D gaussrule)
{
  switch(gaussrule)
  {
  case intrule_line_1point :
  {
    nquad = 1;
    qwgt[0]  =  2.0;
    qxg[0]   =  0.0;
    break;
  }
  case intrule_line_2point :
  {
    nquad = 2;
    qwgt[0]  =  1.0;
    qwgt[1]  =  1.0;
    qxg[0] = -0.5773502691896;
    qxg[1] =  0.5773502691896;
    break;
  }
  case intrule_line_3point:
  {
    nquad = 3;
    qwgt[0]  =  0.5555555555556;
    qwgt[1]  =  0.8888888888889;
    qwgt[2]  =  0.5555555555556;
    const double xi3 = 0.7745966692415;
    qxg[0] = -xi3;
    qxg[1] =  0.0;
    qxg[2] =  xi3;
    break;
  }
  default:
    dserror("unknown 1D integration rule");
  }
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
