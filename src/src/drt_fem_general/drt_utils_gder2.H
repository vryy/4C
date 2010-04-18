/*!
\file drt_utils_gder2.H

\brief calculate second global derivatives w.r.t x/y/z at point r,s,t

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifdef CCADISCRET
#ifndef UTILS_GDER2_H
#define UTILS_GDER2_H

#include "../drt_lib/linalg_fixedsizematrix.H"


namespace DRT
{
namespace UTILS
{

/*---------------------------------------------------------------------*/
/*!
 * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
 */
 /*                                                          gammi 07/07
 |
 | 3 space dimensions:
 |
 | From the six equations
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ------ = -- | --*-- + --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy   ds dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dt^2     dt | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | dt dr     dt | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dt     ds | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 | the matrix (jacobian-bar matrix) system
 |
 | +-                                                                                         -+   +-    -+
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
 | |                                                                                           | * |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
 | +-                                                                                         -+   +-    -+
 |
 |                  +-    -+     +-                           -+
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
 |              =   |      |  -  |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drds |     | drds dx   drds dy   drds dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drdt |     | drdt dx   drdt dy   drdt dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dtds |     | dtds dx   dtds dy   dtds dz |
 |                  +-    -+     +-                           -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 |
 | --------------------------------------------------------------
 | 2 space dimensions:
 |
 | From the three equations
 |
 |              +-             -+
 |  d^2N     d  | dx dN   dy dN |
 |  ----   = -- | --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy |
 |              +-             -+
 |
 |              +-             -+
 |  d^2N     d  | dx dN   dy dN |
 |  ------ = -- | --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy |
 |              +-             -+
 |
 |              +-             -+
 |  d^2N     d  | dx dN   dy dN |
 | -----   = -- | --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy |
 |              +-             -+
 |
 | the matrix system
 |
 | +-                                        -+   +-    -+
 | |   /dx\^2        /dy\^2         dy dx     |   | d^2N |
 | |  | -- |        | ---|        2*--*--     |   | ---- |
 | |   \dr/          \dr/           dr dr     |   | dx^2 |
 | |                                          |   |      |
 | |   /dx\^2        /dy\^2         dy dx     |   | d^2N |
 | |  | -- |        | -- |        2*--*--     | * | ---- |
 | |   \ds/          \ds/           ds ds     |   | dy^2 | =
 | |                                          |   |      |
 | |   dx dx         dy dy      dx dy   dy dx |   | d^2N |
 | |   --*--         --*--      --*-- + --*-- |   | ---- |
 | |   dr ds         dr ds      dr ds   dr ds |   | dxdy |
 | +-                        -+   +-    -+
 |
 |         +-    -+   +-                 -+
 |         | d^2N |   | d^2x dN   d^2y dN |
 |         | ---- |   | ----*-- + ----*-- |
 |         | dr^2 |   | dr^2 dx   dr^2 dy |
 |         |      |   |                   |
 |         | d^2N |   | d^2x dN   d^2y dN |
 |      =  | ---- | - | ----*-- + ----*-- |
 |         | ds^2 |   | ds^2 dx   ds^2 dy |
 |         |      |   |                   |
 |         | d^2N |   | d^2x dN   d^2y dN |
 |         | ---- |   | ----*-- + ----*-- |
 |         | drds |   | drds dx   drds dy |
 |         +-    -+   +-                 -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
void gder2(
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,            DRT::UTILS::DisTypeToDim<DISTYPE>::dim>&                          xjm,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  derxy,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2,DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  deriv2,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,            DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  xyze,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2,      DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>&  derxy2
)
{
  // some numbers already known during compilation
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const int nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;
  const int numderiv2 = DRT::UTILS::DisTypeToNumDeriv2<DISTYPE>::numderiv2;

  // jacobian_bar matrix
  static LINALG::Matrix<numderiv2,numderiv2> bm;

  // calculate elements of jacobian_bar matrix
  switch (nsd)
  {
  case 3:
  {
    bm(0,0) = xjm(0,0)*xjm(0,0);
    bm(1,0) = xjm(1,0)*xjm(1,0);
    bm(2,0) = xjm(2,0)*xjm(2,0);
    bm(3,0) = xjm(0,0)*xjm(1,0);
    bm(4,0) = xjm(0,0)*xjm(2,0);
    bm(5,0) = xjm(2,0)*xjm(1,0);

    bm(0,1) = xjm(0,1)*xjm(0,1);
    bm(1,1) = xjm(1,1)*xjm(1,1);
    bm(2,1) = xjm(2,1)*xjm(2,1);
    bm(3,1) = xjm(0,1)*xjm(1,1);
    bm(4,1) = xjm(0,1)*xjm(2,1);
    bm(5,1) = xjm(2,1)*xjm(1,1);

    bm(0,2) = xjm(0,2)*xjm(0,2);
    bm(1,2) = xjm(1,2)*xjm(1,2);
    bm(2,2) = xjm(2,2)*xjm(2,2);
    bm(3,2) = xjm(0,2)*xjm(1,2);
    bm(4,2) = xjm(0,2)*xjm(2,2);
    bm(5,2) = xjm(2,2)*xjm(1,2);

    bm(0,3) = 2.*xjm(0,0)*xjm(0,1);
    bm(1,3) = 2.*xjm(1,0)*xjm(1,1);
    bm(2,3) = 2.*xjm(2,0)*xjm(2,1);
    bm(3,3) = xjm(0,0)*xjm(1,1)+xjm(1,0)*xjm(0,1);
    bm(4,3) = xjm(0,0)*xjm(2,1)+xjm(2,0)*xjm(0,1);
    bm(5,3) = xjm(1,0)*xjm(2,1)+xjm(2,0)*xjm(1,1);

    bm(0,4) = 2.*xjm(0,0)*xjm(0,2);
    bm(1,4) = 2.*xjm(1,0)*xjm(1,2);
    bm(2,4) = 2.*xjm(2,0)*xjm(2,2);
    bm(3,4) = xjm(0,0)*xjm(1,2)+xjm(1,0)*xjm(0,2);
    bm(4,4) = xjm(0,0)*xjm(2,2)+xjm(2,0)*xjm(0,2);
    bm(5,4) = xjm(1,0)*xjm(2,2)+xjm(2,0)*xjm(1,2);

    bm(0,5) = 2.*xjm(0,1)*xjm(0,2);
    bm(1,5) = 2.*xjm(1,1)*xjm(1,2);
    bm(2,5) = 2.*xjm(2,1)*xjm(2,2);
    bm(3,5) = xjm(0,1)*xjm(1,2)+xjm(1,1)*xjm(0,2);
    bm(4,5) = xjm(0,1)*xjm(2,2)+xjm(2,1)*xjm(0,2);
    bm(5,5) = xjm(1,1)*xjm(2,2)+xjm(2,1)*xjm(1,2);
  }
  break;
  case 2:
  {
    bm(0,0) =                   xjm(0,0)*xjm(0,0);
    bm(0,1) =                   xjm(0,1)*xjm(0,1);
    bm(0,2) =               2.0*xjm(0,0)*xjm(0,1);

    bm(1,0) =                   xjm(1,0)*xjm(1,0);
    bm(1,1) =                   xjm(1,1)*xjm(1,1);
    bm(1,2) =               2.0*xjm(1,1)*xjm(1,0);

    bm(2,0) =                   xjm(0,0)*xjm(1,0);
    bm(2,1) =                   xjm(0,1)*xjm(1,1);
    bm(2,2) = xjm(0,0)*xjm(1,1)+xjm(0,1)*xjm(1,0);
  }
  break;
  case 1:
  {
    bm(0,0) = xjm(0,0)*xjm(0,0);
  }
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd);
  } // switch nsd

  /*------------------ determine 2nd derivatives of coord.-functions */

  /*
  | 3 space dimensions:
  |
  |         0 1 2              0...iel-1
  |        +-+-+-+             +-+-+-+-+        0 1 2
  |        | | | | 0           | | | | | 0     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | 0
  |        | | | | 1           | | | | | 1   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 2           | | | | | 2     +-+-+-+
  |        +-+-+-+       =     +-+-+-+-+       | | | | .
  |        | | | | 3           | | | | | 3     +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 4           | | | | | 4   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 5           | | | | | 5     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | iel-1
  |                                            +-+-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-                  -+
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dr^2   dr^2   dr^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | ds^2   ds^2   ds^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dt^2   dt^2   dt^2 |
  |               yields    xder2  =    |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drds   drds   drds |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drdt   drdt   drdt |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dsdt   dsdt   dsdt |
  |                                     +-                  -+
  |
  |
  | 2 space dimensions:
  |                                             0 1
  |         0 1              0...iel-1         +-+-+
  |        +-+-+             +-+-+-+-+         | | | 0
  |        | | | 0           | | | | | 0       +-+-+
  |        +-+-+             +-+-+-+-+         | | | .
  |        | | | 1     =     | | | | | 1     * +-+-+ .
  |        +-+-+             +-+-+-+-+         | | | .
  |        | | | 2           | | | | | 2       +-+-+
  |        +-+-+             +-+-+-+-+         | | | iel-1
  |                                            +-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-           -+
  |                                     | d^2x   d^2y |
  |                                     | ----   ---- |
  |                                     | dr^2   dr^2 |
  |                                     |             |
  |                                     | d^2x   d^2y |
  |                 yields    xder2  =  | ----   ---- |
  |                                     | ds^2   ds^2 |
  |                                     |             |
  |                                     | d^2x   d^2y |
  |                                     | ----   ---- |
  |                                     | drds   drds |
  |                                     +-           -+
   */

  static LINALG::Matrix<numderiv2,nsd> xder2;
  xder2.MultiplyNT(deriv2,xyze);

  /*
  | 3 space dimensions:
  |
  |        0...iel-1             0 1 2
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 0          | | | | 0
  |        +-+-+-+-+            +-+-+-+            0...iel-1
  |        | | | | | 1          | | | | 1         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 0
  |        | | | | | 2          | | | | 2         +-+-+-+-+
  |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
  |        | | | | | 3          | | | | 3         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 2
  |        | | | | | 4          | | | | 4         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 5          | | | | 5          derxy
  |        +-+-+-+-+            +-+-+-+
  |
  |       chainrulerhs          xder2
  |
  |
  |
  | 2 space dimensions:
  |
  |        0...iel-1             0 1
  |        +-+-+-+-+            +-+-+               0...iel-1
  |        | | | | | 0          | | | 0             +-+-+-+-+
  |        +-+-+-+-+            +-+-+               | | | | | 0
  |        | | | | | 1     =    | | | 1     *       +-+-+-+-+   * (-1)
  |        +-+-+-+-+            +-+-+               | | | | | 1
  |        | | | | | 2          | | | 2             +-+-+-+-+
  |        +-+-+-+-+            +-+-+
  |
  |       chainrulerhs          xder2                 derxy
   */

  derxy2.Multiply(-1.0,xder2,derxy);

  /*
  | 3 space dimensions:
  |
  |        0...iel-1            0...iel-1         0...iel-1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 0          | | | | | 0       | | | | | 0
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 1          | | | | | 1       | | | | | 1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 2          | | | | | 2       | | | | | 2
  |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
  |        | | | | | 3          | | | | | 3       | | | | | 3
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 4          | | | | | 4       | | | | | 4
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 5          | | | | | 5       | | | | | 5
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |
  |       chainrulerhs         chainrulerhs        deriv2
  |
  |
  | 2 space dimensions:
  |
  |        0...iel-1             0...iel-1             0...iel-1
  |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
  |        | | | | | 0           | | | | | 0           | | | | | 0
  |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
  |        | | | | | 1     =     | | | | | 1     +     | | | | | 1
  |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
  |        | | | | | 2           | | | | | 2           | | | | | 2
  |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
  |
  |       chainrulerhs          chainrulerhs             deriv2
   */

  derxy2.Update(1.0,deriv2,1.0);

  /* make LR decomposition and solve system for all right hand sides
   * (i.e. the components of chainrulerhs)
  |
  | 3 space dimensions:
  |
  |          0  1  2  3  4  5         i        i
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 0     | | 0    | | 0
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 1     | | 1    | | 1
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 2     | | 2    | | 2
  |        +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
  |        |  |  |  |  |  |  | 3     | | 3    | | 3
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 4     | | 4    | | 4
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 5     | | 5    | | 5
  |        +--+--+--+--+--+--+       +-+      +-+
  |                                   |        |
  |                                   |        |
  |                                   derxy2[i]|
  |                                    |
  |                                    chainrulerhs[i]
  |
  |   yields
  |
  |                      0...iel-1
  |                      +-+-+-+-+
  |                      | | | | | 0 = drdr
  |                      +-+-+-+-+
  |                      | | | | | 1 = dsds
  |                      +-+-+-+-+
  |                      | | | | | 2 = dtdt
  |            derxy2 =  +-+-+-+-+
  |                      | | | | | 3 = drds
  |                      +-+-+-+-+
  |                      | | | | | 4 = drdt
  |                      +-+-+-+-+
  |                      | | | | | 5 = dsdt
  |                      +-+-+-+-+
  |
  |
  | 2 space dimensions:
  |
  |          0  1  2         i        i
  |        +--+--+--+       +-+      +-+
  |        |  |  |  | 0     | | 0    | | 0
  |        +--+--+--+       +-+      +-+
  |        |  |  |  | 1  *  | | 1 =  | | 1  for i=0...iel-1
  |        +--+--+--+       +-+      +-+
  |        |  |  |  | 2     | | 2    | | 2
  |        +--+--+--+       +-+      +-+
  |                          |        |
  |                          |        |
  |                        derxy2[i]  |
  |                                   |
  |                              chainrulerhs[i]
  |
  |
  |                   0...iel-1
  |                   +-+-+-+-+
  |                   | | | | | 0
  |                   +-+-+-+-+
  |        yields     | | | | | 1
  |                   +-+-+-+-+
  |                   | | | | | 2
  |                   +-+-+-+-+
  |
  |                    derxy2
   */

  LINALG::FixedSizeSerialDenseSolver<numderiv2,numderiv2,numnode> solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(derxy2,derxy2);
  solver.Solve();

  return;

} // DRT::UTILS::gder2


} // namespace UTILS
} // namespace DRT

#endif  // #ifndef UTILS_GDER2_H
#endif  // #ifdef CCADISCRET
