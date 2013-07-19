/*----------------------------------------------------------------------*/
/*!
\file combust3_facemat_twophaseflow.H

\brief system matrix formulation for two-phase flow problem / 3D combustion element

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_FACEMAT_TWOPHASEFLOW_H
#define COMBUST3_FACEMAT_TWOPHASEFLOW_H

#include "combust3_ele.H"
#include "combust_time_integration_element.H"
#include "../drt_geometry/position_array.H"

namespace COMBUST
{

void FaceMatTwoPhase(
     const DRT::ELEMENTS::Combust3IntFace*      ele      ///< the element whose matrix is calculated
)
{
  // hier kommt keine Parameter Liste an, denn die werden uebergeben
  // state groessen fur master und slave werden ebenfalls uebergeben

  // time integration constant
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  return;
}

} // end namespace COMBUST
