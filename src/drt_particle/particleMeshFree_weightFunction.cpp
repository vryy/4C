/*----------------------------------------------------------------------*/
/*!
 \file particleMeshFree_weightFunction.cpp

 \brief weight functions for MeshFree methods

 \level 3

 \maintainer Alessandro Cattabiani
 */

/*----------------------------------------------------------------------*/
/* headers */
#include "particleMeshFree_weightFunction.H"

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------*
 | compute the cubicBspline weight function           cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_CubicBspline::Weight(
  const double &distRel,
  const double &radius
  )
{
  // safety checks
  assert(distRel >= 0);
  assert(radius >= 0);

  const double reszDis = 2 * distRel / radius;

  double weight = 0;
  if (reszDis < 1)
  {
    weight = (1.0 / 6.0) * std::pow(2-reszDis,3) - 4 * std::pow(1-reszDis,3);
  }
  else if (reszDis < 2)
  {
    weight = (1.0 / 6.0) * std::pow(2-reszDis,3);
  }

  // resizing to have an integral = 1
  weight *= Resizer3D();

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the weight function derivative                    cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_CubicBspline::DerivativeWeight(const double &distRel, const double &radius)
{
  // safety checks
  assert(distRel >= 0);
  assert(radius >= 0);

  const double derivativeResizer = 2 / radius;
  const double reszDis = derivativeResizer * distRel;

  double derivativeWeight = 0;
  if (reszDis < 1)
  {
    derivativeWeight = 0.5 * reszDis * (3 * reszDis - 4);
  }
  else if (reszDis < 2)
  {
    derivativeWeight = - 0.5 * (reszDis - 2) * (reszDis - 2);
  }

  // resizing to have an integral = 1
  derivativeWeight *= Resizer3D() * derivativeResizer / distRel;

  return derivativeWeight;
}


/*-----------------------------------------------------------------------------*
 | compute the gradient of the cubicBspline weight function  cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
/*
LINALG::Matrix<3,1> PARTICLE::WeightFunction_CubicBspline::GradientWeight(LINALG::Matrix<3,1> &rRel, const double &radius)
{
  // safety checks
  assert(radius > 0);

  const double distRel = rRel.Norm2();

  // solving the particular case in which two particles perfectly overlap
  if (distRel <= 1e-15)
  {
    dserror("Warning! particles are overlapping! Right now, it is not allowed");
    std::cout << "Warning! particles are overlapping!\n";
    return rRel;
  }

  LINALG::Matrix<3,1> WFGrad(rRel);
  const double derivativeWeight = DerivativeWeight(distRel, radius);
  WFGrad.Scale(derivativeWeight);

  return WFGrad;
}
*/




