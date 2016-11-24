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

// The formula can be found in:
// Numerical simulation of fluid-structure interaction by SPH, DOI: 10.1016/j.compstruc.2007.01.002

double PARTICLE::WeightFunction_CubicBspline::Weight(
  const double &disRel,
  const double &radius
  )
{


  const double rszDisRel = RszDisRel(disRel,radius);

  double weight = 0;
  if (rszDisRel < 1)
  {
    weight = 1 - 1.5 * std::pow(rszDisRel,2) - 0.75 * std::pow(rszDisRel,3);
  }
  else if (rszDisRel < 2)
  {
    weight = 0.25 * std::pow(2-rszDisRel,3);
  }

  // resizing to have an integral = 1
  weight *= Rsz3D(radius);

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the weight function derivative                    cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D

double PARTICLE::WeightFunction_CubicBspline::WeightDerivative(const double &disRel, const double &radius)
{

  const double rszDisRel = RszDisRel(disRel,radius);

  double weightDerivative = 0;
  if (rszDisRel < 1)
  {
    weightDerivative = 6.0 * disRel * (3.0 * disRel - 2.0 * radius) / std::pow(radius,3);
  }
  else if (rszDisRel < 2)
  {
    weightDerivative = - 6.0 * std::pow(radius - disRel,2) / std::pow(radius,3);
  }

  // resizing to have an integral = 1
  // it is divided by disRel so that DW * rRel is the gradient
  weightDerivative *= Rsz3D(radius);

  return weightDerivative;
}


/*----------------------------------------------------------------------*
 | compute the cubicBspline weight function           cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::Weight(
  const double &disRel,
  const double &radius
  )
{
  double weight = 0;
  if (disRel<radius)
  {
    weight = (std::pow(radius/disRel,0.5) - 1) * Rsz3D(radius);
  }

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the weight function derivative                    cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::WeightDerivative(const double &disRel, const double &radius)
{
  double weightDerivative = 0;
  if (disRel<radius)
  {
    weightDerivative = (- std::pow(radius/disRel,0.5) / (2 * disRel) ) * Rsz3D(radius);
  }

  return weightDerivative;
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




