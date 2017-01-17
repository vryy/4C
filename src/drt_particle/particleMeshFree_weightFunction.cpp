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

// The 3D variant can be found in Monaghan2005, Eq. (2.6)
// [see also Eq. (9) in Antoci2007]

double PARTICLE::WeightFunction_CubicBspline::Weight(
  const double &disRel,
  const double &radius
  )
{

  //Attention: The support of our kernel functions is (in 3D) defined by a sphere with radius measured by our variable "radius".
  //In the SLM literature, typically the smoothing length h:=radius/2 is used for defining the kernel functions!!!
  //rszDisRel=q:=||r_{ij}||/h=2*||r_{ij}||/radius
  const double rszDisRel = RszDisRel(disRel,radius);

  double weight = 0;
  if (rszDisRel < 1)
  {
    weight = 1 - 1.5 * std::pow(rszDisRel,2) + 0.75 * std::pow(rszDisRel,3);
  }
  else if (rszDisRel < 2)
  {
    weight = 0.25 * std::pow(2-rszDisRel,3);
  }

  // resizing to have an integral = 1
  weight *= RszDim(radius);

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the cubicBspline weight function derivative       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

double PARTICLE::WeightFunction_CubicBspline::WeightDerivative(const double &disRel, const double &radius)
{

  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} = dW/dq * 1/h * r_{ij}/||r_{ij}||
  // Here, we determine weightDerivative = dW/dq * 1/h while the part e_{ij}:=r_{ij}/||r_{ij}|| is multiplied outside!
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
  weightDerivative *= RszDim(radius);

  return weightDerivative;
}

/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_CubicBspline::RszDim(const double &radius)
{
  switch (PARTICLE_DIM)
  {
  case 3 :
    return 8.0 * M_1_PI / std::pow(radius,3);
  case 2 :
    return 40.0 * M_1_PI / (7.0*std::pow(radius,2));
  case 1 :
    return 4.0 / (3.0*radius);
  default :
    dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }
}



/*----------------------------------------------------------------------*
 | compute the SqrtHyperbola weight function          cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::Weight(
  const double &disRel,
  const double &radius
  )
{
  double weight = 0;
  if (disRel<radius)
  {
    weight = (std::pow(radius/disRel,0.5) - 1) * RszDim(radius);
  }

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola weight function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::WeightDerivative(const double &disRel, const double &radius)
{
  double weightDerivative = 0;
  if (disRel<radius)
  {
    weightDerivative = (- std::pow(radius/disRel,0.5) / (2 * disRel) ) * RszDim(radius);
  }

  return weightDerivative;
}


/*----------------------------------------------------------------------*
 | compute the HyperbolaNoRsz weight function         cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::Weight(
  const double &disRel,
  const double &radius
  )
{
  double weight = 0;
  if (disRel<radius)
  {
    weight = radius/disRel - 1;
  }

  return weight;
}


/*-----------------------------------------------------------------------------*
 | compute the HyperbolaNoRsz weight function derivative     cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::WeightDerivative(const double &disRel, const double &radius)
{
  double weightDerivative = 0;
  if (disRel<radius)
  {
    weightDerivative = - radius/(disRel * disRel);
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




