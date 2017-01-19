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
 | gradient, r_ij = r_i - r_j                         cattabiani 08/16  |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> PARTICLE::WeightFunction_Base::GradW(const LINALG::Matrix<3,1> &rVersor, const double& dw)
{
  LINALG::Matrix<3,1> gradW(rVersor);
  gradW.Scale(dw);

  return gradW;
}


/*----------------------------------------------------------------------*
 | hessian, r_ij = r_i - r_j                          cattabiani 08/16  |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,3> PARTICLE::WeightFunction_Base::HessW(const LINALG::Matrix<3,1> &rVersor, const double& dw, const double& rNorm2, const double& ddw)
{
  // compute the resized first derivative
  const double dw_r = dw / rNorm2;
  // diadic product
  LINALG::Matrix<1,3> rVersorT;
  rVersorT.UpdateT(rVersor);
  LINALG::Matrix<3,3> rVersor2;
  rVersor2.MultiplyNN(rVersor, rVersorT);
  // create the hessian and add the various terms
  LINALG::Matrix<3,3> hessW;
  for (int ii = 0; ii<3; ++ii)
  {
    hessW(ii,ii) = dw_r;
  }
  hessW.Update(ddw - dw_r, rVersor2, 1.0);

  return hessW;
}

/*----------------------------------------------------------------------*
 | compute the cubicBspline w function           cattabiani 08/16  |
 *----------------------------------------------------------------------*/

// The 3D variant can be found in Monaghan2005, Eq. (2.6)
// [see also Eq. (9) in Antoci2007]

double PARTICLE::WeightFunction_CubicBspline::W(
  const double &disRel,
  const double &radius
  )
{

  //Attention: The support of our kernel functions is (in 3D) defined by a sphere with radius measured by our variable "radius".
  //In the SLM literature, typically the smoothing length h:=radius/2 is used for defining the kernel functions!!!
  //rszDisRel=q:=||r_{ij}||/h=2*||r_{ij}||/radius
  const double rszDisRel = RszDisRel(disRel,radius);

  double w = 0;
  if (rszDisRel < 1)
  {
    w = (2.0/3.0) - std::pow(rszDisRel,2) + 0.5 * std::pow(rszDisRel,3);
  }
  else if (rszDisRel < 2)
  {
    w = std::pow(2-rszDisRel,3)/6;
  }

  // resizing to have an integral = 1
  w *= RszDim(radius);

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the cubicBspline w function derivative       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

double PARTICLE::WeightFunction_CubicBspline::DW(const double &disRel, const double &radius)
{

  // In general, we have: grad(W) = dW/dr_{ij} = dW/dq * dq/d||r_{ij}|| * d||r_{ij}||/dr_{ij} = dW/dq * 1/h * r_{ij}/||r_{ij}||
  // Here, we determine DW = dW/dq * 1/h while the part e_{ij}:=r_{ij}/||r_{ij}|| is multiplied outside!
  const double rszDisRel = RszDisRel(disRel,radius);

  double dw = 0;
  if (rszDisRel < 1)
  {
    dw = (- 4 * rszDisRel + 3 * std::pow(rszDisRel,2) ) / radius;
  }
  else if (rszDisRel < 2)
  {
    dw = - std::pow(2 - rszDisRel,2) / radius;
  }

  // resizing to have an integral = 1
  dw *= RszDim(radius);

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the cubicBspline w function derivative       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

// empowered by mathematica:
// https://www.wolframalpha.com/input/?i=Piecewise+%5B%7B%7B+(6+x+(3+x+-+4))%2F8,+0%3C%3Dx%2F2%3C%3D1%2F2%7D,%7B+-(6+(2+-+x)%5E2)%2F8,1%2F2%3Cx%2F2%3C%3D1%7D%7D%5D
// See again Monaghan2005, Eq. (2.6)

double PARTICLE::WeightFunction_CubicBspline::DDW(const double &disRel, const double &radius)
{
  const double rszDisRel = RszDisRel(disRel,radius);

  double ddw = 0;
  if (rszDisRel < 1)
  {
    ddw = (- 8 + 12 * rszDisRel) / std::pow(radius,2);
  }
  else if (rszDisRel < 2)
  {
    ddw = 4 * (2 - rszDisRel) / std::pow(radius,2);
  }

  // resizing to have an integral = 1
  ddw *= RszDim(radius);

  return ddw;
}


/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_CubicBspline::RszDim(const double &radius)
{
  switch (PARTICLE_DIM)
  {
  case 3 :
    return 12.0 * M_1_PI / std::pow(radius,3);
  case 2 :
    return 60.0 * M_1_PI / (7.0*std::pow(radius,2));
  case 1 :
    return 2.0 / radius;
  default :
    dserror("Only the problem dimensions 1, 2 and 3 are possible!");
  }
}



/*----------------------------------------------------------------------*
 | compute the SqrtHyperbola w function          cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::W(
  const double &disRel,
  const double &radius
  )
{
  double w = 0;
  if (disRel<radius)
  {
    w = (std::pow(radius/disRel,0.5) - 1) * RszDim(radius);
  }

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 1 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::DW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel<radius)
  {
    dw = (- radius / (2 * std::pow(radius/disRel,0.5) * disRel)) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 2 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_SqrtHyperbola::DDW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel<radius)
  {
    dw = (- std::pow(radius,2) / (4 * std::pow(radius/disRel,3.0/2.0) * std::pow(disRel,4)) +
        radius / (std::pow(radius/disRel,0.5) * std::pow(disRel,3))) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | rsz in case of different dimensions                       cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/

double PARTICLE::WeightFunction_SqrtHyperbola::RszDim(const double &radius)
{
  double rszDim;
  switch (PARTICLE_DIM)
  {
  case 3 :
  {
    rszDim = 15.0 * M_1_PI / (4.0*std::pow(radius,3));
    break;
  }
  case 2 :
  {
    rszDim =  3.0 * M_1_PI / (std::pow(radius,2));
    break;
  }
  case 1 :
  {
    rszDim = 0.5 / radius;
    break;
  }
  default :
  {
    dserror("Only the problem dimensions 1, 2 and 3 are possible!");
    break;
  }
  }

  return rszDim;
}



/*----------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function         cattabiani 08/16  |
 *----------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::W(
  const double &disRel,
  const double &radius
  )
{
  double w = 0;
  if (disRel<radius)
  {
    w = (radius/disRel - 1) * RszDim(radius);
  }

  return w;
}


/*-----------------------------------------------------------------------------*
 | compute the HyperbolaNoRsz w function derivative     cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::DW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel<radius)
  {
    dw = - radius/(disRel * disRel) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the SqrtHyperbola 2 function derivative      cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
double PARTICLE::WeightFunction_HyperbolaNoRsz::DDW(const double &disRel, const double &radius)
{
  double dw = 0;
  if (disRel<radius)
  {
    dw = 2 * radius/std::pow(disRel, 3) * RszDim(radius);
  }

  return dw;
}


/*-----------------------------------------------------------------------------*
 | compute the gradient of the cubicBspline w function  cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
/*
LINALG::Matrix<3,1> PARTICLE::WeightFunction_CubicBspline::Gradientw(LINALG::Matrix<3,1> &rRel, const double &radius)
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
  const double derivativew = Derivativew(distRel, radius);
  WFGrad.Scale(derivativew);

  return WFGrad;
}
*/




