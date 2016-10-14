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
double PARTICLE::WeightFunction_CubicBspline::ComputeWeight(
  const double distRel,
  const double radius
  )
{
  // safety checks
  assert(distRel >= 0);
  assert(radius >= 0);

  const double norm_const = 1.5 * M_1_PI;
  const double norm_dist_rel = 2 * distRel / radius;

  double weight = 0;
  if (norm_dist_rel< 1)
    weight = norm_const * (2.0 / 3.0 - std::pow(norm_dist_rel,2) + 0.5 * std::pow(norm_dist_rel,3));
  else if (norm_dist_rel< 2)
    weight = norm_const * (1.0 / 6.0) * std::pow(2 - norm_dist_rel,3);

  return weight;
}

/*-----------------------------------------------------------------------------*
 | compute the gradient of the cubicBspline weight function  cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
LINALG::Matrix<3,1> PARTICLE::WeightFunction_CubicBspline::ComputeGradientWeight(
  const LINALG::Matrix<3,1>& rRel,
  const double radius
  )
{
  // safety checks
  assert(radius >= 0);

  const double norm_const = 1.5 * M_1_PI;
  const double rRelNorm = rRel.Norm2();
  const double resizer_temp = 2 / radius;
  const double norm_dist_rel = resizer_temp * rRelNorm;
  const double resizer = resizer_temp / rRelNorm;

  // create an empty vector
  LINALG::Matrix<3,1> weight;

  if (norm_dist_rel< 1)
    weight.Update(norm_const * resizer * (1.5 * std::pow(norm_dist_rel,2) - 2 * norm_dist_rel),rRel);
  else if (norm_dist_rel< 2)
    weight.Update(- 0.5 * norm_const * resizer * std::pow(norm_dist_rel - 2,2),rRel);

  return weight;
}

