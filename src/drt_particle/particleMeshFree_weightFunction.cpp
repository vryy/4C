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
  const double dist_rel,
  const double invRadius
  )
{
  // assertions for the sake of safety
  assert(dist_rel >= 0);
  assert(invRadius >= 0);

  const double norm_const = 3 / (2 * M_PI);
  const double norm_dist_rel = 2 * dist_rel * invRadius;

  double weight = 0;
  if (norm_dist_rel< 1)
  {
    weight = norm_const
        * (2 / 3 + norm_dist_rel * norm_dist_rel
            + 0.5 * norm_dist_rel * norm_dist_rel * norm_dist_rel);
  }
  else if (norm_dist_rel< 2)
  {
    weight = 2 - norm_dist_rel;
    weight = norm_const * (1 / 6) * weight * weight * weight;
  }

  return weight;
}

/*-----------------------------------------------------------------------------*
 | compute the gradient of the cubicBspline weight function  cattabiani 08/16  |
 *-----------------------------------------------------------------------------*/
LINALG::Matrix<3,1> PARTICLE::WeightFunction_CubicBspline::ComputeGradientWeight(
  const LINALG::Matrix<3,1>& r_rel,
  const double invRadius
  )
{

  assert(invRadius >= 0);

  const double norm_const = 3 / (2 * M_PI);
  const double resizer = 2 * invRadius;

  const double norm_dist_rel = resizer * r_rel.Norm2();

  LINALG::Matrix<3,1> weight;

  if (norm_dist_rel< 1)
    weight.Update(norm_const * resizer * ((3 / 2) * norm_dist_rel - 2),r_rel);
  else if (norm_dist_rel< 2)
    weight.Update(- norm_const * resizer * (0.5 * (norm_dist_rel - 2) * (norm_dist_rel - 2) / norm_dist_rel),r_rel);

  return weight;
}

