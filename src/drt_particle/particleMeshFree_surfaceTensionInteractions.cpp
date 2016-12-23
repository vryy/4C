/*
 * particleMeshFree_surfaceTensionInteractions.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: cattabiani
 */


/*----------------------------------------------------------------------*/
/*!
 \file particleMeshFree_surfaceTensionInteractions.cpp

 \brief surface tension - related functions

theory: http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf R1

 \level 3

 \maintainer Alessandro Cattabiani
 */

/*----------------------------------------------------------------------*/
/* headers */
#include "particleMeshFree_surfaceTensionInteractions.H"

#include <math.h>

/*----------------------------------------------------------------------*
 | compute the cohesion coefficient                   cattabiani 12/16  |
 *----------------------------------------------------------------------*/

// The formula can be found in:
// R1-E2

double PARTICLE::SurfaceTensionInteractions::Cohesion(
  const double &disRel,
  const double &radius
  )
{

  double weight = 0;
  if (2 * disRel < radius)
  {
    weight = 2 * std::pow(radius - disRel,3) * std::pow(disRel,3) - std::pow(radius,6)/64;
  }
  else if (disRel < radius)
  {
    weight = std::pow(radius - disRel,3) * std::pow(disRel,3);
  }

  // resizing
  weight *= 32 * M_1_PI / std::pow(radius,9);

  return weight;
}


/*----------------------------------------------------------------------*
 | compute the adhesion coefficient                   cattabiani 12/16  |
 *----------------------------------------------------------------------*/

// The formula can be found in:
// R1-E7

double PARTICLE::SurfaceTensionInteractions::Adhesion(
  const double &disRel,
  const double &radius
  )
{
  double weight = 0;

  if (disRel < radius && disRel > radius/2)
  {
    weight = 0.007 * std::pow(- 4 * disRel * disRel / radius + 6 * disRel - 2 * radius,0.25) / std::pow(radius,3.25);
  }

  return weight;
}


