/*----------------------------------------------------------------------*/
/*!
\file particle_utils.cpp

\brief General functions for the particle dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#include "particle_utils.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/extparticle_mat.H"

/*-----------------------------------------------------------------------------*/
// compute the intersection area of two particles that are in contact. It returns 0 if there is no contact
// TODO: (see reference... equation...)
double PARTICLE::Utils::IntersectionAreaPvsP(const double& radius1, const double& radius2, const double& dis)
{
  //checks
  if (radius1<=0 || radius2<=0 || dis<=0)
    dserror("input parameters are unacceptable");
  if (dis >= radius1 + radius2)
    return 0;

  return M_PI * (- 0.25 * dis * dis
                 - radius1 * radius1 * radius1 * radius1 / (4 * dis * dis)
                 - radius2 * radius2 * radius2 * radius2 / (4 * dis * dis)
                 + 0.5 * ( radius1 * radius1 + radius2 * radius2 + radius1 * radius1 * radius2 * radius2));

}
