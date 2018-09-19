/*----------------------------------------------------------------------*/
/*!
\file particle_utils.H

\brief General functions for the particle dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_UTILS_H_
#define PARTICLE_UTILS_H_
/*----------------------------------------------------------------------*/
/* headers */
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <math.h>
#include "../linalg/linalg_utils.H"
#include <boost/unordered_map.hpp>
#include "particle_algorithm.H"

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

//! forward declarations
namespace MAT
{
  namespace PAR
  {
    class ParticleMat;
    class ExtParticleMat;
  }  // namespace PAR
}  // namespace MAT

namespace PARTICLE
{
  class Utils
  {
   public:
    //! Radius -> volume
    static double Radius2Volume(double radius)
    {
      return (4.0 / 3.0) * M_PI * radius * radius * radius;
    }

    //! Volume -> radius
    static double Volume2Radius(double volume)
    {
      return std::pow((3.0 / 4.0) * M_1_PI * volume, 1.0 / 3.0);
    }

    //! Volume -> effective particle distance
    static double Volume2EffDist(double volume, const INPAR::PARTICLEOLD::WeightFunctionDim WF_DIM)
    {
      // This method is used to determin the effective particle distance based on the volume:
      // dist=(volumme)^(1/WF_DIM) The following cumbersome case destinction is made in order to
      // avoid the use of the expensive std::pow() function:
      double effectiveDistance = 0.0;
      switch (WF_DIM)
      {
        case INPAR::PARTICLEOLD::WF_3D:
        {
          // cubic root
          effectiveDistance = cbrt(volume);
          break;
        }
        case INPAR::PARTICLEOLD::WF_2D:
        {
          // square root
          effectiveDistance = sqrt(volume);
          break;
        }
        case INPAR::PARTICLEOLD::WF_1D:
        {
          // 1D volume = length
          effectiveDistance = volume;
          break;
        }
        default:
        {
          dserror("Only the weight function / problem dimensions 1, 2 and 3 are possible!");
          break;
        }
      }
      return effectiveDistance;
    }

    //! compute the intersection area of two particles that are in contact. It returns 0 if there is
    //! no contact
    static double IntersectionAreaPvsP(
        const double& radius1, const double& radius2, const double& dis);
  };
}  // namespace PARTICLE

#endif /* PARTICLE_UTILS_H_ */
