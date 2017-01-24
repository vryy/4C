/*----------------------------------------------------------------------*/
/*!
\file particle_utils.H

\brief General functions for the particle and meshFree dynamics

\level 3

\maintainer Alessandro Cattabiani
*/

/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLE_UTILS_H_
#define SRC_DRT_PARTICLE_PARTICLE_UTILS_H_
/*----------------------------------------------------------------------*/
/* headers */
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <math.h>
#include "../linalg/linalg_utils.H"

//! forward declarations
namespace MAT
{
  namespace PAR
  {
    class ParticleMat;
    class ExtParticleMat;
  }
}

namespace PARTICLE
{
  class Utils
  {
  public :
    //! Radius -> volume
    static double Radius2Volume(double radius) { return (4.0 / 3.0) * M_PI * radius * radius * radius; }

    //! Volume -> radius
    static double Volume2Radius(double volume) { return std::pow((3.0 / 4.0) * M_1_PI * volume, 1.0 / 3.0); }

    //! inertia, sphere-like particles (I = 2/5 * m * r^2)
    static double ComputeInertia(const double &radius, const double &mass) { return 0.4 * mass * radius * radius; }

    //! Compute the inertia vector
    //! The most updated radius is used (radius_ vs radiusn_)
    //! Inertia-vector -> sphere: I = 2/5 * m * r^2
    static void ComputeInertia(
      const Teuchos::RCP<const Epetra_Vector> radius,
      const Teuchos::RCP<const Epetra_Vector> mass,
      Teuchos::RCP<Epetra_Vector> &inertia,
      bool trg_createInertiaVector = false);

    //! Compute the pressure from the proper speed of sound and the density (small density variation hypotesis)
    static double Density2Pressure(const double& speedOfSound, const double& deltaDensity)
      { return speedOfSound * speedOfSound * deltaDensity; };

    //! Compute pressure for the entire vector
    static void Density2Pressure(
        const Teuchos::RCP<const Epetra_Vector> deltaDensity,
        const Teuchos::RCP<const Epetra_Vector> specEnthalpy,
        Teuchos::RCP<Epetra_Vector> &pressure,
        const MAT::PAR::ExtParticleMat* extParticleMat,
        bool trg_createPressureVector = false);

    //! Compute the temperature from the specEnthalpy
    static double SpecEnthalpy2Temperature(
        const double& specEnthalpy,
        const MAT::PAR::ExtParticleMat* extParticleMat);

    //! Compute the temperature from the specEnthalpy (vector version)
    static void SpecEnthalpy2Temperature(
        Teuchos::RCP<Epetra_Vector> temperature,
        const Teuchos::RCP<const Epetra_Vector>& specEnthalpy,
        const MAT::PAR::ExtParticleMat*);

    //! compute the intersection area of two particles that are in contact. It returns 0 if there is no contact
    static double IntersectionAreaPvsP(const double& radius1, const double& radius2, const double& dis);

    //! compute the Courant-Friedrich-Levy condition - meshfree dynamics - http://dx.doi.org/10.1145/2786784.2786796
    //! section 3.1 - inline equation (no formula number)
    static double CFLcondition(const double& radiusMin, const double& velocityMax) { return 0.8 * radiusMin / velocityMax; }

    //! small wrapping template to store matrices into MultiVector
    template <class matrix>
    static void Assemble(Epetra_MultiVector& V, matrix& Mele, const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      // check
      if (Mele.N() != (unsigned)V.NumVectors() || lm.size() != Mele.M())
      {
        dserror("Received matrix of wrong size!");
      }
      if (lm.size() != lmowner.size())
      {
        dserror("lm and lmowner sizes missmatch!");
      }

      // build up the assembler
      for (size_t ii = 0; ii < Mele.N(); ++ii)
      {
        const Epetra_SerialDenseVector Vele(::View, &(Mele(0,ii)), Mele.M());
        LINALG::Assemble(V, ii, Vele, lm, lmowner);
      }

      return;
    }

    //! small wrapping template to read matrices from MultiVector
    template <class matrix>
    static void ExtractMyValues(const Epetra_MultiVector& global, matrix& local, const std::vector<int>& lm)
    {
      /// safety checks
      if (local.N() != (unsigned)global.NumVectors() || lm.size() != local.M())
      {
        dserror("Received matrix of wrong size!");
      }
      // loop over the rows
      for (size_t ii=0; ii<local.M(); ++ii)
      {
        const int lid = global.Map().LID(lm[ii]);
        // loop over the columns
        for (size_t jj=0; jj<local.N(); ++jj)
        {
          local(ii,jj) = (global[jj])[lid];
        }
      }
    }
  };
}

#endif /* SRC_DRT_PARTICLE_PARTICLE_UTILS_H_ */
