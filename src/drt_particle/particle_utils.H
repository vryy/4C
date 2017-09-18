/*----------------------------------------------------------------------*/
/*!
\file particle_utils.H

\brief General functions for the Particle-MeshFree dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
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

    //! Volume -> effective particle distance
   static double Volume2EffDist(double volume, const INPAR::PARTICLE::WeightFunctionDim WF_DIM)
   {
     //This method is used to determin the effective particle distance based on the volume: dist=(volumme)^(1/WF_DIM)
     //The following cumbersome case destinction is made in order to avoid the use of the expensive std::pow() function:
     double effectiveDistance=0.0;
     switch (WF_DIM)
     {
       case INPAR::PARTICLE::WF_3D :
       {
         //cubic root
         effectiveDistance=cbrt(volume);
         break;
       }
       case INPAR::PARTICLE::WF_2D :
       {
         //square root
         effectiveDistance=sqrt(volume);
         break;
       }
       case INPAR::PARTICLE::WF_1D :
       {
         //1D volume = length
         effectiveDistance=volume;
         break;
       }
       default :
       {
         dserror("Only the weight function / problem dimensions 1, 2 and 3 are possible!");
       }
     }
     return effectiveDistance;
   }

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

    //! Compute the pressure from the proper speed of sound and the density (small density variation hypothesis)
    /// if you change this function you need to touch Inter_pvp_gradAccP too!
    static double Density2Pressure(const double& speedOfSound, const double& restDensity, const double& refdensfac, const double& density, const double& exponent)
      {
        if(exponent==1)
        {
          return speedOfSound * speedOfSound * (density-refdensfac*restDensity);
        }
        else
        {
          const double restPressure=speedOfSound * speedOfSound * restDensity / exponent;
          return  restPressure*(pow(density/restDensity,exponent)-refdensfac);
        }
      };

    //! Compute the pressure from the proper speed of sound and the density (small density variation hypothesis)
    /// if you change this function you need to touch Inter_pvp_gradAccP too! (FAD version)
    static FAD Density2Pressure(const double& speedOfSound, const double& restDensity, const double& refdensfac, const FAD& density, const double& exponent)
      {
        if(exponent==1)
          return speedOfSound * speedOfSound * (density-refdensfac*restDensity);
        else
        {
          const double restPressure=speedOfSound * speedOfSound * restDensity / exponent;
          return  restPressure*(pow(density/restDensity,exponent)-refdensfac);
        }
      };

    //! Compute the density from the proper speed of sound and the pressure (small density variation hypothesis)
    static double Pressure2Density(const double& speedOfSound, const double& restDensity, const double& refdensfac, const double& pressure, const double& exponent)
    {
      if(exponent==1)
        return pressure/(speedOfSound * speedOfSound) + refdensfac*restDensity;
      else
      {
        const double restPressure=speedOfSound * speedOfSound * restDensity / exponent;
        return restDensity*(pow(pressure/restPressure+refdensfac,(double)(1.0/exponent)));
      }
    };

    //! Compute pressure for the entire vector
    /// if you change this function you need to touch Inter_pvp_gradAccP too!
    static void Density2Pressure(
        const double& restDensity,
        const double& refdensfac,
        const Teuchos::RCP<const Epetra_Vector> density,
        const Teuchos::RCP<const Epetra_Vector> specEnthalpy,
        Teuchos::RCP<Epetra_Vector> &pressure,
        const MAT::PAR::ExtParticleMat* extParticleMat,
        bool trg_createPressureVector = false,
        bool solve_thermal_problem = true);

    //! Compute the mechanical energy from the current density and the equation of state:
    //! This method has to be in accordance with the methodes Density2Pressure and Pressure2Density above
    static double Density2Energy(const double& speedOfSound, const double& density, const double& restDensity, const double& refdensfac, const double& mass)
    {
      // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
      // Attention: currently, only the first, pressure-dependent contribution of the thermodynamic energy is implemented!
      // Thus, it is only valid for isentrop problems, i.e.dE/dS=0! Furthermore, it is only considered for the fluid phase so far (since SpeedOfSoundL is used)!
      // From the considered pressure law p_i=c_0^2(rho_i+refdensfac_*restDensity_) and the relation rho_i=m_i/V_i it is possible to gain the energy via integration:
      // E_i=c_0^2 m_i [ln(rho_i/m_i)-refdensfac_*restDensity_/rho_i]+const.  --> The integration constant const. can be determined based on an arbitrary initial condition.
      // In the following, we choose E_i(rho_i=restDensity_)=0 leading to: E_i=c_0^2 m_i [ln(rho_i/restDensity_)+refdensfac_*{1-restDensity_/rho_i}]:

      double density_frac = density/restDensity;
      return mass * speedOfSound * speedOfSound * (log(density_frac)-refdensfac*(1.0-1.0/density_frac));
    };

    static double SpeedOfSound(
        const double& specEnthalpy,
        const MAT::PAR::ExtParticleMat* extParticleMat);

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
    static void Assemble(Epetra_MultiVector& V, matrix& Mele, const std::vector<int>& lm, const double& lmowner, const double& myrank)
    {
      if (lmowner != myrank)
      {
        return;
      }

      const std::vector<int> lmownerV(3,lmowner);

      // check
      if (Mele.N() != (unsigned)V.NumVectors() || lm.size() != Mele.M())
      {
        dserror("Received matrix of wrong size!");
      }
      if (lm.size() != lmownerV.size())
      {
        dserror("lm and lmowner sizes missmatch!");
      }

      // build up the assembler
      for (size_t ii = 0; ii < Mele.N(); ++ii)
      {
        const Epetra_SerialDenseVector Vele(::View, &(Mele(0,ii)), Mele.M());
        LINALG::Assemble(V, ii, Vele, lm, lmownerV);
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

    //! node * dof multiplication vector (usually to multiply the mass for the acceleration)
    //static void MultiplyND(const double coeff, const Epetra_Vector& nodeSV, Epetra_MultiVector& dofSV);

  };
}

#endif /* SRC_DRT_PARTICLE_PARTICLE_UTILS_H_ */
