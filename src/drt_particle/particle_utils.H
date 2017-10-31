/*----------------------------------------------------------------------*/
/*!
\file particle_utils.H

\brief General functions for the particle dynamics

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
      break;
    }
    }
    return effectiveDistance;
  }

  //! Compute the pressure from the proper speed of sound and the density (small density variation hypothesis)
  static double Density2Pressure(const double& speedOfSound, const double& initDensity, const double& refdensfac, const double& density, const double& exponent)
  {
    if(exponent==1)
    {
      return speedOfSound * speedOfSound * (density-refdensfac*initDensity);
    }
    else
    {
      const double restPressure=speedOfSound * speedOfSound * initDensity / exponent;
      return  restPressure*(pow(density/initDensity,exponent)-refdensfac);
    }
  };

  //! Compute the density from the proper speed of sound and the pressure (small density variation hypothesis)
  static double Pressure2Density(const double& speedOfSound, const double& initDensity, const double& refdensfac, const double& pressure, const double& exponent)
  {
    if(exponent==1)
      return pressure/(speedOfSound * speedOfSound) + refdensfac*initDensity;
    else
    {
      const double restPressure=speedOfSound * speedOfSound * initDensity / exponent;
      return initDensity*(pow(pressure/restPressure+refdensfac,(double)(1.0/exponent)));
    }
  };

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

  //! compute the intersection area of two particles that are in contact. It returns 0 if there is no contact
  static double IntersectionAreaPvsP(const double& radius1, const double& radius2, const double& dis);

};
}

#endif /* SRC_DRT_PARTICLE_PARTICLE_UTILS_H_ */
