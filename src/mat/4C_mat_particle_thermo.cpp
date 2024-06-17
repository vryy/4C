/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material thermo

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_thermo.hpp"

#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialThermo::ParticleMaterialThermo(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      initTemperature_(matdata.parameters.Get<double>("INITTEMPERATURE")),
      thermalCapacity_(matdata.parameters.Get<double>("THERMALCAPACITY")),
      invThermalCapacity_((thermalCapacity_ > 0.0) ? (1.0 / thermalCapacity_) : 0.0),
      thermalConductivity_(matdata.parameters.Get<double>("THERMALCONDUCTIVITY")),
      thermalAbsorptivity_(matdata.parameters.Get<double>("THERMALABSORPTIVITY"))
{
  // empty constructor
}

FOUR_C_NAMESPACE_CLOSE
