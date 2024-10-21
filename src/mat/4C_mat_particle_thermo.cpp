#include "4C_mat_particle_thermo.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialThermo::ParticleMaterialThermo(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      initTemperature_(matdata.parameters.get<double>("INITTEMPERATURE")),
      thermalCapacity_(matdata.parameters.get<double>("THERMALCAPACITY")),
      invThermalCapacity_((thermalCapacity_ > 0.0) ? (1.0 / thermalCapacity_) : 0.0),
      thermalConductivity_(matdata.parameters.get<double>("THERMALCONDUCTIVITY")),
      thermalAbsorptivity_(matdata.parameters.get<double>("THERMALABSORPTIVITY"))
{
  // empty constructor
}

FOUR_C_NAMESPACE_CLOSE
