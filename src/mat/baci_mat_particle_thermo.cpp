/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material thermo

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_mat_particle_thermo.hpp"

#include "baci_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialThermo::ParticleMaterialThermo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      initTemperature_(*matdata->Get<double>("INITTEMPERATURE")),
      thermalCapacity_(*matdata->Get<double>("THERMALCAPACITY")),
      invThermalCapacity_((thermalCapacity_ > 0.0) ? (1.0 / thermalCapacity_) : 0.0),
      thermalConductivity_(*matdata->Get<double>("THERMALCONDUCTIVITY")),
      thermalAbsorptivity_(*matdata->Get<double>("THERMALABSORPTIVITY"))
{
  // empty constructor
}

FOUR_C_NAMESPACE_CLOSE
