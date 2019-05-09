/*---------------------------------------------------------------------------*/
/*!
\brief particle material thermo

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_material_thermo.H"

#include "../drt_mat/matpar_bundle.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialThermo::ParticleMaterialThermo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      initTemperature_(matdata->GetDouble("INITTEMPERATURE")),
      thermalCapacity_(matdata->GetDouble("THERMALCAPACITY")),
      invThermalCapacity_((thermalCapacity_ > 0.0) ? (1.0 / thermalCapacity_) : 0.0),
      thermalConductivity_(matdata->GetDouble("THERMALCONDUCTIVITY")),
      thermalAbsorptivity_(matdata->GetDouble("THERMALABSORPTIVITY"))
{
  // empty constructor
}
