/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material thermo

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_PARTICLE_THERMO_HPP
#define FOUR_C_MAT_PARTICLE_THERMO_HPP

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialThermo : virtual public CORE::MAT::PAR::Parameter
    {
     public:
      //! constructor
      ParticleMaterialThermo(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      //! @name material parameters
      //@{

      //! initial temperature
      const double initTemperature_;

      //! thermal capacity
      const double thermalCapacity_;

      //! inverse thermal capacity
      const double invThermalCapacity_;

      //! thermal conductivity
      const double thermalConductivity_;

      //! thermal absorptivity
      const double thermalAbsorptivity_;

      //@}

      //! create material instance of matching type with parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override = 0;
    };

  }  // namespace PAR

}  // namespace MAT

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
