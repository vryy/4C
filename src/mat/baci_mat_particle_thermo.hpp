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
#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_par_parameter.hpp"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialThermo : virtual public Parameter
    {
     public:
      //! constructor
      ParticleMaterialThermo(Teuchos::RCP<MAT::PAR::Material> matdata);

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
      Teuchos::RCP<MAT::Material> CreateMaterial() override = 0;
    };

  }  // namespace PAR

}  // namespace MAT

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
