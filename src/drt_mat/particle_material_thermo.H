/*---------------------------------------------------------------------------*/
/*!
\brief particle material thermo

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_MATERIAL_THERMO_H
#define PARTICLE_MATERIAL_THERMO_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

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

      //! virtual desctructor
      virtual ~ParticleMaterialThermo(){};

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() = 0;
    };

  }  // namespace PAR

}  // namespace MAT

/*---------------------------------------------------------------------------*/
#endif
