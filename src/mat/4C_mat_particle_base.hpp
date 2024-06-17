/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material base

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_PARTICLE_BASE_HPP
#define FOUR_C_MAT_PARTICLE_BASE_HPP

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
namespace Mat
{
  namespace PAR
  {
    class ParticleMaterialBase : virtual public Core::Mat::PAR::Parameter
    {
     public:
      //! constructor
      ParticleMaterialBase(const Core::Mat::PAR::Parameter::Data& matdata);

      //! @name material parameters
      //@{

      //! initial radius
      const double initRadius_;

      //! initial density
      const double initDensity_;

      //@}

      //! create material instance of matching type with parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override = 0;
    };

  }  // namespace PAR

}  // namespace Mat

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
