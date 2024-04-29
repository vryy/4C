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
#include "4C_mat_material.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialBase : virtual public Parameter
    {
     public:
      //! constructor
      ParticleMaterialBase(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! @name material parameters
      //@{

      //! initial radius
      const double initRadius_;

      //! initial density
      const double initDensity_;

      //@}

      //! create material instance of matching type with parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override = 0;
    };

  }  // namespace PAR

}  // namespace MAT

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
