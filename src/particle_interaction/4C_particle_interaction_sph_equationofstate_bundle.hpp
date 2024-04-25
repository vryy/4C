/*---------------------------------------------------------------------------*/
/*! \file
\brief class holding all equation of state handlers
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHEquationOfStateBase;
  class MaterialHandler;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHEquationOfStateBundle final
  {
   public:
    //! constructor
    explicit SPHEquationOfStateBundle(const Teuchos::ParameterList& params);

    //! init equation of state bundle
    void Init(const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial);

    //! setup equation of state bundle
    void Setup();

    //! return pointer to specific equation of state
    inline const PARTICLEINTERACTION::SPHEquationOfStateBase* GetPtrToSpecificEquationOfState(
        PARTICLEENGINE::TypeEnum type_i) const
    {
      return phasetypetoequationofstate_[type_i].get();
    };

   private:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! equation of state handler for all particle types
    std::vector<std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateBase>>
        phasetypetoequationofstate_;

    //! set of particle types of stored equation of state handlers
    std::set<PARTICLEENGINE::TypeEnum> storedtypes_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
