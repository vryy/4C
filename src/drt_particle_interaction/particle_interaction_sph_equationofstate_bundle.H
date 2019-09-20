/*---------------------------------------------------------------------------*/
/*! \file
\brief class holding all equation of state handlers

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_H
#define PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

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

    //! destructor
    ~SPHEquationOfStateBundle(){};

    //! init equation of state bundle
    void Init(const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial);

    //! setup equation of state bundle
    void Setup();

    //! write restart of equation of state bundle
    void WriteRestart(const int step, const double time) const;

    //! read restart of equation of state bundle
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! return pointer to specific equation of state
    inline const PARTICLEINTERACTION::SPHEquationOfStateBase* GetPtrToSpecificEquationOfState(
        PARTICLEENGINE::TypeEnum typeEnum) const
    {
      return phasetypetoequationofstate_[typeEnum].get();
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
#endif
