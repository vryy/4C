/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_equationofstate_bundle.H

\brief class holding all equation of state handlers for smoothed particle hydrodynamics (SPH)
interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_H
#define PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_BUNDLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2018 |
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
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHEquationOfStateBundle
  {
   public:
    //! constructor
    explicit SPHEquationOfStateBundle(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHEquationOfStateBundle(){};

    //! init equation of state bundle
    virtual void Init(std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial);

    //! setup equation of state bundle
    virtual void Setup();

    //! write restart of equation of state bundle
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of equation of state bundle
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get specific equation of state
    virtual std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBase> GetSpecificEquationOfState(
        PARTICLEENGINE::TypeEnum particleType) const;

   protected:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! equation of state handler for all particle types
    std::map<PARTICLEENGINE::TypeEnum, std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBase>>
        phasetypetoequationofstate_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
