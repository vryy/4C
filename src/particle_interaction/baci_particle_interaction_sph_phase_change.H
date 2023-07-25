/*---------------------------------------------------------------------------*/
/*! \file
\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_SPH_PHASE_CHANGE_H
#define BACI_PARTICLE_INTERACTION_SPH_PHASE_CHANGE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_typedefs.H"

#include "baci_inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class SPHEquationOfStateBundle;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHPhaseChangeBase
  {
   public:
    //! constructor
    explicit SPHPhaseChangeBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHPhaseChangeBase() = default;

    //! init phase change handler
    virtual void Init();

    //! setup phase change handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle);

    //! evaluate phase change
    virtual void EvaluatePhaseChange(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const = 0;

   protected:
    //! evaluate phase change from below to above phase
    virtual void EvaluatePhaseChangeFromBelowToAbovePhase(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase,
        std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>>&
            particlestoinsert) const final;

    //! evaluate phase change from above to below phase
    virtual void EvaluatePhaseChangeFromAboveToBelowPhase(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase,
        std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>>&
            particlestoinsert) const final;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! phase below transition value
    PARTICLEENGINE::TypeEnum belowphase_;

    //! phase above transition value
    PARTICLEENGINE::TypeEnum abovephase_;

    //! transition state of phase change
    PARTICLEENGINE::StateEnum transitionstate_;

    //! transition value of phase change
    double transitionvalue_;

    //! hysteresis gap at transition value
    double hysteresisgap_;
  };

  class SPHPhaseChangeOneWayScalarBelowToAbove : public SPHPhaseChangeBase
  {
   public:
    //! constructor
    explicit SPHPhaseChangeOneWayScalarBelowToAbove(const Teuchos::ParameterList& params);

    //! evaluate phase change
    void EvaluatePhaseChange(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const override;
  };

  class SPHPhaseChangeOneWayScalarAboveToBelow : public SPHPhaseChangeBase
  {
   public:
    //! constructor
    explicit SPHPhaseChangeOneWayScalarAboveToBelow(const Teuchos::ParameterList& params);

    //! evaluate phase change
    void EvaluatePhaseChange(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const override;
  };

  class SPHPhaseChangeTwoWayScalar : public SPHPhaseChangeBase
  {
   public:
    //! constructor
    explicit SPHPhaseChangeTwoWayScalar(const Teuchos::ParameterList& params);

    //! evaluate phase change
    void EvaluatePhaseChange(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
