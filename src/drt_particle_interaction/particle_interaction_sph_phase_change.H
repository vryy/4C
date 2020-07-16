/*---------------------------------------------------------------------------*/
/*! \file
\brief phase change handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_PHASE_CHANGE_H
#define PARTICLE_INTERACTION_SPH_PHASE_CHANGE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

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

    //! write restart of phase change handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of phase change handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! evaluate phase change
    virtual void EvaluatePhaseChange() const = 0;

   protected:
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
  };

  class SPHPhaseChangeTwoWayScalar : public SPHPhaseChangeBase
  {
   public:
    //! constructor
    explicit SPHPhaseChangeTwoWayScalar(const Teuchos::ParameterList& params);

    //! init phase change handler
    void Init() override;

    //! setup phase change handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle)
        override;

    //! evaluate phase change
    void EvaluatePhaseChange() const override;

   private:
    //! phase below transition value
    PARTICLEENGINE::TypeEnum belowphase_;

    //! phase above transition value
    PARTICLEENGINE::TypeEnum abovephase_;

    //! transition state of phase change
    PARTICLEENGINE::StateEnum transitionstate_;

    //! transition value of phase change
    double transitionvalue_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
