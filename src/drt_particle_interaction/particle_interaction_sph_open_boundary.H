/*---------------------------------------------------------------------------*/
/*! \file
\brief open boundary handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_OPEN_BOUNDARY_H
#define PARTICLE_INTERACTION_SPH_OPEN_BOUNDARY_H

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
  class SPHKernelBase;
  class MaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHOpenBoundaryBase
  {
   public:
    //! constructor
    explicit SPHOpenBoundaryBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHOpenBoundaryBase(){};

    //! init open boundary handler
    virtual void Init();

    //! setup open boundary handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of open boundary handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of open boundary handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! prescribe open boundary states
    virtual void PrescribeOpenBoundaryStates(const double& evaltime) = 0;

    //! interpolate open boundary states
    virtual void InterpolateOpenBoundaryStates() = 0;

    //! check open boundary phase change
    virtual void CheckOpenBoundaryPhaseChange(const double maxinteractiondistance) final;

   protected:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! states of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh statestorefresh_;

    //! function id of prescribed state
    int prescribedstatefunctid_;

    //! outward normal
    std::vector<double> outwardnormal_;

    //! plane point
    std::vector<double> planepoint_;

    //! fluid phase
    PARTICLEENGINE::TypeEnum fluidphase_;

    //! open boundary phase
    PARTICLEENGINE::TypeEnum openboundaryphase_;
  };

  class SPHOpenBoundaryDirichlet : public SPHOpenBoundaryBase
  {
   public:
    //! constructor
    explicit SPHOpenBoundaryDirichlet(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHOpenBoundaryDirichlet() override{};

    //! init open boundary handler
    void Init() override;

    //! setup open boundary handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! prescribe open boundary states
    void PrescribeOpenBoundaryStates(const double& evaltime) override;

    //! interpolate open boundary states
    void InterpolateOpenBoundaryStates() override;
  };

  class SPHOpenBoundaryNeumann : public SPHOpenBoundaryBase
  {
   public:
    //! constructor
    explicit SPHOpenBoundaryNeumann(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHOpenBoundaryNeumann() override{};

    //! init open boundary handler
    void Init() override;

    //! setup open boundary handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! prescribe open boundary states
    void PrescribeOpenBoundaryStates(const double& evaltime) override;

    //! interpolate open boundary states
    void InterpolateOpenBoundaryStates() override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
