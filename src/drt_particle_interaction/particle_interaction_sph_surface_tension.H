/*---------------------------------------------------------------------------*/
/*! \file
\brief surface tension handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_SURFACE_TENSION_H
#define PARTICLE_INTERACTION_SPH_SURFACE_TENSION_H

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
  class SPHNeighborPairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHSurfaceTensionBase
  {
   public:
    //! constructor
    explicit SPHSurfaceTensionBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHSurfaceTensionBase() = default;

    //! init surface tension handler
    virtual void Init();

    //! setup surface tension handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! set current time
    virtual void SetCurrentTime(const double currenttime);

    //! insert surface tension evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) = 0;

    //! add surface tension contribution to acceleration field
    virtual void AddAccelerationContribution() = 0;

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

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! current time
    double time_;

    //! surface tension ramp function number
    const int surfacetensionrampfctnumber_;

    //! set of fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> fluidtypes_;

    //! set of boundary particle types
    std::set<PARTICLEENGINE::TypeEnum> boundarytypes_;
  };

  class SPHSurfaceTensionContinuumSurfaceForce : public SPHSurfaceTensionBase
  {
   public:
    //! constructor
    explicit SPHSurfaceTensionContinuumSurfaceForce(const Teuchos::ParameterList& params);

    //! init surface tension handler
    void Init() override;

    //! setup surface tension handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! insert surface tension evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! add surface tension contribution to acceleration field
    void AddAccelerationContribution() override;

   private:
    //! compute colorfield gradient
    void ComputeColorfieldGradient() const;

    //! compute unit wall normal
    void ComputeUnitWallNormal() const;

    //! compute wall distance
    void ComputeWallDistance() const;

    //! compute wall correction factor
    void ComputeWallCorrectionFactor();

    //! extrapolate colorfield gradient at triple point
    void ExtrapolateColorfieldGradientAtTriplePoint() const;

    //! compute interface normal
    void ComputeInterfaceNormal() const;

    //! correct normal vector of particles close to triple point
    void CorrectTriplePointNormal() const;

    //! compute surface tension contribution
    void ComputeSurfaceTensionContribution() const;

    //! compute temperature gradient driven contribution
    void ComputeTempGradDrivenContribution() const;

    //! colorfield gradient and wall distance of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh cfgwalldisttorefresh_;

    //! colorfield gradient and interface normal of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh cfgintnormtorefresh_;

    //! constant factor of surface tension coefficient
    const double alpha0_;

    //! static contact angle
    const double staticcontactangle_;

    //! factor of dependence of surface tension coefficient on temperature
    const double alphaT_;

    //! reference temperature for surface tension coefficient
    const double reftemp_;

    //! wall correction factor
    std::vector<std::vector<std::vector<double>>> f_i_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
