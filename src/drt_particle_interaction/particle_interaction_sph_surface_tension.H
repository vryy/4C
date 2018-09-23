/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_surface_tension.H

\brief surface tension handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_SURFACE_TENSION_H
#define PARTICLE_INTERACTION_SPH_SURFACE_TENSION_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_typedefs.H"

#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
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
  class SPHMaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHSurfaceTensionBase
  {
   public:
    //! constructor
    explicit SPHSurfaceTensionBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHSurfaceTensionBase(){};

    //! init surface tension handler
    virtual void Init();

    //! setup surface tension handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::SPHMaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of surface tension handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of surface tension handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current time
    virtual void SetCurrentTime(const double currenttime);

    //! insert surface tension evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) = 0;

    //! add surface tension contribution to acceleration field
    virtual void AddAccelerationContribution() const = 0;

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
    std::shared_ptr<PARTICLEINTERACTION::SPHMaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! current time
    double time_;

    //! surface tension ramp function number
    const int surfacetensionrampfctnumber_;

    //! have boundary particles
    bool haveboundaryparticles_;
  };

  class SPHSurfaceTensionContinuumSurfaceForce : public SPHSurfaceTensionBase
  {
   public:
    //! constructor
    explicit SPHSurfaceTensionContinuumSurfaceForce(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHSurfaceTensionContinuumSurfaceForce(){};

    //! init surface tension handler
    virtual void Init();

    //! insert surface tension evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! add surface tension contribution to acceleration field
    virtual void AddAccelerationContribution() const override;

   private:
    //! compute colorfield gradient
    virtual void ComputeColorfieldGradient() const;

    //! compute wall normal and distance
    virtual void ComputeWallNormalAndDistance() const;

    //! refresh colorfield gradient and wall distance
    virtual void RefreshColorfieldGradientAndWallDistance() const;

    //! extrapolate colorfield gradient at triple point
    virtual void ExtrapolateColorfieldGradientAtTriplePoint() const;

    //! compute interface normal
    virtual void ComputeInterfaceNormal() const;

    //! correct normal vector of particles close to triple point
    virtual void CorrectTriplePointNormal() const;

    //! refresh colorfield gradient and interface normal
    virtual void RefreshColorfieldGradientAndInterfaceNormal() const;

    //! compute curvature and acceleration contribution
    virtual void ComputeCurvatureAndAddAccelerationContribution() const;

    //! surface tension coefficient
    const double surfacetensioncoefficient_;

    //! static contact angle
    const double staticcontactangle_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
