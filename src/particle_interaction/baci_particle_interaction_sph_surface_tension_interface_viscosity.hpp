/*---------------------------------------------------------------------------*/
/*! \file
\brief interface viscosity handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_INTERFACE_VISCOSITY_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_INTERFACE_VISCOSITY_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_inpar_particle.hpp"
#include "baci_particle_engine_enums.hpp"
#include "baci_particle_engine_typedefs.hpp"

BACI_NAMESPACE_OPEN

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
  class SPHKernelBase;
  class MaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
  class SPHArtificialViscosity;
}  // namespace PARTICLEINTERACTION

namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialSPHFluid;
  }
}  // namespace MAT

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHInterfaceViscosity
  {
   public:
    //! constructor
    explicit SPHInterfaceViscosity(const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 12/2020
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~SPHInterfaceViscosity();

    //! init interface viscosity handler
    void Init();

    //! setup interface viscosity handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! compute interface viscosity contribution
    void ComputeInterfaceViscosityContribution() const;

   private:
    //! init artificial viscosity handler
    void InitArtificialViscosityHandler();

    //! compute interface viscosity contribution (particle contribution)
    void ComputeInterfaceViscosityParticleContribution() const;

    //! compute interface viscosity contribution (particle-boundary contribution)
    void ComputeInterfaceViscosityParticleBoundaryContribution() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! artificial viscosity handler
    std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity> artificialviscosity_;

    //! pointer to fluid material of particle types
    std::vector<const MAT::PAR::ParticleMaterialSPHFluid*> fluidmaterial_;

    //! liquid particle type
    PARTICLEENGINE::TypeEnum liquidtype_;

    //! gas particle type
    PARTICLEENGINE::TypeEnum gastype_;

    //! set of fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> fluidtypes_;

    //! set of boundary particle types
    std::set<PARTICLEENGINE::TypeEnum> boundarytypes_;

    //! artificial viscosity on liquid-gas interface
    const double artvisc_lg_int_;

    //! artificial viscosity on solid-liquid interface
    const double artvisc_sl_int_;

    //! transition reference temperature
    const double trans_ref_temp_;

    //! transition temperature difference for interface viscosity evaluation
    const double trans_dT_intvisc_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
