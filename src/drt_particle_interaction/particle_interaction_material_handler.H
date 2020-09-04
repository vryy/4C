/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material handler for particle simulations
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_MATERIAL_HANDLER_H
#define PARTICLE_INTERACTION_MATERIAL_HANDLER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_mat/particle_material_base.H"
#include "../drt_mat/particle_material_thermo.H"
#include "../drt_mat/particle_material_sph_fluid.H"
#include "../drt_mat/particle_material_sph_boundary.H"
#include "../drt_mat/particle_material_dem.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class MaterialHandler final
  {
   public:
    //! constructor
    explicit MaterialHandler(const Teuchos::ParameterList& params);

    //! destructor
    ~MaterialHandler() = default;

    //! init particle material handler
    void Init();

    //! setup particle material handler
    void Setup();

    //! return pointer to particle material parameter
    inline const MAT::PAR::ParticleMaterialBase* GetPtrToParticleMatParameter(
        PARTICLEENGINE::TypeEnum type_i) const
    {
      return phasetypetoparticlematpar_[type_i];
    }

    //! get particle types of stored particle material parameters
    inline const std::set<PARTICLEENGINE::TypeEnum> GetParticleTypes() const
    {
      return storedtypes_;
    };

   private:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! relate particle types to particle material parameters
    std::vector<const MAT::PAR::ParticleMaterialBase*> phasetypetoparticlematpar_;

    //! set of particle types of stored particle material parameters
    std::set<PARTICLEENGINE::TypeEnum> storedtypes_;
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
