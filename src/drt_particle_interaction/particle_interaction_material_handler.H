/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_material_handler.H

\brief particle material handler for particle simulations

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
#ifndef PARTICLE_INTERACTION_MATERIAL_HANDLER_H
#define PARTICLE_INTERACTION_MATERIAL_HANDLER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_mat/particle_material_base.H"
#include "../drt_mat/particle_material_sph.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class BaseMaterialHandler
  {
   public:
    //! constructor
    explicit BaseMaterialHandler(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~BaseMaterialHandler(){};

    //! init particle material handler
    virtual void Init();

    //! setup particle material handler
    virtual void Setup();

    //! write restart of particle material handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of particle material handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

   protected:
    //! map particle types to particle material parameters
    virtual void MapParticleTypeToMatParameter(
        PARTICLEENGINE::TypeEnum particleType, int matID) = 0;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;
  };

  class SPHMaterialHandler : public BaseMaterialHandler
  {
   public:
    //! constructor
    explicit SPHMaterialHandler(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHMaterialHandler(){};

    //! get reference to particle material parameter map
    const std::map<PARTICLEENGINE::TypeEnum, const MAT::PAR::ParticleMaterialSPH*>&
    GetRefToParticleMatParMap() const
    {
      return phasetypetoparticlematpar_;
    };

    //! return pointer to particle material parameter
    const MAT::PAR::ParticleMaterialSPH* GetPtrToParticleMatParameter(
        PARTICLEENGINE::TypeEnum particleType) const;

   protected:
    //! map particle types to particle material parameters
    virtual void MapParticleTypeToMatParameter(
        PARTICLEENGINE::TypeEnum particleType, int matID) override;

    //! map relating particle types to particle material parameters
    std::map<PARTICLEENGINE::TypeEnum, const MAT::PAR::ParticleMaterialSPH*>
        phasetypetoparticlematpar_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
