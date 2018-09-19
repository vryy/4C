/*----------------------------------------------------------------------*/
/*!
\file particle_sph_rendering.H

\brief rendering routine for smoothed particle hydrodynamics

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/
#ifndef particle_sph_RENDERING_H
#define particle_sph_RENDERING_H

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_particle_old.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace PARTICLE
{
  class Algorithm;
  class ParticleSPHInteractionHandler;
  class WeightFunction_Base;
}  // namespace PARTICLE

namespace DRT
{
  class ResultTest;
}  // namespace DRT

/*----------------------------------------------------------------------*
 | rendering routine for smoothed particle hydrodynamics                |
 *----------------------------------------------------------------------*/
namespace PARTICLE
{
  class Rendering
  {
   public:
    //! constructor
    Rendering(Teuchos::RCP<PARTICLE::Algorithm> particleAlgorithm,
        Teuchos::RCP<PARTICLE::ParticleSPHInteractionHandler> interHandler,
        Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler);

    //! virtual destructor
    virtual ~Rendering() { ; }

    //! update the rendering vectors
    virtual void UpdateRenderingVectors(Teuchos::RCP<DRT::Discretization> pDiscret,
        Teuchos::RCP<const Epetra_Vector> pDis, Teuchos::RCP<const Epetra_Vector> pVel,
        Teuchos::RCP<const Epetra_Vector> pAcc, Teuchos::RCP<const Epetra_Vector> pVelMod,
        Teuchos::RCP<const Epetra_Vector> pDensity, Teuchos::RCP<const Epetra_Vector> pRadius,
        Teuchos::RCP<const Epetra_Vector> pPressure, Teuchos::RCP<const Epetra_Vector> pMass);

    //! output rendering state
    virtual void OutputState();

    //! clear rendering state
    virtual void ClearState();

    //! read-only rendering discretization
    virtual Teuchos::RCP<const DRT::Discretization> GetRenderingDiscret() { return discret_; }

    //! read-only rendering velocity
    virtual Teuchos::RCP<const Epetra_Vector> GetRenderingVelocity() { return vel_; }

    //! read-only rendering acceleration
    virtual Teuchos::RCP<const Epetra_Vector> GetRenderingAcceleration() { return acc_; }

    //! read-only rendering modified velocity
    virtual Teuchos::RCP<const Epetra_Vector> GetRenderingVelocityMod() { return velmod_; }

    //! read-only rendering density
    virtual Teuchos::RCP<const Epetra_Vector> GetRenderingDensity() { return density_; }

    //! read-only rendering pressure
    virtual Teuchos::RCP<const Epetra_Vector> GetRenderingPressure() { return pressure_; }

    //! create result test for rendering
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

   private:
    //! get neighboring rendering nodes
    std::list<DRT::Node*> GetNeighboringRenderingNodes(const int binId);

    //! rendering discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! pointer to underlying particle algorithm
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //! pointer to interaction handler
    Teuchos::RCP<PARTICLE::ParticleSPHInteractionHandler> interHandler_;

    //! pointer to weight function handler
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    //! type of rendering
    INPAR::PARTICLEOLD::RenderingType renderingType_;

    //! output rendering on discretization and/or in matlab format
    INPAR::PARTICLEOLD::RenderingOutput renderingOutput_;

    //! consider boundary particles in rendering routine
    INPAR::PARTICLEOLD::RenderingBdryPart renderingBdryParticle_;

    //! average rendering vectors over time steps
    int avrgRendering_;

    //! bins to rendering nodes map
    std::map<int, std::list<DRT::Node*>> binsToRenderingNodes_;

    //! rendering vectors
    Teuchos::RCP<Epetra_Vector> vel_;  //! rendering velocity

    Teuchos::RCP<Epetra_Vector> acc_;  //! rendering acceleration

    Teuchos::RCP<Epetra_Vector> velmod_;  //! rendering modified velocity

    Teuchos::RCP<Epetra_Vector> density_;  //! rendering density

    Teuchos::RCP<Epetra_Vector> pressure_;  //! rendering pressure
  };

}  // namespace PARTICLE

#endif  // particle_sph_RENDERING_H
