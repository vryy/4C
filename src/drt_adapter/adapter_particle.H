/*----------------------------------------------------------------------*/
/*!
\file adapter_particle.H

\brief Particle field adapter

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/13 |
 *----------------------------------------------------------------------*/
#ifndef ADAPTER_PARTICLE_H
#define ADAPTER_PARTICLE_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/13 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>


/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 12/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class ResultTest;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace PARTICLE
{
  class Algorithm;
}


/*----------------------------------------------------------------------*
 | adapting adapter                                         ghamm 12/13 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{

  class Particle
  {
  public:

    /// virtual to get polymorph destruction
    virtual ~Particle() { ; }

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent initial accelerations
    virtual void DetermineMassDampConsistAccel() = 0;

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;

    //! initialize particle problem
    virtual void Init() = 0;

    //! update at time step end
    virtual void Update() = 0;

    //! output results
    virtual void OutputStep(bool forced_writerestart = false) = 0;

    /// apply interface loads
    virtual void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce) = 0;

    //! State vectors are updated according to the new distribution of particles
    virtual void UpdateStatesAfterParticleTransfer() = 0;

    //! particle algorithm is set from outside
    virtual void SetParticleAlgorithm(Teuchos::RCP<PARTICLE::Algorithm> algorithm) = 0;

    //! Calculate all output quantities depending on the constitutive model
    virtual void PrepareOutput() = 0;

    /// create result test for particle problem
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;

    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    //! Access to rcp pointer to dofrowmap of discretization
    virtual const Epetra_Map* DofRowMapView() = 0;

    //! dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    //! Access output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() = 0;

    //! Read restart values
    virtual void ReadRestart(const int step) = 0;

    //! Read and set restart state
    virtual void ReadRestartState() = 0;

    //! known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() = 0;

    //! unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    //! unknown velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    //! unknown accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() = 0;

    //! Return velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    //! Return radii of particles
    virtual Teuchos::RCP<const Epetra_Vector> Radius() = 0;

    //! Return mass of particles
    virtual Teuchos::RCP<const Epetra_Vector> Mass() = 0;

    //! Return density
    virtual const double ParticleDensity() = 0;

    //@}

    //! @name Extract methods
    //@{

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() = 0;

    //! extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispn() = 0;

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVeln() = 0;

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() = 0;

    //! Return radii of particles
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius() = 0;

    //! extract mass vector
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessMass() = 0;

    //@}

    //! @name Time step helpers
    //@{

    /// Get the current time
    virtual double TimeOld() const = 0;

    /// Set time and step
    virtual void SetTimeStep(const double time, const int step) = 0;

    //@}

  };  // class Particle


  /// particle field solver
  class ParticleBaseAlgorithm
  {
  public:

    /// constructor
    explicit ParticleBaseAlgorithm(
      const Teuchos::ParameterList& prbdyn,
      Teuchos::RCP<DRT::Discretization> actdis
      );

    /// virtual destructor to support polymorph destruction
    virtual ~ParticleBaseAlgorithm() { ; }

    /// return particle field
    Teuchos::RCP<Particle> ParticleField() { return particle_; }

  private:

    /// setup particle time integrator
    void SetupTimInt(
      const Teuchos::ParameterList& prbdyn,
      Teuchos::RCP<DRT::Discretization> actdis
      );

    /// particle field
    Teuchos::RCP<Particle> particle_;

  };  // class ParticleBaseAlgorithm

}  // namespace ADAPTER


/*----------------------------------------------------------------------*/
#endif // ADAPTER_PARTICLE_H
