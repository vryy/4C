/*----------------------------------------------------------------------*/
/*!
\file adapter_particle.H

\brief Particle field adapter

\level 1

\maintainer Georg Hammerl
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

namespace LINALG
{
  class MapExtractor;
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

    /// start new time step
    virtual void PrepareTimeStep() = 0;

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

    //! Access to rcp pointer to noderowmap of discretization
    virtual const Epetra_Map* NodeRowMapView() = 0;

    //! node map of particles
    virtual Teuchos::RCP<const Epetra_Map> NodeRowMap() = 0;

    //! Access output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() = 0;

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

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

    //! unknown temperatures \f$T_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturen() = 0;

    //! unknown pressure \f$P_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressuren() = 0;

    //! unknown pressure \f$D_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densityn() = 0;

    //! Return velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    //! Return accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() = 0;

    //! Return temperatures \f$T_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturenp() = 0;

    //! Return pressure \f$P_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressurenp() = 0;

    //! Return density \f$D_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densitynp() = 0;

    //! unknown angular velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVeln() = 0;

    //! unknown angular accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccn() = 0;

    //! Return radii of particles \f$R_{0}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius0() = 0;

    //! Return radii of particles \f$R_{n}\f$ or \f$R_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius() = 0;

    //! Return latent heat of particles \f$L_{n}\f$ or \f$L_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> LatentHeat() = 0;

    //! Return first time derivative of radii \f$\dot{R}_{n}\f$ or \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> RadiusDot() = 0;

    //! Return mass of particles
    virtual Teuchos::RCP<const Epetra_Vector> Mass() = 0;

    //! Return inertia of particles
    virtual Teuchos::RCP<const Epetra_Vector> Inertia() = 0;

    //! Return initial density
    virtual double InitDensity() = 0;

    //! Return initial radius
    virtual double InitRadius() = 0;

    //! Return dismembering radius for melting particles
    virtual double DismemberRadius() = 0;

    //! Return specific heat -constant pressure - solid state of particles
    virtual double CPS() = 0;

    //! Return specific heat -constant pressure - liquid state of particles
    virtual double CPL() = 0;

    //! Return specific latent heat - solid <-> liquid
    virtual double MaxSpecificLatentHeat() = 0;

    //! transition temperature - solid <-> liquid
    virtual double SL_TransitionTemperature() = 0;

    //! does collision handler exist?
    virtual bool HaveCollHandler() = 0;

    //@}

    //! @name Extract methods
    //@{

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() = 0;

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() = 0;

    //! extract accelerations
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp() = 0;

    //! extract temperatures
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessTemperaturenp() = 0;

    //! extract pressures
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessPressurenp() = 0;

    //! extract densities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensitynp() = 0;

    //! extract SL latent heat
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSL_latentHeat() = 0;

    //! extract initial radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius0() = 0;

    //! extract radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius() = 0;

    //! extract latent heat
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessLatentHeat() = 0;

    //! wxtract first time derivative of radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusDot() = 0;

    //! extract mass vector
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessMass() = 0;

    //! extract inertia vector
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessInertia() = 0;

    //@}

    //! @name Time step helpers
    //@{

    /// Get the current time
    virtual double TimeOld() const = 0;

    /// Return target time \f$t_{n+1}\f$
    virtual double Time() const = 0;

    /// Return current step number $n$
    virtual int StepOld() const = 0;

    //! Return current step number $n+1$
    virtual int Step() const = 0;

    /// Set time and step
    virtual void SetTimeStep(const double time, const int step) = 0;

    /// Get time step size \f$\Delta t_n\f$
    virtual double Dt() const = 0;

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
