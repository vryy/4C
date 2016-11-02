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

namespace MAT
{
  namespace PAR
  {
    class ParticleMat;
    class ExtParticleMat;
  }
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
    virtual void SetExternalDerivativeChangers(Teuchos::RCP<Epetra_MultiVector> iforce) = 0;

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

    //! Read-only displacements \f$D_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() = 0;

    //! Read-only velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    //! Read-only accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() = 0;

    //! Read-only angular-velocities \f$angV_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVeln() = 0;

    //! Read-only angular-accelerations \f$angA_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccn() = 0;

    //! Read-only radii \f$R_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radiusn() = 0;

    //! Read-only time derivative of radii \f$\dot{R}_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> RadiusDotn() = 0;

    //! Read-only densities \f$\rho_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densityn() = 0;

    //! Read-only time derivative of densities \f$\dot{\rho}_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> DensityDotn() = 0;

    //! Read-only specific enthalpies \f$H_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyn() = 0;

    //! Read-only time derivative of specific enthalpies \f$\dot{H}_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyDotn() = 0;

    //! Read-only specific temperatures \f$T_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturen() = 0;



    //! Read-only displacements \f$D_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    //! Read-only velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    //! Read-only accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() = 0;

    //! Read-only angular-velocities \f$angV_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVelnp() = 0;

    //! Read-only angular-accelerations \f$angA_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccnp() = 0;

    //! Read-only radii \f$R_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radiusnp() = 0;

    //! Read-only time derivative of radii \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> RadiusDotnp() = 0;

    //! Read-only densities \f$\rho_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densitynp() = 0;

    //! Read-only time derivative of densities \f$\dot{\rho}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> DensityDotnp() = 0;

    //! Read-only specific enthalpies \f$H_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpynp() = 0;

    //! Read-only time derivative of specific enthalpies \f$\dot{H}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyDotnp() = 0;

    //! Read-only temperatures \f$T_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturenp() = 0;

    //! Read-only radii \f$R_0\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius0() = 0;

    //! Read-only masses \f$M_{n}\f$ or \f$M_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Mass() = 0;

    //! Read-only inertiae \f$I_{n}\f$ or \f$I_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Inertia() = 0;

    //! Read-only pressures \f$P_{n}\f$ or \f$P_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressure() = 0;



    //! does collision handler exist?
    virtual bool HaveCollHandler() = 0;

    //! does interaction handler exist?
    virtual bool HaveInterHandler() = 0;

    //! check state vector
    virtual void CheckStateVector(std::string vecName, const Teuchos::RCP<const Epetra_Vector> vec, bool trg_showVec) = 0;

    //@}

    //! @name Extract methods
    //@{

    //! Extract displacements \f$D_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() = 0;

    //! Extract velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() = 0;

    //! Extract accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp() = 0;

    //! Extract angular-velocities \f$angV_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAngVelnp() = 0;

    //! Extract angular-accelerations \f$angA_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAngAccnp() = 0;

    //! Extract radii \f$R_{n}\f$ for legacy code
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusn() = 0;

    //! Extract radii \f$R_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusnp() = 0;

    //! Extract time derivative radii \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusDotnp() = 0;

    //! Extract densities \f$\rho_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensitynp() = 0;

    //! Extract time derivative densities \f$\dot{\rho}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensityDotnp() = 0;

    //! Extract specific enthalpies \f$H_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSpecEnthalpynp() = 0;

    //! Extract time derivative specific enthalpies \f$\dot{H}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSpecEnthalpyDotnp() = 0;

    //! extract radii \f$R_{0}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius0() = 0;

    //! Extract masses \f$M_{n}\f$ or \f$M_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessMass() = 0;

    //! Extract inertiae \f$I_{n}\f$ or \f$I_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessInertia() = 0;

    //! Extract pressures \f$P_{n}\f$ or \f$P_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessPressure() = 0;

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
