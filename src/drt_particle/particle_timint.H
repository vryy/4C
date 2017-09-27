/*----------------------------------------------------------------------*/
/*!
\file particle_timint.H

\brief Time integration for particle dynamics

\level 1

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
/* definitions */
#ifndef PARTICLE_TIMINT_H
#define PARTICLE_TIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_utils.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_inpar/inpar_particle.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>
#include <Epetra_Time.h>
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_io/io.H"

#include <string>
#include <iostream>
#include <fstream>

// forward declarations
namespace DRT
{
  class Discretization;
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
  class ParticleCollisionHandlerBase;
  class ParticleMeshFreeInteractionHandler;
}

namespace MAT
{
  namespace PAR
  {
    class ParticleMat;
    class ExtParticleMat;
  }
}

/*----------------------------------------------------------------------*/
//! Particle Dynamics
namespace PARTICLE
{

   /// Front-end for particle dynamics for integrating in time.
  class TimInt : public ADAPTER::Particle
  {

  public:

    //! @name Life
    //@{

    //! Print logo
    void Logo();

    //! Constructor
    TimInt
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimInt() { ; }

    //@}

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations
    virtual void DetermineMassDampConsistAccel();

    //! Apply Dirichlet boundary conditions on provided state vectors
    virtual void ApplyDirichletBC
    (
      const double time,  //!< at time
      Teuchos::RCP<Epetra_Vector> dis,  //!< displacements
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> vel,  //!< velocities
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> acc,  //!< accelerations
                                       //!< (may be Teuchos::null)
      bool recreatemap  //!< recreate mapextractor/toggle-vector
                        //!< which stores the DOF IDs subjected
                        //!< to Dirichlet BCs
                        //!< This needs to be true if the bounded DOFs
                        //!< have been changed.
    );

    //! acceleration is applied from given forces
    virtual void ComputeAcc(
      Teuchos::RCP<Epetra_Vector> f_contact,
      Teuchos::RCP<Epetra_Vector> m_contact,
      Teuchos::RCP<Epetra_Vector> acc,
      Teuchos::RCP<Epetra_Vector> ang_acc
      );

    //! determine acceleration in meshfree case
    virtual void DetermineMeshfreeDensAndAcc(Teuchos::RCP<Epetra_Vector> acc,
                                             Teuchos::RCP<Epetra_Vector> accmod,
                                             Teuchos::RCP<Epetra_Vector> velConv,
                                             Teuchos::RCP<Epetra_Vector> acc_A,
                                             const double time,
                                             const double dt);

    //! Prepare time step
    virtual void PrepareTimeStep();

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;

    /// initialize time integration
    virtual void Init();

    /// update at time step end
    virtual void Update() = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible.
    virtual void UpdateStepState();

    //! Update time and step counter
    virtual void UpdateStepTime();

    //! State vectors are updated according to the new distribution of particles
    virtual void UpdateStatesAfterParticleTransfer();

    //! particle algorithm is set from outside
    virtual void SetParticleAlgorithm(Teuchos::RCP<PARTICLE::Algorithm> algorithm) { particle_algorithm_ = algorithm; }

    //! set vector of particle forces on struct dis from outside
    virtual void SetFstructure(Teuchos::RCP<Epetra_FEVector> f_structure) { f_structure_ = f_structure ; }

    //! Calculate all output quantities depending on the constitutive model
    //  (and, hence, on a potential material history)
    virtual void PrepareOutput();

    //! Calculate kinetic, internal and external energy
    virtual void DetermineEnergy();

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    //@}


    //! @name Output
    //@{

    //! Output to file
    //! This routine always prints the last converged state, i.e.
    //! \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    virtual void OutputStep(bool forced_writerestart = false);

    //! Write restart
    virtual void OutputRestart
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );
    //! Output displacements, velocities and accelerations
    //! and more system vectors
    virtual void OutputState
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Energy output
    virtual void OutputEnergy();

    //! Meshfree rendering output
    virtual void PerformMeshfreeRendering(bool clearstate, bool writeoutput);

    //! output particle statistics
    void OutputParticleStatistics();

    //! Check whether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (not energyfile_.is_null()) return true;
      else return false;
    }

    //! Attach file handle for energy file #energyfile_
    virtual void AttachEnergyFile();

    //! attach file handle for particle statistics
    void AttachParticleStatisticsFile();
    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    virtual Teuchos::RCP<DRT::Discretization> Discretization(){ return discret_; }

    //! Access to rcp pointer to dofrowmap of discretization
    virtual const Epetra_Map* DofRowMapView();

    //! dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    //! Access to rcp pointer to noderowmap of discretization
    virtual const Epetra_Map* NodeRowMapView();

    //! node map of particles
    virtual Teuchos::RCP<const Epetra_Map> NodeRowMap();

    //! Access output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
    {
      return output_;
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return dbcmaps_;
    }

    //! Read restart values
    virtual void ReadRestart
    (
      const int step //!< restart step
    );

    //! Read and set restart state
    virtual void ReadRestartState();

    //! Read-only displacements \f$D_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn(){ return (*dis_)(0); }

    //! Read-only velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() { return (*vel_)(0); }

    //! Read-only accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { return (*acc_)(0); }

    //! Read-only modified velocities \f$Vm_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velmodn() { return (*velmod_)(0); }

    //! Read-only modified accelerations \f$Am_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accmodn() { return (*accmod_)(0); }

    //! Read-only angular-velocities \f$angV_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVeln() { return (*angVel_)(0); }

    //! Read-only angular-accelerations \f$angA_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccn() { return (*angAcc_)(0); }

    //! Read-only radii \f$R_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radiusn() { return (*radius_)(0); }

    //! Read-only densities \f$\rho_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densityn() { return (*density_)(0); }

    //! Read-only time derivative of densities \f$\dot{\rho}_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> DensityDotn() { return (*densityDot_)(0); }

    //! Read-only specific enthalpies \f$H_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyn() { return (*specEnthalpy_)(0); }

    //! Read-only time derivative of specific enthalpies \f$\dot{H}_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyDotn() { return (*specEnthalpyDot_)(0); }

    //! Read-only displacements \f$D_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp(){ return disn_; }

    //! Read-only velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() { return veln_; }

    //! Read-only accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { return accn_; }

    //! Read-only modified velocities \f$Vm_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velmodnp() { return velmodn_; }

    //! Read-only modified accelerations \f$Am_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accmodnp() { return accmodn_; }

    //! Read-only angular-velocities \f$angV_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVelnp() { return angVeln_; }

    //! Read-only angular-accelerations \f$angA_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccnp() { return angAccn_; }

    //! Read-only radii \f$\rho_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radiusnp() { return radiusn_; }

    //! Read-only densities \f$\rho_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densitynp() { return densityn_; }

    //! Read-only time derivative of densities \f$\dot{\rho}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> DensityDotnp() { return densityDotn_; }

    //! Read-only specific enthalpies \f$H_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpynp() { return specEnthalpyn_; }

    //! Read-only time derivative of specific enthalpies \f$\dot{H}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> SpecEnthalpyDotnp() { return specEnthalpyDotn_; }

    //! Read-only radii \f$R_0\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius0() {return radius0_; }

    //! Read-only time derivative of radii \f$\dot{R}_{n}\f$ or \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> RadiusDot() { return radiusDot_; }

    //! Read-only masses \f$M_{n}\f$ or \f$M_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Mass() { return mass_; }

    //! Read-only inertiae \f$I_{n}\f$ or \f$I_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Inertia() { return inertia_; }

    //! Read-only temperatures \f$T_{n}\f$ or \f$T_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperature() { return temperature_; }

    //! Read-only pressures \f$P_{n}\f$ or \f$P_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressure() { return pressure_; }

    //! Read-only force on wall discretization \f$f_{n}\f$ or \f$f_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_FEVector> Fstructure(){ return f_structure_; }

    //! does collision handler exist?
    virtual bool HaveCollHandler() {return (collhandler_!=Teuchos::null); }

    //! does interaction handler exist?
    virtual bool HaveInterHandler() {return (interHandler_!=Teuchos::null); }

    //! check state vectors (debugging)
    virtual void CheckStateVector(std::string vecName, const Teuchos::RCP<const Epetra_Vector> vec, bool trg_showVec = false);

    //! return kinetic energy
    double KineticEnergy() const {return kinergy_;};

    //! return internal energy
    double InternalEnergy() const {return intergy_;};

    //! return external energy
    double ExternalEnergy() const {return extergy_;};

    //! return maximum particle-particle or particle-wall penetration
    double MaximumPenetration() const;
    //@}

    //! @name Extract methods
    //@{

    //! Extract displacements \f$D_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp(){ return disn_; }

    //! Extract velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp(){ return veln_; }

    //! Extract accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp() { return accn_; }

    //! Extract modified velocities \f$Vm_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelmodnp() { return velmodn_; }

    //! Extract modified accelerations \f$Am_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccmodnp() { return accmodn_; }

    //! Extract angular-velocities \f$angV_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAngVelnp(){ return angVeln_; }

    //! Extract angular-accelerations \f$angA_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAngAccnp() { return angAccn_; }

    //! Extract radii \f$R_{n}\f$ for legacy code
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusn() { return (*radius_)(0); }

    //! Extract radii \f$R_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusnp() { return radiusn_; }

    //! Extract densities \f$\rho_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensitynp() { return densityn_; }

    //! Extract time derivative densities \f$\dot{\rho}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensityDotnp() { return densityDotn_; }

    //! Extract specific enthalpies \f$H_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSpecEnthalpynp(){ return specEnthalpyn_; }

    //! Extract time derivative specific enthalpies \f$\dot{H}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSpecEnthalpyDotnp(){ return specEnthalpyDotn_; }

    //! extract radii \f$R_{0}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius0() { return radius0_; }

    //! Extract time derivative radii \f$\dot{R}_{n}\f$ or \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusDot() { return radiusDot_; }

    //! Extract masses \f$M_{n}\f$ or \f$M_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessMass() { return mass_; }

    //! Extract masses \f$I_{n}\f$ or \f$I_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessInertia() { return inertia_; }

    //! Extract pressures \f$P_{n}\f$ or \f$P_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessPressure() { return pressure_; }

    //! Extract force on wall discretization \f$f_{n}\f$ or \f$f_{n+1}\f$
    virtual Teuchos::RCP<Epetra_FEVector> WriteAccessFstructure(){ return f_structure_; }

    //! Forces \f$F_{n}\f$ or \f$F_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessFifc() { return fifc_; }

    //@}

    //! @name Time step helpers
    //@{

    //! Return current time map
    virtual double TimeOld() const { return (*time_)[0]; }

    //! Return target time \f$t_{n+1}\f$
    virtual double Time() const { return timen_; }

    //! Return current step number $n$
    virtual int StepOld() const { return step_; }

    //! Return current step number $n+1$
    virtual int Step() const { return stepn_; }

    //! Set current time and step number map
    virtual void SetTimeStep(const double time, const int step)
    {
      (*time_)[0] = time;
      step_ = step;
      timen_ = (*time_)[0] + (*dt_)[0];  // set target time to initial time plus step size
      stepn_ = step_ + 1;
      return;
    }

    //! Get time step size \f$\Delta t_n\f$
    virtual double Dt() const { return (*dt_)[0]; }

    //! overloadable function to apply external actions (i.e. forces, accelerations, etc) and to erase outdated vectors
    virtual void UpdateExtActions(bool init = false);

    virtual void GravityForces(Teuchos::RCP<Epetra_Vector> stateVector, const double extMulti = 1.0);

    virtual void GravityAcc(Teuchos::RCP<Epetra_Vector> stateVector, const double extMulti, const double time);

    //@}

  protected:

    //! @name General purpose algorithm members
    //@{

    //! Set initial fields
    virtual void SetInitialFields();

    //! Determine boundary particles
    virtual void DetermineBdryParticles();

    /* Update exporter objects if layout has changed*/
    virtual void UpdateExportersIfNecessary(const Epetra_BlockMap& oldnodemap, const Epetra_BlockMap& olddofmap);

    /* Update TimIntMStep state vector with the new map from discret_*/
    virtual void UpdateStateVectorMap(Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > &stateVector, bool trg_nodeVectorType = false);

    /* Update state vector with the new map from discret_*/
    virtual void UpdateStateVectorMap(Teuchos::RCP<Epetra_Vector> &stateVector, bool trg_nodeVectorType = false);

    /* Wrapper update state multi vector with the new map from discret_*/
    virtual void UpdateStateVectorMap(Teuchos::RCP<Epetra_MultiVector> &stateVector, bool trg_nodeVectorType = false);

    //! Read-only temperatures \f$T_{n}\f$
    virtual void UpdateTemperaturen();

    //! Read-only temperatures \f$T_{n+1}\f$
    virtual void UpdateTemperaturenp();

    //! Read-only temperatures \f$P_{n+1}\f$
    virtual void UpdatePressure();

    //! initialization of vector for visualization of the particle orientation
    virtual void InitializeOrientVector();

    //! update of vector for visualization of the particle orientation
    virtual void RotateOrientVector(double dt);

    //! wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
    virtual void WriteVector(const std::string name, Teuchos::RCP<Epetra_Vector> vec, const bool isdof = true);

    //! wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
    virtual void WriteVector(const std::string name, Teuchos::RCP<Epetra_MultiVector> vec, const bool isdof = true);

    //! wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
    virtual void WriteVector(const std::string name, Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > vec, const bool isdof = true);

    //! smart wrapper, the update happens only when vectors exist
    inline virtual void UpdateStateVector(const Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> >& stateVector,
        const Teuchos::RCP<Epetra_Vector>& stateVectorn)
    { if (stateVector != Teuchos::null && stateVectorn != Teuchos::null) { stateVector->UpdateSteps(*stateVectorn); } }


    //@}

    //! @name Protected Compute Methods
    //@{

    //@}

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int myrank_;  //!< ID of actual processor in parallel
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                  //!< containing non-overlapping
                                                  //!< map of global DOFs on Dirichlet
                                                  //!< boundary conditions
    Teuchos::RCP<Epetra_Vector> dbcdofs_;         //!< vector of dirichlet bc dofs
    //@}

    //! @name Printing and output
    //@{
    const Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    const bool printlogo_;  //!< true: enjoy your cuppa
    const int printscreen_;  //!< print infos to standard out every n steps
    FILE* errfile_;  //!< error file handle
    const bool printerrfile_;  //!< print infos to error file #errfile_
    const int writerestartevery_;  //!< write restart every given step;
                             //!< if 0, restart is not written
    const bool writestate_;  //!< write state on/off
    const bool writevelacc_;  //!< write velocity and acceleration on/off
    const int writeresultsevery_;  //!< write state/stress/strain every given step
    const int writeenergyevery_;  //!< write system energy every given step
    const int writerenderingevery_; //!< write meshfree rendering every given step
    const int avrgrenderingsteps_; //!< write meshfree rendering every given step
    const int writeparticlestatsevery_;  //!< write particle statistics every given step
    Teuchos::RCP<std::ofstream> energyfile_;  //!< output file for energy
    Teuchos::RCP<std::ofstream> particlestatsfile_;  //!< output file for particle statistics
    bool writeorientation_;  //!< write orientation of particle on/off

    double kinergy_; //!< kinetic energy
    double intergy_; //!< internal energy
    double bpintergy_; //!< internal energy of boundary particles (different definition of density/pressure)
    double extergy_; //!< external energy
    LINALG::Matrix<3,1> linmomentum_; //!< linear momentum
    //@}

    //! @name General control parameters
    //@{
    Teuchos::RCP<TIMINT::TimIntMStep<double> > time_;  //!< time \f$t_{n}\f$ of last converged step
    double timen_;  //!< target time \f$t_{n+1}\f$
    Teuchos::RCP<TIMINT::TimIntMStep<double> > dt_;  //!< time step size \f$\Delta t\f$
    double timemax_;  //!< final time \f$t_\text{fin}\f$
    int stepmax_;  //!< final step \f$N\f$
    int step_;  //!< time step index \f$n\f$
    int stepn_;  //!< time step index \f$n+1\f$
    int restart_;  //!< restart step
    //@}

    //! @name Global state vectors
    //@{

    //! Displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > dis_;
    //! Velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > vel_;
    //! Accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > acc_;
    //! Modified velocities \f${Vm}_{n}, Vm_{n-1}, ...\f$ resembling particle convection velocities
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > velmod_;
    //! Modified accelerations \f${Am}_{n}, Am_{n-1}, ...\f$ allowing to distinguish physical velocities and particle convection velocities
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > accmod_;
    //! Angular-velocities \f${angV}_{n}, angV_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > angVel_;
    //! Angular-accelerations \f${angA}_{n}, angA_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > angAcc_;
    //! Radii \f$R_{n}, R_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > radius_;
    //! Densities \f$\rho_{n}, \rho_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > density_;
    //! Time derivative of densities \f$\dot{\rho}_{n}, \dot{\rho}_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > densityDot_;
    //! Specific enthalpies \f$H_{n}, H_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > specEnthalpy_;
    //! Time derivative of specific enthalpies \f$\dot{H}_{n}, \dot{H}_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > specEnthalpyDot_;

    //! Displacements \f${D}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disn_;
    //! Velocities \f${V}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;
    //! Accelerations \f${A}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;
    //! Modified velocities \f${Vm}_{n+1}\f$ resembling particle convection velocities
    Teuchos::RCP<Epetra_Vector> velmodn_;
    //! Modified accelerations \f${Am}_{n+1}\f$ allowing to distinguish physical velocities and particle convection velocities
    Teuchos::RCP<Epetra_Vector> accmodn_;
    //! Angular-velocities \f${angV}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> angVeln_;
    //! Angular-acceleration \f${angA}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> angAccn_;
    //! Radii \f$R_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> radiusn_;
    //! Densities \f${\rho}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> densityn_;
    //! Time derivative of densities \f$\dot{\rho}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> densityDotn_;
    //! Specific enthalpies \f${H}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> specEnthalpyn_;
    //! Time derivative of specific enthalpies \f$\dot{H}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> specEnthalpyDotn_;

    //! additional external forces
    Teuchos::RCP<Epetra_Vector> fifc_;
    //! Vector for visualization of particle orientation
    Teuchos::RCP<Epetra_Vector> orient_;

    //! Radii \f$R_{0}\f$
    Teuchos::RCP<Epetra_Vector> radius0_;
    //! Time derivative of radii \f$\dot{R}_{n}\f$ or \f$\dot{R}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> radiusDot_;
    //! Masses \f$M_{n}\f$ or \f$M_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> mass_;
    //! Inertiae \f$I_{n}\f$ or \f$I_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> inertia_;
    //! Temperature \f$T_{n}\f$ or \f$T_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> temperature_;
    //! Pressures \f$P_{n}\f$ or \f$P_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> pressure_;
    //! Particle force on wall discretization \f$f_{n}\f$ or \f$f_{n+1}\f$
    Teuchos::RCP<Epetra_FEVector> f_structure_;

    //! color field
    Teuchos::RCP<Epetra_Vector> colorField_;
    //! free-surface particle type
    Teuchos::RCP<Epetra_Vector> fspType_;

    //! initial density for meshFree-pressure dynamics
    double initDensity_;
    //! rest density for meshFree-pressure dynamics defining the equation of state: pressure=c_0^2(density-refdensfac_*restDensity_)  with the speed of sound c_0^2=BULK_MODULUS/restDensity_
    double restDensity_;
    /// factor defining the equation of state: pressure=c_0^2(density-refdensfac_*restDensity_)  with the speed of sound c_0^2=BULK_MODULUS/restDensity_;
    double refdensfac_;

    /// global number of boundary particles over all processors
    int global_num_boundaryparticles_;


    //@}

    //! @name Global state vector transfer helpers
    //@{

    //! exporter for transferring dof based map from old to new layout
    Teuchos::RCP<Epetra_Export> dofmapexporter_;
    //! exporter for transferring node based map from old to new layout
    Teuchos::RCP<Epetra_Export> nodemapexporter_;

    //@}

    //! @name control flags for special features
    //@{

    //! flag for random distribution of particle radii
    const INPAR::PARTICLE::RadiusDistribution radiusdistribution_;

    //! bool to control whether radius is constant or variable in time
    const bool variableradius_;

    //! number of curve governing radius change
    const int radiuschangefunct_;

    //@}

    //! @name helper classes
    //@{

    //! particle algorithm that is connected to the time integration
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //! collision handler for particles
    Teuchos::RCP<PARTICLE::ParticleCollisionHandlerBase> collhandler_;

    //! interaction handler for particles
    Teuchos::RCP<PARTICLE::ParticleMeshFreeInteractionHandler> interHandler_;

    //@}

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_H
