/*----------------------------------------------------------------------*/
/*!
\file particle_timint.H

\brief Time integration for particle dynamics

\level 2

\maintainer Georg Hammerl
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef PARTICLE_TIMINT_H
#define PARTICLE_TIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_adapter/adapter_particle.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_inpar/inpar_cavitation.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Time.h>

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
    virtual void UpdateStepState() = 0;

    //! Update time and step counter
    virtual void UpdateStepTime();

    //! State vectors are updated according to the new distribution of particles
    virtual void UpdateStatesAfterParticleTransfer();

    //! Set initial fields in structure (e.g. initial velocities)
    virtual void SetInitialFields();

    //! particle algorithm is set from outside
    virtual void SetParticleAlgorithm(Teuchos::RCP<PARTICLE::Algorithm> algorithm) { particle_algorithm_ = algorithm; }

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

    //! Check whether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (not energyfile_.is_null()) return true;
      else return false;
    }

    //! Attach file handle for energy file #energyfile_
    virtual void AttachEnergyFile();

    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    {
      return discret_;
    }

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

    //! known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn(){ return (*dis_)(0); }

    //! unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp(){ return disn_; }

    //! unknown velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() { return (*vel_)(0); }

    //! unknown accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { return (*acc_)(0); }

    //! unknown temperatures \f$T_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturen() { return (*temperature_)(0); }

    //! unknown pressure \f$P_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressuren() { return (*pressure_)(0); }

    //! unknown pressure \f$D_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densityn() { return (*density_)(0); }

    //! Return velocities \f$V_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() { return veln_; }

    //! Return accelerations \f$A_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { return accn_; }

    //! Return temperatures \f$T_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Temperaturenp() { return temperaturen_; }

    //! Return pressure \f$P_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Pressurenp() { return pressuren_; }

    //! Return density \f$D_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Densitynp() { return densityn_; }

    //! unknown angular velocities \f$V_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngVeln() { return (*ang_vel_)(0); }

    //! unknown angular accelerations \f$A_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> AngAccn() { return (*ang_acc_)(0); }

    //! Return radii of particles \f$R_{0}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius0() {return radius0_; }

    //! Return radii of particles \f$R_{n}\f$ or \f$R_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Radius() { return radius_; }

    //! Return latent heat of particles \f$L_{n}\f$ or \f$L_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> LatentHeat() { return SL_latent_heat_; }

    //! Return first time derivative of radii \f$\dot{R}_{n}\f$ or \f$\dot{R}_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> RadiusDot() { return radiusdot_; }

    //! Return mass of particles
    virtual Teuchos::RCP<const Epetra_Vector> Mass() { return mass_; }

    //! Return inertia of particles
    virtual Teuchos::RCP<const Epetra_Vector> Inertia() { return inertia_; }

    //! Return initial density
    virtual const double initDensity() { return initDensity_; }

    //! Return initial radius
    virtual const double initRadius() { return initRadius_; }

    //! extract heat source - constant pressure - solid state
    virtual const double CPS() { return CPS_; }

    //! extract heat source - constant pressure - liquid state
    virtual const double CPL() { return CPL_; }

    //! Return specific latent heat - solid <-> liquid
    virtual const double MaxSpecificLatentHeat() {return SL_latent_heat_max_; }

    //! transition temperature - solid <-> liquid
    virtual const double SL_TransitionTemperature() {return SL_transitionTemperature_; }

    //! does collision handler exist?
    virtual bool HaveCollHandler() {return (collhandler_!=Teuchos::null); }


    //@}

    //! @name Extract methods
    //@{

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp(){ return disn_; }

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp(){ return veln_; }

    //! Return accelerations
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp() { return accn_; }

    //! extract temperatures
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessTemperaturenp(){ return temperaturen_; }

    //! extract pressure
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessPressurenp(){ return pressuren_; }

    //! Return density
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDensitynp() { return densityn_; }

    //! extract SL latent heat
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSL_latentHeat(){ return SL_latent_heat_; }

    //! extract initial radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius0() { return radius0_; }

    //! extract radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadius() { return radius_; }

    //! extract radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessLatentHeat() { return SL_latent_heat_; }

    //! extract fist time derivative of radii
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessRadiusDot() { return radiusdot_; }

    //! extract mass vector
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessMass() { return mass_; }

    //! extract inertia vector
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessInertia() { return inertia_; }

    //@}

    //! Set forces due to interface loads, the force is expected external-force-like
    virtual void SetForceInterface
    (
      Teuchos::RCP<Epetra_MultiVector> iforce  ///< the force on interface
    );

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

    //@}

  protected:

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int myrank_;  //!< ID of actual processor in parallel
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                  //!< containing non-overlapping
                                                  //!< map of global DOFs on Dirichlet
                                                  //!< boundary conditions
    //@}

    //! @name Printing and output
    //@{
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    bool printlogo_;  //!< true: enjoy your cuppa
    int printscreen_;  //!< print infos to standard out every n steps
    FILE* errfile_;  //!< error file handle
    bool printerrfile_;  //!< print infos to error file #errfile_
    int writerestartevery_;  //!< write restart every given step;
                             //!< if 0, restart is not written
    bool writestate_;  //!< write state on/off
    bool writevelacc_;  //!< write velocity and acceleration on/off
    int writeresultsevery_;  //!< write state/stress/strain every given step
    int writeenergyevery_;  //!< write system energy every given step
    Teuchos::RCP<std::ofstream> energyfile_;  //!< outputfile for energy
    bool writeorientation_;  //!< write orientation of particle on/off

    double kinergy_; //!< kinetic energy
    double intergy_; //!< internal energy
    double extergy_; //!< external energy
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

    //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > dis_;
    //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > vel_;
    //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > acc_;
    //! global temperatures \f${T}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > temperature_;
    //! global pressure \f${P}_{n+1}\f$ at \f$p_{n+1}\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > pressure_;
    //! global pressure \f${D}_{n+1}\f$ at \f$d_{n+1}\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > density_;
    //! global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disn_;
    //! global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;
    //! global accelerations \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;
    //! global temperatures \f${T}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>  temperaturen_;
    //! global pressure \f${P}_{n+1}\f$ at \f$p_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>  pressuren_;
    //! global density \f${D}_{n+1}\f$ at \f$d_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>  densityn_;
    //! radius of each particle at \f$t_{n}\f$ or \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> SL_latent_heat_;
    //! radius of each particle at \f$t_{0}\f$
    Teuchos::RCP<Epetra_Vector> radius_;
    //! latent heat of each particle at \f$t_{n}\f$ or \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> radius0_;
    //! first derivative of time of radius of each particle at \f$t_{n}\f$ or $\f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> radiusdot_;
    //! mass of each particle at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> mass_;
    //! moment of inertia of each particle at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> inertia_;
    //! global angular-velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ang_vel_;
    //! global angular-accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ang_acc_;
    //! global angular-velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ang_veln_;
    //! global angular-acceleration \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ang_accn_;
    //! global vector for visualization of particle orientation
    Teuchos::RCP<Epetra_Vector> orient_;
    //! initial density
    double initDensity_;
    //! initial radius
    double initRadius_;
    //! specific heat - constant pressure - solid state
    double CPS_;
    //! specific heat - constant pressure - liquid state
    double CPL_;
    //! specific latent heat - solid <-> liquid
    double SL_latent_heat_max_;
    //! transition temperature - solid <-> liquid
    double SL_transitionTemperature_;
    //! volumetric thermal expansion coefficient - solid state - usually 3 times the linear coefficient
    double S_thermalExpansion_;
    //! volumetric thermal expansion coefficient - liquid state
    double L_thermalExpansion_;
    //! volumetric thermal expansion coefficient - solid <-> liquid - Beware! It is linked to latent heat instead of temperature!
    double SL_thermalExpansion_;

    //! additional external forces
    Teuchos::RCP<Epetra_Vector> fifc_;


    //@}

    //! @name control flags for special features
    //@{

    const bool variableradius_;   //!< bool to control whether radius is constant or variable

    //@}

    //! @name helper classes
    //@{

    //! particle algorithm that is connected to the time integration
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //! collision handler for particles
    Teuchos::RCP<PARTICLE::ParticleCollisionHandlerBase> collhandler_;

    //@}

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_H
