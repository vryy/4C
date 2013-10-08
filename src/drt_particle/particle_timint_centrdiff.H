/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.H
\brief Particle time integration with central difference scheme 2nd order (explicit),
              also known as Velocity-Verlet algorithm

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_CENTRDIFF_H
#define PARTICLE_TIMINT_CENTRDIFF_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_structure/strtimint_centrdiff.H"


//forward declarations
namespace PARTICLE
{
  class Algorithm;
  class ParticleCollisionHandler;
}

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Central differences: 2nd order accurate,
   *                             explicit time integrator,
   *                             1-step method
   * \author hammerl
   * \date 09/12
   */
  class TimIntCentrDiff : public STR::TimIntCentrDiff
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntCentrDiff
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntCentrDiff(){ ; }

    //! Empty constructor
    TimIntCentrDiff(){ ; }

    //! Copy constructor
    TimIntCentrDiff(const TimIntCentrDiff& old){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    virtual void DetermineMassDampConsistAccel();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Nothing to update on element level
    void UpdateStepElement(){ return; }

    //! State vectors are updated according to the new distribution of particles
    void UpdateStatesAfterParticleTransfer();

    //! Write restart which is equivalent to standard output for particles
    virtual void OutputRestart(bool& datawritten);

    //! output displacements, velocities and accelerations
    virtual void OutputState(bool& datawritten);

    //! Read and set restart state
    virtual void ReadRestartState();

    //! Read restart forces
    virtual void ReadRestartForce() { return; }

    //! particle algorithm is set from outside
    void SetParticleAlgorithm(Teuchos::RCP<PARTICLE::Algorithm> algorithm) { particle_algorithm_ = algorithm; }

    //@}

    //! @name Vector access
    //@{

    //! Return radii of particles
    Teuchos::RCP<Epetra_Vector> ExtractRadius() { return radius_; }

    //! Return mass of particles
    Teuchos::RCP<Epetra_Vector> ExtractMass() { return mass_; }

    //! Return density
    double ParticleDensity() { return density_; }

    //@}

  protected:

    //! apply initial condition for particle radius
    void ApplyInitialRadiusCondition();

    //! acceleration is applied from given forces
    void ComputeAcc(
      Teuchos::RCP<Epetra_Vector> f_contact,
      Teuchos::RCP<Epetra_Vector> m_contact,
      Teuchos::RCP<Epetra_Vector> acc,
      Teuchos::RCP<Epetra_Vector> ang_acc
      );

    //! update step
    virtual void UpdateStepState();

    //! states are given to the collision handler
    void SetStatesForCollision();

    //! initialization of vector for visualization of the particle orientation
    void InitializeOrientVector();

    //! update of vector for visualization of the particle orientation
    void RotateOrientVector(double dt);
    
    //! energy of the particle system is calculated
    virtual void DetermineEnergy();

    //! Attach file handle for energy file #energyfile_
    virtual void AttachEnergyFile(std::string name = "");

    //! @name Global properties at \f$t_{n}\f$ and \f$t_{n+1}\f$
    //@{

    Teuchos::RCP<Epetra_Vector> radius_;  //!< radius of each particle at \f$t_{n}\f$

    Teuchos::RCP<Epetra_Vector> mass_;   //!< mass of each particle at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> inertia_;   //!< moment of inertia of each particle at \f$t_{n+1}\f$

    Teuchos::RCP<STR::TimIntMStep<Epetra_Vector> > ang_vel_;   //!< global angular-velocities \f${V}_{n}, V_{n-1}, ...\f$

    Teuchos::RCP<STR::TimIntMStep<Epetra_Vector> > ang_acc_;   //!< global angular-accelerations \f${A}_{n}, A_{n-1}, ...\f$

    Teuchos::RCP<Epetra_Vector> ang_veln_;   //!< global angular-velocitie \f${V}_{n+1}\f$ at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> ang_accn_;   //!< global angular-acceleration \f${A}_{n+1}\f$ at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> orient_;   //!< global vector for visualization of particle orientation

    double density_; //!< density of each particle at \f$t_{n+1}\f$

    //@}

    //! @name helper classes
    //@{

    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;   //!< particle algorithm that is connected to the time integration

    Teuchos::RCP<PARTICLE::ParticleCollisionHandler> collhandler_;   //!< collision handler for particles

    //@}

  };  // class TimIntCentrDiff

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_CENTRDIFF_H
