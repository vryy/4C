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

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Nothing to update on element level
    void UpdateStepElement(){ return; }

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    //! Make damping and mass matrix --> not yet needed
    virtual void DetermineMassDampConsistAccel();

    //! State vectors are updated according to the new distribution of particles
    void UpdateStatesAfterParticleTransfer();

    //! Write restart which is equivalent to standard output for particles
    virtual void OutputRestart(bool& datawritten);

    //! output displacements, velocities and accelerations
    virtual void OutputState(bool& datawritten);

    //! Read and set restart state
    virtual void ReadRestartState();

    //! Read restart forces
    virtual void ReadRestartForce();

    //@}

    //! @name Vector access
    //@{

    //! Return radii of particles
    Teuchos::RCP<Epetra_Vector> ExtractRadiusnp() { return radiusn_; }

    //! Return density
    double ParticleDensity() { return density_; }

    //@}


  protected:

    //! @name Global properties at \f$t_{n}\f$ and \f$t_{n+1}\f$
    //@{
    Teuchos::RCP<Epetra_Vector> radius_;  //!< radius of each particle at \f$t_{n}\f$

    Teuchos::RCP<Epetra_Vector> radiusn_;  //!< radius of each particle at \f$t_{n+1}\f$

    double density_; //!< density of each particle at \f$t_{n+1}\f$

    //@}


  };  // class TimIntCentrDiff

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_CENTRDIFF_H
