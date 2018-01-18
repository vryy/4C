/*----------------------------------------------------------------------*/
/*!
\file particle_timint_rk.H
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

\level 3
<pre>
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_RK_H
#define PARTICLE_TIMINT_RK_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Runge-Kutta:         2nd/4th order accurate,
   *                             explicit time integrator
   * \author rasthofer
   * \date 01/14
   */
  class TimIntRK : public PARTICLE::TimIntExpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntRK
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntRK(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of level-set vectors
    virtual void Init();

    //! do time integration of single time step
    virtual int IntegrateStep();

    //! state vectors are updated according to the new distribution of particles
    virtual void UpdateStatesAfterParticleTransfer();

    //! write restart
    virtual void OutputRestart
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! output displacements, velocities and accelerations
    //! and more system vectors
    virtual void OutputState
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! read and set restart state
    virtual void ReadRestartState();

    //@}

    //! @name Access methods
    //@{

    //! return sign of particles
    virtual Teuchos::RCP<const Epetra_Vector> Sign() { return sign_; }

    //@}

    //! @name Extract methods
    //@{

    //! return sign of particles
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessSign() { return sign_; }

    //@}

  protected:

    //! Runge-Kutta scheme 2nd order
    void Integrate_RK_Second();

    //! Runge-Kutta scheme 4th order
    void Integrate_RK_Fourth();

  private:

    //! output displacement
    virtual void OutputDisplacement() const;

    //! @name General method parameters
    //@{

    INPAR::PARTICLE::DynamicType rk_scheme_;

    //@}

    //! @name Global state vectors
    //@{

    //! sign of each particle
    Teuchos::RCP<Epetra_Vector> sign_;

    //@}

  };  // class TimIntRK

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_RK_H
