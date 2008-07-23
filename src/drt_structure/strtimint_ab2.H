/*----------------------------------------------------------------------*/
/*!
\file strtimint_ab2.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_AB2_H
#define STRTIMINT_AB2_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Adams-Bashforth2: 2nd order accurate, 
   *                          explicit time integrator,
   *                          linear 2-step method
   * \author bborn
   * \date 06/08
   */
  class TimIntAB2 : public TimIntExpl
  {

    //! A friendly relation to adaptive partner
    friend class TimAdaAB2;

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntAB2
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      //const Teuchos::ParameterList& ab2params,  //!< AB2 flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    //! Empty constructor
    TimIntAB2() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntAB2(const TimIntAB2& old) : TimIntExpl(old) { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    void ResizeMStep();

    //@}

    //! @name Actions
    //@{

    //! Do time integration of single step
    void IntegrateStep();

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step 
    //! becomes impossible. We are ready and keen awating the next time step.
    void UpdateStep();

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum TimInt::NameEnum MethodName() { return TimInt::name_ab2; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() { return 2; }

    //! Give local order of accuracy of displacement part
    int MethodOrderOfAccuracyDis() { return 2; }

    //! Give local order of accuracy of velocity part
    int MethodOrderOfAccuracyVel() { return 2; }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis()
    {
      const double dt = (*dt_)[0];
      const double dto = (*dt_)[-1];
      return (2.*dt + 3.*dto)/(12.*dt);
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel()
    {
      return MethodLinErrCoeffDis();
    }

    //@}

  protected:

    //! Initialise
    void Init();

    //! @name Global forces at \f$t_{n+1}\f$
    //@{
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fviscn_;  //!< Rayleigh viscous forces
                                          //!< \f$C \cdot V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> frimpn_;  //!< time derivative of
                                          //!< linear momentum
                                          //!< (temporal rate of impulse)
                                          //!< \f$\dot{P}_{n+1} = M \cdot \dot{V}_{n+1}\f$
    //@}

  };  // class TimIntAB2

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_AB2_H
#endif  // #ifdef CCADISCRET
