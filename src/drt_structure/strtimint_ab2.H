/*----------------------------------------------------------------------*/
/*!
\file strtimint_ab2.H
\brief Structural time integration with Adams-Bashforth 2nd order (explicit)
\level 1

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
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
      const Teuchos::ParameterList& timeparams,   //!< time parameters
      const Teuchos::ParameterList& ioparams,     //!< ioflags
      const Teuchos::ParameterList& sdynparams,   //!< input parameters
      const Teuchos::ParameterList& xparams,      //!< extra flags
      //const Teuchos::ParameterList& ab2params,  //!< AB2 flags
      Teuchos::RCP<DRT::Discretization> actdis,   //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,        //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    TimIntAB2() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntAB2(const TimIntAB2& old) : TimIntExpl(old) { ; }

    //! Destructor
    virtual ~TimIntAB2() { ; };

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    virtual void Init(
        const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<LINALG::Solver> solver );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    //@}

    //! @name Actions
    //@{

    //! Resize \p TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next time step.
    void UpdateStepState();

    //! Update Element
    void UpdateStepElement();

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum INPAR::STR::DynamicType MethodName() const
    {
      return INPAR::STR::dyna_ab2;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() const { return 2; }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const { return 2; }

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const { return 2; }

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const
    {
      const double dt = (*dt_)[0];
      const double dto = (*dt_)[-1];
      return (2.*dt + 3.*dto)/(12.*dt);
    }

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const
    {
      return MethodLinErrCoeffDis();
    }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      return fextn_;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      dserror("FextNew() not available in AB2");
      return Teuchos::null;
    }

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);

    //@}


  protected:

    //! @name Global forces at \f$t_{n+1}\f$
    //@{
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fviscn_; //!< Rayleigh viscous forces
                                         //!< \f$C \cdot V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fcmtn_;  //!< contact or meshtying forces
                                         //!< \f$F_{cmt;n+1}\f$
    Teuchos::RCP<Epetra_Vector> frimpn_; //!< time derivative of
                                         //!< linear momentum
                                         //!< (temporal rate of impulse)
                                         //!< \f$\dot{P}_{n+1} = M \cdot \dot{V}_{n+1}\f$
    //@}

  };  // class TimIntAB2

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_AB2_H
