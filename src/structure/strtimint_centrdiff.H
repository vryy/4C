/*----------------------------------------------------------------------*/
/*! \file
\brief Structural time integration with central difference scheme 2nd order (explicit)
\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_CENTRDIFF_H
#define STRTIMINT_CENTRDIFF_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{
  /*====================================================================*/
  /*!
   * \brief Central differences: 2nd order accurate,
   *                             explicit time integrator,
   *                             1-step method
   *
   *  The displacement update reads
   *  \f[
   *  d_{n+1} = d_n + \Delta t \dot{d}_n + \frac{1}{2}\Delta t^2 \ddot{d}_n
   *  \f]
   *
   *  The velocity update reads
   *  \f[
   *  \dot{d}_{n+1} = \dot{d}_n
   *                + \frac{1}{2}\Delta t \ddot{d}_n
   *                + \frac{1}{2}\Delta t \ddot{d}_{n+1}
   *  \f]
   *
   * \author popp
   * \date 01/12
   */
  class TimIntCentrDiff : public TimIntExpl
  {
   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntCentrDiff(const Teuchos::ParameterList& timeparams,  //!< time params
        const Teuchos::ParameterList& ioparams,                //!< ioflags
        const Teuchos::ParameterList& sdynparams,              //!< input parameters
        const Teuchos::ParameterList& xparams,                 //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,              //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,                   //!< the solver
        Teuchos::RCP<LINALG::Solver> contactsolver,            //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output          //!< the output
    );

    //! Empty constructor
    TimIntCentrDiff() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntCentrDiff(const TimIntCentrDiff& old) : TimIntExpl(old) { ; }

    //! Destructor
    virtual ~TimIntCentrDiff() { ; }

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
    virtual void Init(const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver);

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
    enum INPAR::STR::DynamicType MethodName() const { return INPAR::STR::dyna_centrdiff; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() const { return 1; }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const { return 2; }

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const { return 2; }

    /*! Return linear error coefficient of displacements
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx \frac{1}{6}\Delta t_n^3 \dddot{d_n} + HOT(\Delta t_n^4)
     *  \f]
     */
    virtual double MethodLinErrCoeffDis() const { return 1. / 6.; }

    /*! Return linear error coefficient of velocities
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx -\frac{1}{12}\Delta t_n^3 d_n^{(4)} + HOT(\Delta t_n^4)
     *  \f]
     */
    virtual double MethodLinErrCoeffVel() const { return -1. / 12.; }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() { return fextn_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      dserror("FextNew() not available in CentrDiff");
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
    Teuchos::RCP<Epetra_Vector> fextn_;   //!< external force
                                          //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fintn_;   //!< internal force
                                          //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fviscn_;  //!< Rayleigh viscous forces
                                          //!< \f$C \cdot V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fcmtn_;   //!< contact or meshtying forces
                                          //!< \f$F_{cmt;n+1}\f$
    Teuchos::RCP<Epetra_Vector> frimpn_;  //!< time derivative of
                                          //!< linear momentum
                                          //!< (temporal rate of impulse)
                                          //!< \f$\dot{P}_{n+1} = M \cdot \dot{V}_{n+1}\f$
    //@}

  };  // class TimIntCentrDiff

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_CENTRDIFF_H
