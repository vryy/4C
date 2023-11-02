/*----------------------------------------------------------------------*/
/*! \file
\brief Structural time integration with forward Euler (explicit)
\level 2
*/

/*----------------------------------------------------------------------*/
#ifndef BACI_STRUCTURE_TIMINT_EXPLEULER_H
#define BACI_STRUCTURE_TIMINT_EXPLEULER_H

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_structure_timint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{
  /*====================================================================*/
  /*!
   * \brief forward Euler: 1st order accurate,
   *                       explicit time integrator,
   * \author bborn
   * \date 06/08
   */
  class TimIntExplEuler : public TimIntExpl
  {
   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntExplEuler(const Teuchos::ParameterList& timeparams,  //!< time params
        const Teuchos::ParameterList& ioparams,                //!< ioflags
        const Teuchos::ParameterList& sdynparams,              //!< input parameters
        const Teuchos::ParameterList& xparams,                 //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,              //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,             //!< the solver
        Teuchos::RCP<CORE::LINALG::Solver> contactsolver,      //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output          //!< the output
    );

    //! Empty constructor
    TimIntExplEuler() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntExplEuler(const TimIntExplEuler& old) : TimIntExpl(old) { ; }

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
    void Init(const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<CORE::LINALG::Solver> solver) override;

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
    void Setup() override;

    //@}

    //! @name Actions
    //@{

    //! Resize \p TimIntMStep<T> multi-step quantities
    void ResizeMStep() override;

    //! Do time integration of single step
    int IntegrateStep() override;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next time step.
    void UpdateStepState() override;

    //! Update Element
    void UpdateStepElement() override;

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum INPAR::STR::DynamicType MethodName() const override { return INPAR::STR::dyna_expleuler; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() const override { return 1; }

    //! Give local order of accuracy of displacement part
    int MethodOrderOfAccuracyDis() const override { return 1; }

    //! Give local order of accuracy of velocity part
    int MethodOrderOfAccuracyVel() const override { return 1; }

    /*! \brief Return linear error coefficient of displacements
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx \frac{1}{2}\Delta t_n^2 \ddot{d_n} + HOT(\Delta t_n^3)
     *  \f]
     */
    double MethodLinErrCoeffDis() const override { return 0.5; }

    /*! \brief Return linear error coefficient of velocities
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx \frac{1}{2}\Delta t_n^2 \dddot{d_n} + HOT(\Delta t_n^3)
     *  \f]
     */
    double MethodLinErrCoeffVel() const override { return 0.5; }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override { return fextn_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() override
    {
      dserror("FextNew() not available in AB2");
      return Teuchos::null;
    }

    //! Read and set restart for forces
    void ReadRestartForce() override;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) override;

    //@}


   protected:
    bool modexpleuler_;  //!< modified explicit Euler equation (veln_ instead of vel_ for calc of
                         //!< disn_), default: true

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

  };  // class TimIntExplEuler

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // STRUCTURE_TIMINT_EXPLEULER_H
