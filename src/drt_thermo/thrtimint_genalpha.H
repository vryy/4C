/*----------------------------------------------------------------------*/
/*!
\file thrtimint_genalpha.H
\brief Thermal time integration with generalised-alpha
\level 1
\maintainer Christoph Meier
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 05/13 |
 *----------------------------------------------------------------------*/
#ifndef THRTIMINT_GENALPHA_H
#define THRTIMINT_GENALPHA_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 05/13 |
 *----------------------------------------------------------------------*/
#include "thrtimint_impl.H"


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                     dano 05/13 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  /*!
   * \brief Generalised-alpha time integration
   *
   * References
   * - [1] J Chung and GM Hulbert, A time integration algorithm for structural
   *   dynamics with improved numerical dissipation: the generalized-alpha method
   *   Journal of Applied Mechanics, 60:371-375, 1993.
   *
   * temporal discretisation of a first order ODE according to
   * - [2] KE Jansen, CH Whiting and GM Hulbert, A generalized-alpha
   *   method for integrating the filtered Navier-Stokes equations with a
   *   stabilized finite element method, Computer Methods in Applied Mechanics
   *   and Engineering, 190:305-319, 2000.
   *
   *
   * \author danowski
   * \date 06/13
   */
  class TimIntGenAlpha : public TimIntImpl
  {
   public:
    //! verify if given coefficients are in admissable range
    //! prints also info to STDOUT
    void VerifyCoeff();

    //! calculate coefficients for given spectral radius
    void CalcCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntGenAlpha(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,           //!< input parameters
        const Teuchos::ParameterList& xparams,              //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,           //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,                //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output       //!< the output
    );

    //! Destructor
    // ....

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::THR::DynamicType MethodName() const { return INPAR::THR::dyna_genalpha; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() { return 1; }

    //! Give linear order of accuracy of temperature part
    int MethodOrderOfAccuracy() { return (fabs(MethodLinErrCoeff()) < EPS6) ? 2 : 1; }

    // TODO 2013-07-05 check the calculation of the factor again
    //! Return linear error coefficient of temperatures
    double MethodLinErrCoeff()
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0 / 2.0 - gamma_ + alphaf_ - alpham_;
    }

    //! Consistent predictor with constant temperatures
    //! and consistent temperature rates and temperatures
    void PredictConstTempConsistRate();

    //! Evaluate ordinary internal force, its tangent at state
    void ApplyForceTangInternal(const double time,  //!< evaluation time
        const double dt,                            //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,     //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,    //!< residual temperatures
        Teuchos::RCP<Epetra_Vector> fcap,           //!< capacity force
        Teuchos::RCP<Epetra_Vector> fint,           //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang     //!< tangent matrix
    );

    //! Evaluate ordinary internal force
    void ApplyForceInternal(const double time,    //!< evaluation time
        const double dt,                          //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,   //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
        Teuchos::RCP<Epetra_Vector> fint          //!< internal force
    );

    //! Evaluate a convective boundary condition
    // (nonlinear --> add term to tangent)
    void ApplyForceExternalConv(const double time,  //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> tempn,    //!< temperature state T_n
        const Teuchos::RCP<Epetra_Vector> temp,     //!< temperature state T_n+1
        Teuchos::RCP<Epetra_Vector> fext,           //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang     //!< tangent matrix
    );

    //! Create force residual #fres_ and its tangent #tang_
    void EvaluateRhsTangResidual();

    //! Determine characteristic norm for temperatures
    //! \author lw (originally)
    double CalcRefNormTemperature();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #tempn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    void UpdateIterIncrementally();

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #tempn_ and #raten_ by
    //! incrementing these vector proportional to the residual
    //! temperatures #tempi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual temperatures #tempi_ are blanked at these DOFs.
    void UpdateIterIteratively();

    //! Update step
    void UpdateStepState();

    //! Update Element
    void UpdateStepElement();

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() { return fextn_; }

    //@}

    //! @name Generalised-alpha specific methods
    //@{
    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();
    //@}

   protected:
    //! equal operator is NOT wanted
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    //! copy constructor is NOT wanted
    TimIntGenAlpha(const TimIntGenAlpha& old);

    //! @name set-up
    //@{
    //! mid-average type more at #MidAverageEnum
    enum INPAR::THR::MidAverageEnum midavg_;
    //@}

    //! @name Key coefficients
    //! Please note, to obtain a second-order accurate scheme, you need
    //! to follow the following formulas in which \f$\rho_\infty\f$ is the
    //! spectral radius.
    //! \f[ \alpha_m = 1/2 (3 - \rho_\infty)/(\rho_\infty + 1) \f]
    //! \f[ \alpha_f = 1/(\rho_\infty + 1) \f]
    //! \f[ \gamma = 1/2 + \alpha_m - \alpha_f \f]
    //! The spectral radius is responsible for the magnitude of numerical
    //! dissipation introduced.
    //! For instance
    //!
    //! Without numerical dissipation at \f$\rho_\infty=1\f$: 2nd order mid-point rule
    //! \f[ \alpha_f=0.5, \alpha_m=0.5, \gamma=0.5 \f]
    //!
    //! Strong numerical dissipation at \f$\rho_\infty=0.5\f$: default
    //! \f[ \alpha_f=2/3, \alpha_m=5/6, \gamma=2/3 \f]
    //!
    //! Maximal numerical dissipation at \f$\rho_\infty=0.0\f$: BDF2
    //! \f[ \alpha_f=1, \alpha_m=3/2, \gamma=1 \f]
    //@{
    double gamma_;    //!< factor (0,1]
    double alphaf_;   //!< factor [0,1]
    double alpham_;   //!< factor [0,1.5]
    double rho_inf_;  //!< factor[0,1]
    //@}

    //! @name Global mid-state vectors
    //@{
    Teuchos::RCP<Epetra_Vector> tempm_;  //!< mid-temperatures
                                         //!< \f$T_m = T_{n+\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> ratem_;  //!< mid-temperature rates
                                         //!< \f$R_m = R_{n+\alpha_m}\f$
    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;   //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force at \f$t_{n+\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;   //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force \f$t_{n+\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fcap_;  //!< capacity force \f$C\cdot\Theta_n\f$ at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector>
        fcapm_;  //!< capacity force \f$C\cdot\Theta_{n+\alpha_m}\f$ at \f$t_{n+\alpha_m}\f$
    Teuchos::RCP<Epetra_Vector>
        fcapn_;  //!< capacity force \f$C\cdot\Theta_{n+1}\f$ at \f$t_{n+1}\f$

    //@}

  };  // class TimIntGenAlpha

}  // namespace THR


#endif  // #ifndef THRTIMINT_GENALPHA_H
