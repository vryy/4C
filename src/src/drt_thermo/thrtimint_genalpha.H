/*----------------------------------------------------------------------*/
/*!
\file thrtimint_genalpha.H

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef THRTIMINT_GENALPHA_H
#define THRTIMINT_GENALPHA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "thrtimint_impl.H"

/*----------------------------------------------------------------------*/
/* belongs to thermal dynamics namespace */
namespace THR
{

  /*====================================================================*/
  /*!
   * \brief Generalized-alpha time integration
   *
   * References
   * - [1] NM Newmark, A method of computation for structural dynamics,
   *   Journal of Earthquake Mechanic Division, 85(EM3):67-94, 1959.
   * - [2] HM Hilber, TJR Hughes and RL Taylor, Improved numerical
   *   dissipation for time integration algorithms, Earthquake
   *   Engineering and Structural Mechanics, 5:283-292, 1977.
   *
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntGenAlpha : public TimIntImpl
  {

//    //! Zienkiewicz-Xie auxiliar scheme is friend
//    friend class TimAdaZienXie;

  public:

    //! verify if given coefficients are in admissable range;
    //! prints also info to STDOUT
    void VerifyCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntGenAlpha
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& tdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
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
    enum INPAR::THR::DynamicType MethodName() const
    {
      return INPAR::THR::dyna_genalpha;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() { return 1; }

    //! Give linear order of accuracy of temperature part
    int MethodOrderOfAccuracy()
    {
//      // Ueberpruefe diesen Term!!!!!!!!!!! 02.10.2009
//      //return fabs(1./2.-theta_)<EPS10 ? 2 : 1;
      return 1;
    }
//
//    //! Return linear error coefficient of temperatures
    double MethodLinErrCoeff()
    {
//      // Ueberpruefe diesen Term!!!!!!!!!!! 02.10.2009
//      //return 1./2. - theta_;
      return 1./2.;
    }

    //! Consistent predictor with constant temperatures
    //! and consistent temperature rates and temperatures
    void PredictConstTempConsistRate();

    //! Evaluate ordinary internal force, its tangent at state
    void ApplyForceTangInternal(
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
      const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
      Teuchos::RCP<LINALG::SparseMatrix> tang  //!< tangent matrix
    );

    //! Evaluate ordinary internal force
    void ApplyForceInternal(
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
      const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
      Teuchos::RCP<Epetra_Vector> fint  //!< internal force
    );

    //! Create force residual #fres_ and ist tangent #tang_
    void EvaluateRhsTangResidual();

    //! Determine characteristic norm for temperatures
    //! \author lw (originally)
    double CalcRefNormTemperature();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_ and #acc_
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

    //! Read and set restart for forces
    void ReadRestartForce();
    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      return fext_;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      return fextn_;
    }

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
    //! \f[ \alpha_m = (2*\rho_\infty - 1)/(\rho_\infty + 1) \f]
    //! \f[ \alpha_f = \rho_\infty/(\rho_\infty + 1) \f]
    //! \f[ \beta = 1/4*(1 - \alpha_m + \alpha_f)^2 \mbox{(max. damp. of high-freq. modes)} \f]
    //! \f[ \gamma = 1/2 - \alpha_m + \alpha_f \f]
    //! The spectral radius is responsible for the magnitude of
    //! numerical dissipation introduced.
    //! For instance
    //! Without numerical dissipation at \f$\rho_\infty=1\f$
    //! \f[ \beta=0.25, \gamma=0.5, \alpha_m=0.5, \alpha_f=0.5 \f]
    //! Medium dissipation at \f$\rho_\infty=0.8\f$
    //! \f[ \beta=25/81, \gamma=11/18, \alpha_m=1/3, \alpha_f=4/9  \f]
    //! Strong numerical dissipation at \f$\rho_\infty=0.5\f$
    //! \f[ \beta=4/9, \gamma=10/12, \alpha_m=0, \alpha_f=1/3 \f]
    //@{
    double beta_;  //!< factor (0,1/2]
    double gamma_;  //!< factor (0,1]
    double alphaf_;  //!< factor [0,1)
    double alpham_;  //!< factor [0,1)
    //@}

    //! @name Global mid-state vectors
    //@{
    Teuchos::RCP<Epetra_Vector> tempm_;  //!< mid-temperatures
                                        //!< \f$T_m = T_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> ratem_;  //!< mid-temperature rates
                                        //!< \f$R_m = R_{n+1-\alpha_f}\f$
    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    //@}

  };  // class TimIntGenAlpha

}  // namespace THR


#endif  // #ifndef THRTIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
