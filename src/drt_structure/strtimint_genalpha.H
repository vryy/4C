/*----------------------------------------------------------------------*/
/*!
\file strtimint_genalpha.H
\brief Structural time integration with generalised-alpha

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_GENALPHA_H
#define STRTIMINT_GENALPHA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_impl.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Generalised-alpha time integration
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

    //! Zienkiewicz-Xie auxiliary scheme is friend
    friend class TimAdaZienXie;

  public:

    //! verify if given coefficients are in admissible range;
    //! prints also info to STDOUT
    void VerifyCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntGenAlpha
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
    // ....

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const
    {
      return INPAR::STR::dyna_genalpha;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps()
    {
      return 1;
    }

    //! Give linear order of accuracy of displacement part
    int MethodOrderOfAccuracyDis()
    {
      return (fabs(MethodLinErrCoeffDis2()) < EPS6) ? 3 : 2;
    }

    //! Give linear order of accuracy of velocity part
    int MethodOrderOfAccuracyVel()
    {
      return (fabs(MethodLinErrCoeffVel1()) < EPS6) ? 2 : 1;
    }

    //! return time integration factor
    virtual const double TimIntParam() const { return alphaf_; }

    //! return Newmark parameter \f$beta\f$
    virtual const double GetNewmarkBeta() const { return beta_; }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis()
    {
      if (MethodOrderOfAccuracyDis() == 2)
        return MethodLinErrCoeffDis2();
      else
        return MethodLinErrCoeffDis3();
    }

    //! 2nd order linear error coefficient of displacements
    double MethodLinErrCoeffDis2()
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0/6.0 - beta_ + alphaf_/2.0 - alpham_/2.0;
    }

    //! 3rd order linear error coefficient of displacements
    double MethodLinErrCoeffDis3()
    {
      // at least true for am<1/2 and large enough n->infty
      return 1./24.
             - beta_/2.*(1.-2*alphaf_+2.*alpham_)
             - 1./4.*(alphaf_-alpham_)*(1.-2.*alpham_);
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel()
    {
      if (MethodOrderOfAccuracyVel() == 1)
        return MethodLinErrCoeffVel1();
      else
        return MethodLinErrCoeffVel2();
    }

    //! 1st order linear error coefficient of velocities
    double MethodLinErrCoeffVel1()
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0/2.0 - gamma_ + alphaf_ - alpham_;
    }

    //! 2nd order linear error coefficient of velocities
    double MethodLinErrCoeffVel2()
    {
      // at least true for am<1/2 and large enough n->infty
      return 1./6.
             - gamma_/2.*(1.-2*alphaf_+2.*alpham_)
             - 1./2.*(alphaf_-alpham_)*(1.-2.*alpham_);
    }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void PredictConstDisConsistVelAcc();

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! \author mayr.mt
    void PredictConstVelConsistAcc();

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! \author mayr
    void PredictConstAcc();

    //! Create force residual #fres_ and its stiffness #stiff_
    void EvaluateForceStiffResidual(bool predict = false);

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    void EvaluateForceStiffResidualRelax();

    //! Determine characteristic norm for displacements
    //! \author lw (originally)
    double CalcRefNormDisplacement();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_
    //! from scratch by using the newly updated #disn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstDisConsistVelAcc)
    void UpdateIterIncrementally();

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void UpdateIterIteratively();

    //! Update step
    void UpdateStepState();

    //! Update element
    void UpdateStepElement();

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);

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

    virtual Teuchos::RCP<Epetra_Vector> Velaf() { return velm_; }

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
    enum INPAR::STR::MidAverageEnum midavg_;
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
    Teuchos::RCP<Epetra_Vector> dism_;  //!< mid-displacements
                                        //!< \f$D_m = D_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> velm_;  //!< mid-velocities
                                        //!< \f$V_m = V_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> accm_;  //!< mid-accelerations
                                        //!< \f$A_m = A_{n+1-\alpha_m}\f$
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

    Teuchos::RCP<Epetra_Vector> finert_;  //!< inertia force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> finertm_;  //!< inertia mid-force
    Teuchos::RCP<Epetra_Vector> finertn_;  //!< inertia force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fviscm_;  //!< viscous force

    //@}

  };  // class TimIntGenAlpha

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_GENALPHA_H
