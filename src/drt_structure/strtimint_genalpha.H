/*----------------------------------------------------------------------*/
/*!
\brief Structural time integration with generalised-alpha

\level 1

\maintainer Matthias Mayr
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

    //! calculate coefficients from given spectral radius
    void CalcCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntGenAlpha(const Teuchos::ParameterList& timeparams,  //!< time params
        const Teuchos::ParameterList& ioparams,               //!< ioflags
        const Teuchos::ParameterList& sdynparams,             //!< input parameters
        const Teuchos::ParameterList& xparams,                //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,             //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,                  //!< the solver
        Teuchos::RCP<LINALG::Solver> contactsolver,           //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output         //!< the output
    );

    //! Destructor
    // ....

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

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const { return INPAR::STR::dyna_genalpha; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() const { return 1; }

    //! Give linear order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const
    {
      return (fabs(MethodLinErrCoeffDis2()) < EPS6) ? 3 : 2;
    }

    //! Give linear order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const
    {
      return (fabs(MethodLinErrCoeffVel1()) < EPS6) ? 2 : 1;
    }

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const
    {
      if (MethodOrderOfAccuracyDis() == 2)
        return MethodLinErrCoeffDis2();
      else
        return MethodLinErrCoeffDis3();
    }

    //! 2nd order linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis2() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0 / 6.0 - beta_ + alphaf_ / 2.0 - alpham_ / 2.0;
    }

    //! 3rd order linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis3() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1. / 24. - beta_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 4. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const
    {
      if (MethodOrderOfAccuracyVel() == 1)
        return MethodLinErrCoeffVel1();
      else
        return MethodLinErrCoeffVel2();
    }

    //! 1st order linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel1() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0 / 2.0 - gamma_ + alphaf_ - alpham_;
    }

    //! 2nd order linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel2() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1. / 6. - gamma_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 2. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! return time integration factor
    virtual double TimIntParam() const { return alphaf_; }

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
    void EvaluateForceStiffResidual(Teuchos::ParameterList& params) final;

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    void EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params);

    //! Evaluate residual #fres_
    void EvaluateForceResidual();

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
    Teuchos::RCP<Epetra_Vector> Fext() { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() { return fextn_; }

    //@}

    //! @name Generalised-alpha specific methods
    //@{
    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();
    //@}

    //! Build total residual vector and effective tangential stiffness
    //! matrix in case of nonlinear, rotational inertia effects
    void BuildResStiffNLMassRot(Teuchos::RCP<Epetra_Vector> fres_,
        Teuchos::RCP<Epetra_Vector> fextn_, Teuchos::RCP<Epetra_Vector> fintn_,
        Teuchos::RCP<Epetra_Vector> finertn_, Teuchos::RCP<LINALG::SparseOperator> stiff_,
        Teuchos::RCP<LINALG::SparseOperator> mass_);

    //! Check, if there are solely beam elements in the whole discretization
    bool SolelyBeam3Elements(Teuchos::RCP<DRT::Discretization> actdis);

   protected:
    //! equal operator is NOT wanted
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    //! copy constructor is NOT wanted
    TimIntGenAlpha(const TimIntGenAlpha& old);

    //! @name set-up
    //@{
    //! mid-average type more at MidAverageEnum
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
    double beta_;     //!< factor (0,1/2]
    double gamma_;    //!< factor (0,1]
    double alphaf_;   //!< factor [0,1)
    double alpham_;   //!< factor [-1,1)
    double rho_inf_;  //!< factor[0,1]
    //@}

    //! @name Global mid-state vectors
    //@{

    //! mid-displacements \f$D_m = D_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> dism_;

    //! mid-velocities \f$V_m = V_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> velm_;

    //! mid-accelerations \f$A_m = A_{n+1-\alpha_m}\f$
    Teuchos::RCP<Epetra_Vector> accm_;

    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;   //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;   //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> finert_;   //!< inertia force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> finertm_;  //!< inertia mid-force
    Teuchos::RCP<Epetra_Vector> finertn_;  //!< inertia force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fviscm_;  //!< viscous force

    Teuchos::RCP<Epetra_Vector> fint_str_;  //!< pure structural global internal force at \f$t_n\f$
                                            //!< i.e. no condensation of EAS,pasticity,...
    //@}

  };  // class TimIntGenAlpha

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_GENALPHA_H
