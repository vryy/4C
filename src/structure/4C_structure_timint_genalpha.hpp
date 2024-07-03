/*----------------------------------------------------------------------*/
/*! \file
\brief Structural time integration with generalised-alpha

\level 1

*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_TIMINT_GENALPHA_HPP
#define FOUR_C_STRUCTURE_TIMINT_GENALPHA_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_structure_timint_impl.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace Solid
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
        Teuchos::RCP<Core::FE::Discretization> actdis,        //!< current discretisation
        Teuchos::RCP<Core::LinAlg::Solver> solver,            //!< the solver
        Teuchos::RCP<Core::LinAlg::Solver> contactsolver,     //!< the solver for contact meshtying
        Teuchos::RCP<Core::IO::DiscretizationWriter> output   //!< the output
    );

    //! Destructor
    // ....

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    void init(const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams, Teuchos::RCP<Core::FE::Discretization> actdis,
        Teuchos::RCP<Core::LinAlg::Solver> solver) override;

    /*! \brief Setup all class internal objects and members

     setup() is not supposed to have any input arguments !

     Must only be called after init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    void setup() override;

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() override { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum Inpar::Solid::DynamicType MethodName() const override
    {
      return Inpar::Solid::dyna_genalpha;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() const override { return 1; }

    //! Give linear order of accuracy of displacement part
    int method_order_of_accuracy_dis() const override
    {
      return (fabs(method_lin_err_coeff_dis2()) < 1e-6) ? 3 : 2;
    }

    //! Give linear order of accuracy of velocity part
    int method_order_of_accuracy_vel() const override
    {
      return (fabs(method_lin_err_coeff_vel1()) < 1e-6) ? 2 : 1;
    }

    //! Return linear error coefficient of displacements
    double method_lin_err_coeff_dis() const override
    {
      if (method_order_of_accuracy_dis() == 2)
        return method_lin_err_coeff_dis2();
      else
        return method_lin_err_coeff_dis3();
    }

    //! 2nd order linear error coefficient of displacements
    virtual double method_lin_err_coeff_dis2() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0 / 6.0 - beta_ + alphaf_ / 2.0 - alpham_ / 2.0;
    }

    //! 3rd order linear error coefficient of displacements
    virtual double method_lin_err_coeff_dis3() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1. / 24. - beta_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 4. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! Return linear error coefficient of velocities
    double method_lin_err_coeff_vel() const override
    {
      if (method_order_of_accuracy_vel() == 1)
        return method_lin_err_coeff_vel1();
      else
        return method_lin_err_coeff_vel2();
    }

    //! 1st order linear error coefficient of velocities
    virtual double method_lin_err_coeff_vel1() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1.0 / 2.0 - gamma_ + alphaf_ - alpham_;
    }

    //! 2nd order linear error coefficient of velocities
    virtual double method_lin_err_coeff_vel2() const
    {
      // at least true for am<1/2 and large enough n->infty
      return 1. / 6. - gamma_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 2. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! return time integration factor
    double TimIntParam() const override
    {
      check_is_init();
      return alphaf_;
    }

    //! return time integration factor-alpham
    virtual double TimIntParamAlpham() const
    {
      check_is_init();
      return alpham_;
    }

    //! return time integration factor-beta
    virtual double TimIntParamBeta() const
    {
      check_is_init();
      return beta_;
    }

    //! return time integration factor-gamma
    virtual double TimIntParamGamma() const
    {
      check_is_init();
      return gamma_;
    }


    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void predict_const_dis_consist_vel_acc() override;

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! \author mayr.mt
    void predict_const_vel_consist_acc() override;

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! \author mayr
    void PredictConstAcc() override;

    //! Create force residual #fres_ and its stiffness #stiff_
    void evaluate_force_stiff_residual(Teuchos::ParameterList& params) final;

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with solve_relaxation_linear
    void evaluate_force_stiff_residual_relax(Teuchos::ParameterList& params) override;

    //! Evaluate residual #fres_
    void evaluate_force_residual() override;

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce() override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_ from scratch by using the
    //! newly updated #disn_.
    //! This method is necessary for certain predictors (like #predict_const_dis_consist_vel_acc)
    void update_iter_incrementally() override;

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void update_iter_iteratively() override;

    //! Update step
    void UpdateStepState() override;

    //! Update element
    void UpdateStepElement() override;

    //! Read and set restart for forces
    void ReadRestartForce() override;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<Core::IO::DiscretizationWriter> output) override;

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() override { return fextn_; }

    //@}

    //! @name Generalised-alpha specific methods
    //@{
    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();
    //@}

    //! Build total residual vector and effective tangential stiffness
    //! matrix in case of nonlinear, rotational inertia effects
    void build_res_stiff_nl_mass_rot(Teuchos::RCP<Epetra_Vector> fres_,
        Teuchos::RCP<Epetra_Vector> fextn_, Teuchos::RCP<Epetra_Vector> fintn_,
        Teuchos::RCP<Epetra_Vector> finertn_, Teuchos::RCP<Core::LinAlg::SparseOperator> stiff_,
        Teuchos::RCP<Core::LinAlg::SparseOperator> mass_);

    //! Check, if there are solely beam elements in the whole discretization
    bool SolelyBeam3Elements(Teuchos::RCP<Core::FE::Discretization> actdis);

   protected:
    //! equal operator is NOT wanted
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    //! copy constructor is NOT wanted
    TimIntGenAlpha(const TimIntGenAlpha& old);

    //! @name set-up
    //@{
    //! mid-average type more at MidAverageEnum
    enum Inpar::Solid::MidAverageEnum midavg_;
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

}  // namespace Solid

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
