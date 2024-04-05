/*----------------------------------------------------------------------*/
/*! \file
\brief Structural time integration with generalised energy-momentum method
\level 1
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_TIMINT_GEMM_HPP
#define FOUR_C_STRUCTURE_TIMINT_GEMM_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_structure_timint_impl.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{
  /*====================================================================*/
  /*!
   * \brief Generalised energy-momentum method (GEMM)
   *
   * <h3>References</h3>
   * - [1] D Kuhl and MA Crisfield, Energy-conserving and decaying algorithms in non-linear
   * structural dynamics, International Journal for Numerical Methods in Engineering, 45(5):569-599,
   * 1999.
   * - [2] JC Simo and N Tarnow, The discrete energy-momentum method, Conserving algorithms for
   * nonlinear elastodynamics, Journal of Applied Mathematics and Physics, 43:757-792, 1992.
   *
   * \author bborn
   * \date 07/08
   */
  class TimIntGEMM : public TimIntImpl
  {
    //! Zienkiewicz-Xie auxiliary scheme is friend
    friend class TimAdaZienXie;

   public:
    //! @name Lifetime
    //@{

    //! Constructor
    TimIntGEMM(const Teuchos::ParameterList& timeparams,   //!< time parameters
        const Teuchos::ParameterList& ioparams,            //!< ioflags
        const Teuchos::ParameterList& sdynparams,          //!< input parameters
        const Teuchos::ParameterList& xparams,             //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,          //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,         //!< the solver
        Teuchos::RCP<CORE::LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output      //!< the output
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

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() override { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const override { return INPAR::STR::dyna_gemm; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() const override { return 1; }

    //! Give linear order of accuracy of displacement part
    int MethodOrderOfAccuracyDis() const override
    {
      return (fabs(MethodLinErrCoeffDis2()) < 1e-6) ? 3 : 2;
    }

    //! Give linear order of accuracy of velocity part
    int MethodOrderOfAccuracyVel() const override
    {
      return (fabs(MethodLinErrCoeffVel1()) < 1e-6) ? 2 : 1;
    }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const override
    {
      if (MethodOrderOfAccuracyDis() == 2)
        return MethodLinErrCoeffDis2();
      else
        return MethodLinErrCoeffDis3();
    }

    //! 2nd order linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis2() const
    {
      // at least true for am<1/2 and large enough n->infty and xi=0
      return 1.0 / 6.0 - beta_ + alphaf_ / 2.0 - alpham_ / 2.0;
    }

    //! 3rd order linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis3() const
    {
      // at least true for am<1/2 and large enough n->infty and xi=0
      return 1. / 24. - beta_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 4. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const override
    {
      if (MethodOrderOfAccuracyVel() == 1)
        return MethodLinErrCoeffVel1();
      else
        return MethodLinErrCoeffVel2();
    }

    //! 1st order linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel1() const
    {
      // at least true for am<1/2 and large enough n->infty and xi=0
      return 1.0 / 2.0 - gamma_ + alphaf_ - alpham_;
    }

    //! 2nd order linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel2() const
    {
      // at least true for am<1/2 and large enough n->infty and xi=0
      return 1. / 6. - gamma_ / 2. * (1. - 2 * alphaf_ + 2. * alpham_) -
             1. / 2. * (alphaf_ - alpham_) * (1. - 2. * alpham_);
    }

    //! return time integration factor
    double TimIntParam() const override { return alphaf_; }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void PredictConstDisConsistVelAcc() override;

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! \author mayr
    void PredictConstVelConsistAcc() override;

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! \author mayr.mt
    void PredictConstAcc() override;

    //! Create force residual #fres_ and its stiffness #stiff_
    void EvaluateForceStiffResidual(Teuchos::ParameterList& params) final;

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    void EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params) override;

    //! Evaluate residual #fres_
    void EvaluateForceResidual() override;

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce() override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_ from scratch by using the
    //! newly updated #disn_.
    //! This method is necessary for certain predictors (like #PredictConstDisConsistVelAcc)
    void UpdateIterIncrementally() override;

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void UpdateIterIteratively() override;

    //! Update step
    void UpdateStepState() override;

    //! Update Element
    void UpdateStepElement() override;

    //! Read and set restart for forces
    void ReadRestartForce() override;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) override;

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() override { return fextn_; }

    //@}

    //! @name GEMM-specific methods
    //@{

    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();

    //! Evaluate ordinary internal force, its stiffness at mid-state
    void ApplyForceStiffInternalMid(const double time,    //!< evaluation time
        const double dt,                                  //!< step size
        const Teuchos::RCP<Epetra_Vector> dis,            //!< displacement state at t_n
        const Teuchos::RCP<Epetra_Vector> disn,           //!< displacement state at t_{n+1}
        const Teuchos::RCP<Epetra_Vector> disi,           //!< residual displacements
        const Teuchos::RCP<Epetra_Vector> vel,            // velocity state
        Teuchos::RCP<Epetra_Vector> fint,                 //!< internal force
        Teuchos::RCP<CORE::LINALG::SparseOperator> stiff  //!< stiffness matrix
    );

    //! Evaluate ordinary internal force at mid-state
    void ApplyForceInternalMid(const double time,  //!< evaluation time
        const double dt,                           //!< step size
        const Teuchos::RCP<Epetra_Vector> dis,     //!< displacement state at t_n
        const Teuchos::RCP<Epetra_Vector> disn,    //!< displacement state at t_{n+1}
        const Teuchos::RCP<Epetra_Vector> disi,    //!< residual displacements
        const Teuchos::RCP<Epetra_Vector> vel,     // velocity state
        Teuchos::RCP<Epetra_Vector> fint           //!< internal force
    );

    //@}

   protected:
    //! equal operator is NOT wanted
    TimIntGEMM operator=(const TimIntGEMM& old);

    //! copy constructor is NOT wanted
    TimIntGEMM(const TimIntGEMM& old);

    //! @name Key coefficients
    //@{
    double beta_;    //!< factor (0,0.5]
    double gamma_;   //!< factor (0,1]
    double alphaf_;  //!< factor [0,1)
    double alpham_;  //!< factor [0,1)
    double xi_;      //!< factor [0,1)
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
    // Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
    // Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;   //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> finertm_;  //!< inertia force
    Teuchos::RCP<Epetra_Vector> fviscm_;   //!< viscous force

    //@}

  };  // class TimIntGEMM

}  // namespace STR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_TIMINT_GEMM_H
