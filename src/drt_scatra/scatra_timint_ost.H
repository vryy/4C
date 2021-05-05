/*----------------------------------------------------------------------*/
/*! \file
\brief One-Step-Theta time-integration scheme

\level 1


*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_OST_H
#define SCATRA_TIMINT_OST_H

#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

namespace SCATRA
{
  class TimIntOneStepTheta : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<LINALG::Solver> solver,                      //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,              //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,         //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,            //!< output writer
        const int probnum = 0                                     //!< global problem number
    );

    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old) = delete;

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old) = delete;

    /// initialize time integration scheme
    void Init() override;

    /// setup time integration scheme
    void Setup() override;

    //! operator for manipulations before call to \ref Solve() ; May be overridden by subclass.
    void PreSolve() override{};

    //! operator for manipulations after call to \ref Solve() ; May be overridden by subclass.
    void PostSolve() override{};

    /// Print information about current time step to screen (reimplementation for OST)
    void PrintTimeStepInfo() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() override{};

    /// compute values at the interior of the elements (required for hdg)
    void ComputeInteriorValues() override{};

    ///  compute scalar time derivative
    void ComputeTimeDerivative() override;

    ///  compute scalar time derivate parameters of the input voltage
    void ComputeTimeDerivPot0(const bool init) override{};

    /// update the solution after convergence of the nonlinear iteration.
    /// current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    /// routine to return scalar field phi at time step n+alpha_F
    Teuchos::RCP<Epetra_Vector> Phiaf() override { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phiam() override { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phidtam() override { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi at time step n+1
    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() override
    {
      Teuchos::RCP<Teuchos::ParameterList> timeparams;
      timeparams = Teuchos::rcp(new Teuchos::ParameterList());
      timeparams->set("using stationary formulation", false);
      timeparams->set("using generalized-alpha time integration", false);
      timeparams->set("total time", time_);
      timeparams->set("time factor", theta_ * dta_);
      timeparams->set("alpha_F", 1.0);
      return timeparams;
    }

    //! set state on micro scale in multi-scale simulations
    void SetState(Teuchos::RCP<Epetra_Vector> phin,  //!< micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector> phinp,           //!< micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector>
            phidtn,  //!< time derivative of micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector>
            phidtnp,  //!< time derivative of micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector> hist,               //!< micro-scale history vector
        Teuchos::RCP<IO::DiscretizationWriter> output,  //!< micro-scale discretization writer
        const std::vector<double>&
            phinp_macro,   //!< values of state variables at macro-scale Gauss point
        const int step,    //!< time step
        const double time  //!< time
    );

    //! clear state on micro scale in multi-scale simulations
    void ClearState();

    void PreCalcInitialTimeDerivative() override;

    void PostCalcInitialTimeDerivative() override;

   protected:
    /// set time parameter for element evaluation
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
    void CalcInitialTimeDerivative() override;

    /// set part of residual vector belonging to previous time step
    void SetOldPartOfRighthandside() override;

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    void ExplicitPredictor() const override;

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual() override;

    /// AVM3-based scale separation
    void AVM3Separation() override;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override;

    // dynamic Vreman model
    void DynamicComputationOfCv() override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    /// write additional data required for restart
    void OutputRestart() const override;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const override { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;

  };  // class TimIntOneStepTheta

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
