/*----------------------------------------------------------------------*/
/*! \file
\brief One-Step-Theta time-integration scheme

\level 1


*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_OST_H
#define SCATRA_TIMINT_OST_H

#include "scatra_timint_implicit.H"

#include "linalg_utils_sparse_algebra_create.H"

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

    void Init() override;

    void Setup() override;

    void PreSolve() override{};

    void PostSolve() override{};

    void PrintTimeStepInfo() override;

    void ComputeIntermediateValues() override{};

    void ComputeInteriorValues() override{};

    void ComputeTimeDerivative() override;

    void ComputeTimeDerivPot0(const bool init) override{};

    void Update(const int num = 0) override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    Teuchos::RCP<Epetra_Vector> Phiaf() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> Phiam() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> Phidtam() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    }

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
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    void CalcInitialTimeDerivative() override;

    void SetOldPartOfRighthandside() override;

    void ExplicitPredictor() const override;

    void AddNeumannToResidual() override;

    void AVM3Separation() override;

    void DynamicComputationOfCs() override;

    void DynamicComputationOfCv() override;

    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    void OutputRestart() const override;

    double ResidualScaling() const override { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;
  };

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
