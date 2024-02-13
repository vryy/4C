/*----------------------------------------------------------------------*/
/*! \file
\brief Generalized-alpha time-integration scheme

\level 1


*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_GENALPHA_HPP
#define BACI_SCATRA_TIMINT_GENALPHA_HPP

#include "baci_config.hpp"

#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_scatra_timint_implicit.hpp"

BACI_NAMESPACE_OPEN

namespace SCATRA
{
  class TimIntGenAlpha : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntGenAlpha(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<CORE::LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Setup() override;

    void PrintTimeStepInfo() override
    {
      if (myrank_ == 0)
      {
        printf(
            "\nTIME: %11.4E/%11.4E  DT = %11.4E  %s(a_F=%3.2f | a_M=%3.2f | gamma=%3.2f) STEP = "
            "%4d/%4d\n",
            time_, maxtime_, dta_, MethodTitle().c_str(), alphaF_, alphaM_, gamma_, step_,
            stepmax_);
      }
    }

    void ComputeIntermediateValues() override;

    void ComputeInteriorValues() override{};

    void ComputeTimeDerivative() override;

    void ComputeTimeDerivPot0(const bool init) override{};

    void Update() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    Teuchos::RCP<Epetra_Vector> Phiaf() override { return phiaf_; }

    Teuchos::RCP<Epetra_Vector> Phiafnp() override { return phiaf_; }

    Teuchos::RCP<Epetra_Vector> Phiam() override { return phiam_; }

    Teuchos::RCP<Epetra_Vector> Phidtam() override { return phidtam_; }

    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phiaf_, *fsphiaf_);
      return fsphiaf_;
    };

    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() override
    {
      Teuchos::RCP<Teuchos::ParameterList> timeparams;
      timeparams = Teuchos::rcp(new Teuchos::ParameterList());
      timeparams->set("using stationary formulation", false);
      timeparams->set("using generalized-alpha time integration", true);
      timeparams->set("total time", time_ - (1 - alphaF_) * dta_);
      timeparams->set("time factor", genalphafac_ * dta_);
      timeparams->set("alpha_F", alphaF_);
      return timeparams;
    }

    void PreCalcInitialTimeDerivative() override;

    void PostCalcInitialTimeDerivative() override;


   protected:
    /// don't want = operator and cctor
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    /// copy constructor
    TimIntGenAlpha(const TimIntGenAlpha& old);

    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    void SetElementTimeParameterBackwardEuler() const override;

    void CalcInitialTimeDerivative() override;

    void SetOldPartOfRighthandside() override;

    void ExplicitPredictor() const override;

    void AddNeumannToResidual() override;

    void AVM3Separation() override;

    void DynamicComputationOfCs() override;

    void DynamicComputationOfCv() override;

    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    void WriteRestart() const override;

    double ResidualScaling() const override { return 1.0 / (dta_ * genalphafac_); }

    /// scalar at time n+alpha_F and n+alpha_M
    Teuchos::RCP<Epetra_Vector> phiaf_;
    Teuchos::RCP<Epetra_Vector> phiam_;

    /// scalar time derivative at time n+alpha_M
    Teuchos::RCP<Epetra_Vector> phidtam_;

    /// fine-scale part at time n+alpha_F
    Teuchos::RCP<Epetra_Vector> fsphiaf_;

    /// time factors for generalized-alpha time integration
    double alphaM_;
    double alphaF_;
    double gamma_;
    double genalphafac_;
  };  // class TimIntGenAlpha

}  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif  // BACI_SCATRA_TIMINT_GENALPHA_H
