/*----------------------------------------------------------------------*/
/*! \file

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2



*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_ELCH_SCHEME_H
#define BACI_SCATRA_TIMINT_ELCH_SCHEME_H

#include "baci_scatra_timint_ost.H"
#include "baci_scatra_timint_bdf2.H"
#include "baci_scatra_timint_genalpha.H"
#include "baci_scatra_timint_stat.H"
#include "baci_scatra_timint_elch.H"
#include "baci_scatra_timint_elch_scl.H"

namespace SCATRA
{
  class ScaTraTimIntElchOST : public ScaTraTimIntElch, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchOST(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Init() override;

    void Setup() override;

    void Update(const int num = 0) override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override;

   protected:
    void OutputRestart() const override;

    void ElectrodeKineticsTimeUpdate() override;

    void ExplicitPredictor() const override;

    void ComputeTimeDerivPot0(const bool init) override;

    void SetOldPartOfRighthandside() override;
  };

  class ScaTraTimIntElchSCLOST : public ScaTraTimIntElchSCL, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchSCLOST(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    void Init() override;

    void PostCalcInitialPotentialField() override;

    void PreCalcInitialPotentialField() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void Setup() override;

    void Update(const int num = 0) override;

   protected:
    void ComputeTimeDerivPot0(const bool init) override{};

    void ElectrodeKineticsTimeUpdate() override{};

    void ExplicitPredictor() const override;

    void SetOldPartOfRighthandside() override;

    void OutputRestart() const override;
  };


  class ScaTraTimIntElchBDF2 : public ScaTraTimIntElch, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchBDF2(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Init() override;

    void Setup() override;

    void Update(const int num = 0) override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override{};

   protected:
    void OutputRestart() const override;

    void ElectrodeKineticsTimeUpdate() override;

    void ComputeTimeDerivPot0(const bool init) override;

    void SetOldPartOfRighthandside() override;
  };


  class ScaTraTimIntElchGenAlpha : public ScaTraTimIntElch, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Init() override;

    void Setup() override;

    void Update(const int num = 0) override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override;

   protected:
    void OutputRestart() const override;

    void ElectrodeKineticsTimeUpdate() override;

    void ComputeTimeDerivPot0(const bool init) override;
  };


  class ScaTraTimIntElchStationary : public ScaTraTimIntElch, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Init() override;

    void Setup() override;

    void Update(const int num = 0) override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override{};

   protected:
    void OutputRestart() const override;

    void ElectrodeKineticsTimeUpdate() override
    {
      dserror("Galvanostatic-BC is not implemented for the stationary time-integration scheme");
    };

    void ComputeTimeDerivPot0(const bool init) override;
  };
}  // namespace SCATRA

#endif  // SCATRA_TIMINT_ELCH_SCHEME_H
