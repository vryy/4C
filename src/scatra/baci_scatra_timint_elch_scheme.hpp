/*----------------------------------------------------------------------*/
/*! \file

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2



*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_ELCH_SCHEME_HPP
#define FOUR_C_SCATRA_TIMINT_ELCH_SCHEME_HPP

#include "baci_config.hpp"

#include "baci_scatra_timint_bdf2.hpp"
#include "baci_scatra_timint_elch.hpp"
#include "baci_scatra_timint_elch_scl.hpp"
#include "baci_scatra_timint_genalpha.hpp"
#include "baci_scatra_timint_ost.hpp"
#include "baci_scatra_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN

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

    void Update() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override;

   protected:
    void WriteRestart() const override;

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

    void Update() override;

   protected:
    void ComputeTimeDerivPot0(const bool init) override{};

    void ElectrodeKineticsTimeUpdate() override{};

    void ExplicitPredictor() const override;

    void SetOldPartOfRighthandside() override;

    void WriteRestart() const override;
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

    void Update() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override{};

   protected:
    void WriteRestart() const override;

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

    void Update() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override;

   protected:
    void WriteRestart() const override;

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

    void Update() override;

    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    void PreCalcInitialPotentialField() override;

    void PostCalcInitialPotentialField() override{};

   protected:
    void WriteRestart() const override;

    void ElectrodeKineticsTimeUpdate() override
    {
      FOUR_C_THROW(
          "Galvanostatic-BC is not implemented for the stationary time-integration scheme");
    };

    void ComputeTimeDerivPot0(const bool init) override;
  };
}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
