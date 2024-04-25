/*----------------------------------------------------------------------*/
/*! \file

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        Cardiac-monodomain-specific implementation (class TimIntCardiacMonodomain)

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPP
#define FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_bdf2.hpp"
#include "4C_scatra_timint_cardiac_monodomain.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_ost.hpp"

FOUR_C_NAMESPACE_OPEN


namespace SCATRA
{
  class TimIntCardiacMonodomainOST : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainOST(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    void WriteRestart() const override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    //! do not calculate initial scalar time derivatives for ep
    void CalcInitialTimeDerivative() override { return; };

  };  // class TimIntCardiacMonodomainOST


  class TimIntCardiacMonodomainBDF2 : public virtual TimIntCardiacMonodomain,
                                      public virtual TimIntBDF2
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainBDF2(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    void WriteRestart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void CalcInitialTimeDerivative() override { return; };

  };  // class TimIntCardiacMonodomainBDF2


  class TimIntCardiacMonodomainGenAlpha : public virtual TimIntCardiacMonodomain,
                                          public virtual TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! Setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

   protected:
    void WriteRestart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void CalcInitialTimeDerivative() override { return; };

  };  // class TimIntCardiacMonodomainGenAlpha
}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
