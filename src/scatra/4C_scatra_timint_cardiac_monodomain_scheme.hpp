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


namespace ScaTra
{
  class TimIntCardiacMonodomainOST : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainOST(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

   protected:
    void write_restart() const override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainOST


  class TimIntCardiacMonodomainBDF2 : public virtual TimIntCardiacMonodomain,
                                      public virtual TimIntBDF2
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainBDF2(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

   protected:
    void write_restart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainBDF2


  class TimIntCardiacMonodomainGenAlpha : public virtual TimIntCardiacMonodomain,
                                          public virtual TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainGenAlpha(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! Setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

   protected:
    void write_restart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainGenAlpha
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
