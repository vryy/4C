/*--------------------------------------------------------------------------*/
/*! \file

\brief solution algorithm for stationary problems

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_TIMINT_STAT_HPP
#define FOUR_C_LUBRICATION_TIMINT_STAT_HPP


#include "4C_config.hpp"

#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace LUBRICATION
{
  class TimIntStationary : public TimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void Init() override;

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    /// read restart data
    void read_restart(int step) override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #prenp_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #predict_const_temp_consist_rate)
    void update_iter_incrementally() override;

   protected:
    /// don't want = operator and cctor
    TimIntStationary operator=(const TimIntStationary& old);

    /// copy constructor
    TimIntStationary(const TimIntStationary& old);

    /// set time parameter for element evaluation
    void set_element_time_parameter() const override;

    //! set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    /// add actual Neumann loads with time factor
    void add_neumann_to_residual() override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const override { return 1.0; }

  };  // class TimIntStationary

}  // namespace LUBRICATION

FOUR_C_NAMESPACE_CLOSE

#endif
