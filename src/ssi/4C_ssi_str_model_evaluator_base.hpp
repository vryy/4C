/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for scalar-structure interaction

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SSI_STR_MODEL_EVALUATOR_BASE_HPP
#define FOUR_C_SSI_STR_MODEL_EVALUATOR_BASE_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid::MODELEVALUATOR
{
  class BaseSSI : public Generic
  {
   public:
    bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override { return true; }

    bool assemble_jacobian(
        Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override
    {
      return true;
    }

    void determine_energy() override {}

    void determine_optional_quantity() override {}

    void determine_stress_strain() override;

    bool evaluate_force() override { return true; }

    bool evaluate_force_stiff() override { return true; }

    bool evaluate_stiff() override { return true; }

    [[nodiscard]] Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override
    {
      FOUR_C_THROW("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override
    {
      FOUR_C_THROW("Not implemented!");
      return Teuchos::null;
    }

    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> get_mechanical_stress_state() const override
    {
      return mechanical_stress_state_;
    }

    void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override {}

    void post_evaluate() override {}

    void post_output() override {}

    void predict(const Inpar::Solid::PredEnum& pred_type) override {}

    void pre_evaluate() override {}

    void read_restart(Core::IO::DiscretizationReader& ioreader) override {}

    void reset(const Epetra_Vector& x) override {}

    void reset_step_state() override {}

    void run_post_compute_x(
        const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
    {
    }

    void run_post_iterate(const ::NOX::Solver::Generic& solver) override {}

    void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
        const NOX::Nln::Group& curr_grp) override
    {
    }

    void setup() override;

    [[nodiscard]] Inpar::Solid::ModelType type() const override
    {
      return Inpar::Solid::model_basic_coupling;
    }

    void update_step_element() override {}

    void update_step_state(const double& timefac_n) override {}

    void write_restart(
        Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override
    {
    }

   private:
    //! mechanical stress state
    Teuchos::RCP<Epetra_Vector> mechanical_stress_state_;
  };
}  // namespace Solid::MODELEVALUATOR
FOUR_C_NAMESPACE_CLOSE

#endif