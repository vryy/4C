/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all spring dashpot terms


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_SPRINGDASHPOT_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_SPRINGDASHPOT_HPP

#include "4C_config.hpp"

#include "4C_constraint_springdashpot.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Solid
{
  namespace ModelEvaluator
  {
    class SpringDashpot : public Generic
    {
     public:
      //! constructor
      SpringDashpot();

      void setup() override;

      //! derived
      Inpar::Solid::ModelType type() const override { return Inpar::Solid::model_springdashpot; }

      //! derived
      void reset(const Epetra_Vector& x) override;

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override{};

      //! derived
      void post_evaluate() override{};

      //! derived
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void predict(const Inpar::Solid::PredEnum& pred_type) override{};

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override{};

      //! derived
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
      }

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override{};

      //! derived
      void update_step_state(const double& timefac_n) override;

      //! derived
      void update_step_element() override{};

      //! derived
      void determine_stress_strain() override{};

      //! derived
      void determine_energy() override{};

      //! derived
      void determine_optional_quantity() override{};

      //! derived
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void reset_step_state() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! [derived]
      void post_output() override;

     private:
      //! all spring dashpot instances
      std::vector<Teuchos::RCP<CONSTRAINTS::SpringDashpot>> springs_;

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! structural velocity at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> velnp_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_spring_ptr_;

      //! spring forces at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fspring_np_ptr_;
    };

  }  // namespace ModelEvaluator
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
