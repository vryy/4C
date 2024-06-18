/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all constraint terms


\level 3
*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_LAGPENCONSTRAINT_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_LAGPENCONSTRAINT_HPP

#include "4C_config.hpp"

#include "4C_constraint_manager.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace LAGPENCONSTRAINT
{
  class NoxInterface;
  class NoxInterfacePrec;
}  // namespace LAGPENCONSTRAINT

namespace STR
{
  namespace MODELEVALUATOR
  {
    class LagPenConstraint : public Generic
    {
     public:
      //! constructor
      LagPenConstraint();


      void setup() override;

      //! derived
      Inpar::STR::ModelType Type() const override { return Inpar::STR::model_lag_pen_constraint; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override;

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override { return; };

      //! derived
      void post_evaluate() override { return; };

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
      void Predict(const Inpar::STR::PredEnum& pred_type) override { return; };

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void update_step_state(const double& timefac_n) override;

      //! derived
      void update_step_element() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void determine_energy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void reset_step_state() override;

      //! [derived]
      void post_output() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      const Teuchos::RCP<CONSTRAINTS::ConstrManager>& StrategyPtr();

      //! Return the NOX::Nln::CONSTRAINT::Interface::Required member object
      const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface>& nox_interface_ptr();

      //! Return the NOX::Nln::CONSTRAINT::Interface::Preconditioner member object
      const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec>& NoxInterfacePrecPtr();

     protected:
      //! Returns the underlying contact strategy object
      CONSTRAINTS::ConstrManager& strategy();
      const CONSTRAINTS::ConstrManager& strategy() const;

     private:
      //! all constraint instances
      Teuchos::RCP<CONSTRAINTS::ConstrManager> constrman_;  //!< Constraint manager

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_constr_ptr_;

      //! constraint contributions to the structural rhs at \f%t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fstrconstr_np_ptr_;

      //! pointer to the NOX::Nln::CONSTRAINT::Interface::Required object
      Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface> noxinterface_ptr_;

      //! pointer to the NOX::Nln::CONSTRAINT::Interface::Preconditioner object
      Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec> noxinterface_prec_ptr_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
