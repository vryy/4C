/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all 0D cardiovascular model terms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_STRUCTURE_NEW_MODEL_EVALUATOR_HPP
#define FOUR_C_CARDIOVASCULAR0D_STRUCTURE_NEW_MODEL_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d_manager.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Cardiovascular0D : public Generic
    {
     public:
      //! constructor
      Cardiovascular0D();


      void Setup() override;

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_cardiovascular0d; }

      //! derived
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

      //! derived
      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void DetermineEnergy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void ResetStepState() override;

      //! [derived]
      void PostOutput() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

     private:
      Teuchos::RCP<UTILS::Cardiovascular0DManager> cardvasc0dman_;  //!< Cardiovascular0D manager

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! cardio contributions to the structural stiffness matrix
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_cardio_ptr_;

      //! structural rhs contributions of the cardio model at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fstructcardio_np_ptr_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
