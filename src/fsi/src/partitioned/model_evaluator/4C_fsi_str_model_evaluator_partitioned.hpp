/*-----------------------------------------------------------*/
/*! \file

\brief Model evaluator for structure part of partitioned fsi


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FSI_STR_MODEL_EVALUATOR_PARTITIONED_HPP
#define FOUR_C_FSI_STR_MODEL_EVALUATOR_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class Structure;
}  // namespace ADAPTER

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedFSI : public Generic
    {
     public:
      //! constructor
      PartitionedFSI();


      //! get pointer to force vector at time level n+1 (full structural map).
      //! interface part is inserted in ADAPTER::FSIStructureWrapper.
      const Teuchos::RCP<Epetra_Vector>& get_interface_force_np_ptr()
      {
        return interface_force_np_ptr_;
      };

      //! setup class variables [derived]
      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //! @{

      //! [derived]
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_partitioned_coupling; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override { return; };

      //! [derived]
      bool evaluate_force() override { return true; };

      //! [derived]
      bool evaluate_stiff() override { return true; };

      //! [derived] not needed in partitioned scheme
      bool evaluate_force_stiff() override { return true; };

      //! derived
      void pre_evaluate() override { return; };

      //! derived
      void post_evaluate() override { return; };

      //! derived
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override
      {
        return true;
      };

      //! [derived]
      void write_restart(
          CORE::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override
      {
        return;
      };

      //! [derived]
      void read_restart(CORE::IO::DiscretizationReader& ioreader) override { return; };

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! recover condensed Lagrange multipliers
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
        return;
      };

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! [derived]
      void UpdateStepState(const double& timefac_n) override;

      //! [derived]
      void UpdateStepElement() override { return; };

      //! [derived]
      void determine_stress_strain() override { return; };

      //! [derived]
      void DetermineEnergy() override { return; };

      //! [derived]
      void determine_optional_quantity() override { return; };

      //! [derived]
      void OutputStepState(CORE::IO::DiscretizationWriter& iowriter) const override { return; };

      //! derived
      void ResetStepState() override { return; };

      //! [derived]
      void PostOutput() override { return; };

      //! [derived]
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! @}

      /*! \brief linear structure solve with just a interface load
       *
       * The very special solve done in steepest descent relaxation
       * calculation (and matrix free Newton Krylov).
       *
       * \note Can only be called after a valid structural solve. */
      Teuchos::RCP<const Epetra_Vector> solve_relaxation_linear(
          Teuchos::RCP<ADAPTER::Structure> structure);

      virtual void setup_multi_map_extractor();

      //! set flag true if in order to request a relaxation solve
      void set_is_relaxation_solve(bool trueorfalse) { is_relaxationsolve_ = trueorfalse; };

      //! Returns the global input/output data container
      const STR::TIMINT::BaseDataIO& GetInOutput() const;

     private:
      //! fsi interface force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> interface_force_np_ptr_;

      //! true if relaxation solve is requested
      bool is_relaxationsolve_;

    };  // class PartitionedFSI

  }  // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
