/*---------------------------------------------------------------------------*/
/*! \file
\brief model evaluator for structure part of partitioned pasi
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PASI_STR_MODEL_EVALUATOR_PARTITIONED_HPP
#define FOUR_C_PASI_STR_MODEL_EVALUATOR_PARTITIONED_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace ADAPTER
{
  class Structure;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedPASI : public Generic
    {
     public:
      //! constructor
      PartitionedPASI();

      //! setup class variables [derived]
      void Setup() override;

      //! get pointer to force vector at time level n+1 (full structural map)
      //! interface part is inserted in ADAPTER::PASIStructureWrapper
      const Teuchos::RCP<Epetra_Vector>& get_interface_force_np_ptr()
      {
        return interface_force_np_ptr_;
      };

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

      //! [derived]
      void pre_evaluate() override { return; };

      //! [derived]
      void post_evaluate() override { return; };

      //! [derived]
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override
      {
        return true;
      };

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override
      {
        return;
      };

      //! [derived]
      void read_restart(IO::DiscretizationReader& ioreader) override { return; };

      //! [derived]
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! recover condensed Lagrange multipliers
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
        return;
      };

      //! [derived]
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override { return; };

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
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override { return; };

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! [derived]
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

     private:
      //! pasi interface force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> interface_force_np_ptr_;
    };

  }  // namespace MODELEVALUATOR

}  // namespace STR

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
