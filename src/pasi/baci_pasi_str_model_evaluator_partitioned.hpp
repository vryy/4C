/*---------------------------------------------------------------------------*/
/*! \file
\brief model evaluator for structure part of partitioned pasi
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PASI_STR_MODEL_EVALUATOR_PARTITIONED_HPP
#define BACI_PASI_STR_MODEL_EVALUATOR_PARTITIONED_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_structure_new_model_evaluator_generic.hpp"

BACI_NAMESPACE_OPEN

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
      const Teuchos::RCP<Epetra_Vector>& GetInterfaceForceNpPtr()
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
      bool EvaluateForce() override { return true; };

      //! [derived]
      bool EvaluateStiff() override { return true; };

      //! [derived] not needed in partitioned scheme
      bool EvaluateForceStiff() override { return true; };

      //! [derived]
      void PreEvaluate() override { return; };

      //! [derived]
      void PostEvaluate() override { return; };

      //! [derived]
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      bool AssembleJacobian(
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
      void ReadRestart(IO::DiscretizationReader& ioreader) override { return; };

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
      void DetermineStressStrain() override { return; };

      //! [derived]
      void DetermineEnergy() override { return; };

      //! [derived]
      void DetermineOptionalQuantity() override { return; };

      //! [derived]
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override { return; };

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! [derived]
      void ResetStepState() override { return; };

      //! [derived]
      void PostOutput() override { return; };

      //! [derived]
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! @}

     private:
      //! pasi interface force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> interface_force_np_ptr_;
    };

  }  // namespace MODELEVALUATOR

}  // namespace STR

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
