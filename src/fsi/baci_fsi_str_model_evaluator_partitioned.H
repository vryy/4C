/*-----------------------------------------------------------*/
/*! \file

\brief Model evaluator for structure part of partitioned fsi


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_FSI_STR_MODEL_EVALUATOR_PARTITIONED_H
#define BACI_FSI_STR_MODEL_EVALUATOR_PARTITIONED_H

#include "baci_config.H"

#include "baci_structure_new_model_evaluator_generic.H"

BACI_NAMESPACE_OPEN

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
      const Teuchos::RCP<Epetra_Vector>& GetInterfaceForceNpPtr()
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
      bool EvaluateForce() override { return true; };

      //! [derived]
      bool EvaluateStiff() override { return true; };

      //! [derived] not needed in partitioned scheme
      bool EvaluateForceStiff() override { return true; };

      //! derived
      void PreEvaluate() override { return; };

      //! derived
      void PostEvaluate() override { return; };

      //! derived
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
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
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

      //! derived
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

      //! derived
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

      /*! \brief linear structure solve with just a interface load
       *
       * The very special solve done in steepest descent relaxation
       * calculation (and matrix free Newton Krylov).
       *
       * \note Can only be called after a valid structural solve. */
      Teuchos::RCP<const Epetra_Vector> SolveRelaxationLinear(
          Teuchos::RCP<ADAPTER::Structure> structure);

      virtual void SetupMultiMapExtractor();

      //! set flag true if in order to request a relaxation solve
      void SetIsRelaxationSolve(bool trueorfalse) { is_relaxationsolve = trueorfalse; };

      //! Returns the global input/output data container
      const STR::TIMINT::BaseDataIO& GetInOutput() const;

     private:
      //! fsi interface force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> interface_force_np_ptr_;

      //! true if relaxation solve is requested
      bool is_relaxationsolve;

    };  // class PartitionedFSI

  }  // namespace MODELEVALUATOR
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // FSI_STR_MODEL_EVALUATOR_PARTITIONED_H
