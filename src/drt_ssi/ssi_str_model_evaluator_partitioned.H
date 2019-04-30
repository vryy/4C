/*-----------------------------------------------------------*/
/*!
\file ssi_str_model_evaluator_partitioned.H

\brief Model evaluator for structure part of partitioned ssi.

\maintainer Christoph Schmidt

\level 3
*/
/*-----------------------------------------------------------*/


#ifndef SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_
#define SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"

namespace ADAPTER
{
  class Structure;
}  // namespace ADAPTER

namespace SSI
{
  class SSI_Part;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedSSI : public Generic
    {
     public:
      //! constructor
      PartitionedSSI(const Teuchos::RCP<const SSI::SSI_Part>
              ssi_part  //!< partitioned algorithm for scalar-structure interaction
      );

      //! destructor
      virtual ~PartitionedSSI() = default;

      //! setup class variables [derived]
      virtual void Setup();

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
      void PostEvaluate() { return; };

      //! derived
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override
      {
        return;
      };

      //! [derived]
      void ReadRestart(IO::DiscretizationReader& ioreader) override { return; };

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! derived
      virtual void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;

      //! recover condensed Lagrange multipliers
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
        return;
      };

      //! derived
      void RunPostIterate(const NOX::Solver::Generic& solver) override { return; };

      //! [derived]
      void UpdateStepState(const double& timefac_n) override;

      //! [derived]
      void UpdateStepElement() override { return; };

      //! [derived]
      virtual void DetermineStressStrain() { return; };

      //! [derived]
      virtual void DetermineEnergy() { return; };

      //! [derived]
      virtual void DetermineOptionalQuantity() override { return; };

      //! [derived]
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const { return; };

      //! derived
      virtual void ResetStepState() { return; };

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
      //! partitioned algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSI_Part> ssi_part_;
    };  // class PartitionedSSI

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_ */
