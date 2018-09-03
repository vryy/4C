/*-----------------------------------------------------------*/
/*!
\file ssi_str_model_evaluator_partitioned.H

\brief Model evaluator for structure part of partitioned ssi.

\maintainer Andreas Rauch

\date Nov 11, 2016

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
      virtual ~PartitionedSSI(){};

      //! setup class variables [derived]
      virtual void Setup();

      //! @name Functions which are derived from the base generic class
      //! @{
      //! [derived]
      INPAR::STR::ModelType Type() const { return INPAR::STR::model_partitioned_coupling; }

      //! reset class variables (without jacobian) [derived]
      virtual void Reset(const Epetra_Vector& x) { return; };

      //! [derived]
      virtual bool EvaluateForce() { return true; };

      //! [derived]
      virtual bool EvaluateStiff() { return true; };

      //! [derived] not needed in partitioned scheme
      virtual bool EvaluateForceStiff() { return true; };

      //! derived
      virtual void PostEvaluate(){/* currently unused */};

      //! derived
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
      {
        return;
      };

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader) { return; };

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! derived
      virtual void RunPreComputeX(
          const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp);

      //! recover condensed Lagrange multipliers
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
      {
        return;
      };

      //! derived
      virtual void RunPostIterate(const NOX::Solver::Generic& solver) { return; };

      //! [derived]
      virtual void UpdateStepState(const double& timefac_n);

      //! [derived]
      virtual void UpdateStepElement() { return; };

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
      virtual void PostOutput() { return; };

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
      //! @}

     private:
      //! partitioned algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSI_Part> ssi_part_;
    };  // class PartitionedSSI

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_ */
