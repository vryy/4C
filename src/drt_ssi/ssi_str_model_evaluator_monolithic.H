/*----------------------------------------------------------------------*/
/*!
\file ssi_str_model_evaluator_monolithic.H

\brief structural model evaluator for monolithic scalar-structure interaction

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#ifndef SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_
#define SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"

// forward declarations
class Epetra_MultiVector;

namespace SSI
{
  class SSI_Mono;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class MonolithicSSI : public Generic
    {
     public:
      //! constructor
      MonolithicSSI(const Teuchos::RCP<const SSI::SSI_Mono>
              ssi_mono  //!< monolithic algorithm for scalar-structure interaction
      );

      //! destructor
      virtual ~MonolithicSSI() = default;

      //! assemble right-hand side vector
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override
      {
        return true;
      };

      //! assemble Jacobian
      bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const override
      {
        return true;
      };

      //! calculate energy
      void DetermineEnergy() override { return; };

      // calculate stresses and strains
      void DetermineStressStrain() override;

      // calculate an optional quantity
      void DetermineOptionalQuantity() override { return; };

      //! evaluate current right-hand side vector
      bool EvaluateForce() override { return true; };

      //! evaluate current right-hand-side vector and tangential stiffness matrix
      bool EvaluateForceStiff() override { return true; };

      //! evaluate current tangential stiffness matrix
      bool EvaluateStiff() override { return true; };

      //! return dofrowmap
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! return solution
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override
      {
        dserror("Not implemented!");
        return Teuchos::null;
      };

      //! return previous solution
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override
      {
        dserror("Not implemented!");
        return Teuchos::null;
      };

      //! output state
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! pre-evaluate
      void PreEvaluate() override { return; };

      //! post-evaluate
      void PostEvaluate() override { return; };

      //! post-output
      void PostOutput() override { return; };

      //! predict solution
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! read model-specific restart
      void ReadRestart(IO::DiscretizationReader& ioreader) override { return; };

      //! recover condensed solution variables
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override
      {
        return;
      };

      //! reset model-specific variables (without Jacobian)
      void Reset(const Epetra_Vector& x) override { return; };

      //! reset state
      void ResetStepState() override { return; };

      //! post-iterate
      void RunPostIterate(const NOX::Solver::Generic& solver) override { return; };

      //! pre-compute
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! set up model evaluator
      void Setup() override;

      //! return nodal stresses
      const Epetra_MultiVector& Stresses() const { return *stresses_; };

      //! return type of model evaluator
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_monolithic_coupling; };

      //! update element
      void UpdateStepElement() override { return; };

      //! update state
      void UpdateStepState(const double& timefac_n) override { return; };

      //! write model-specific restart
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

     private:
      //! nodal stresses
      Teuchos::RCP<Epetra_MultiVector> stresses_;

      //! monolithic algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono_;
    };  // class MonolithicSSI
  }     // namespace MODELEVALUATOR
}  // namespace STR
#endif
