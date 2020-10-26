/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for monolithic scalar-structure interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_
#define SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"

// forward declarations
class Epetra_MultiVector;

namespace SSI
{
  class SSIMono;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class MonolithicSSI : public Generic
    {
     public:
      //! constructor
      MonolithicSSI(const Teuchos::RCP<const SSI::SSIMono>
              ssi_mono  //!< monolithic algorithm for scalar-structure interaction
      );

      //! destructor
      ~MonolithicSSI() override = default;

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
      void DetermineEnergy() override{};

      // calculate stresses and strains
      void DetermineStressStrain() override;

      // calculate an optional quantity
      void DetermineOptionalQuantity() override{};

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
      void PreEvaluate() override{};

      //! post-evaluate
      void PostEvaluate() override{};

      //! post-output
      void PostOutput() override{};

      //! predict solution
      void Predict(const INPAR::STR::PredEnum& pred_type) override{};

      //! read model-specific restart
      void ReadRestart(IO::DiscretizationReader& ioreader) override{};

      //! recover condensed solution variables
      void RunPostComputeX(const Epetra_Vector& xold, const Epetra_Vector& dir,
          const Epetra_Vector& xnew) override{};

      //! reset model-specific variables (without Jacobian)
      void Reset(const Epetra_Vector& x) override{};

      //! reset state
      void ResetStepState() override{};

      //! post-iterate
      void RunPostIterate(const NOX::Solver::Generic& solver) override{};

      //! pre-compute
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override{};

      //! set up model evaluator
      void Setup() override;

      //! return nodal stresses
      const Epetra_MultiVector& Stresses() const { return *stresses_; };

      //! return type of model evaluator
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_monolithic_coupling; };

      //! update element
      void UpdateStepElement() override{};

      //! update state
      void UpdateStepState(const double& timefac_n) override{};

      //! write model-specific restart
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

     private:
      //! nodal stresses
      Teuchos::RCP<Epetra_MultiVector> stresses_;

      //! monolithic algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSIMono> ssi_mono_;
    };  // class MonolithicSSI
  }     // namespace MODELEVALUATOR
}  // namespace STR
#endif
