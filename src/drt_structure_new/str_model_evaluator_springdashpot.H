/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all spring dashpot terms

\maintainer Amadeus Gebauer

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_SPRINGDASHPOT_H
#define STR_MODEL_EVALUATOR_SPRINGDASHPOT_H

#include "str_model_evaluator_generic.H"
#include "../drt_constraint/springdashpot.H"
#include <Teuchos_RCP.hpp>

// forward declarations
namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class SpringDashpot : public Generic
    {
     public:
      //! constructor
      SpringDashpot();

      //! destructor
      virtual ~SpringDashpot() = default;

      virtual void Setup() override;

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_springdashpot; }

      //! derived
      virtual void Reset(const Epetra_Vector& x) override;

      //! derived
      virtual bool EvaluateForce() override;

      //! derived
      virtual bool EvaluateStiff() override;

      //! derived
      virtual bool EvaluateForceStiff() override;

      //! derived
      void PreEvaluate() override { return; };

      //! derived
      void PostEvaluate() override { return; };

      //! derived
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew);

      //! derived
      void RunPostIterate(const NOX::Solver::Generic& solver) override { return; };

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void DetermineStressStrain() override;

      //! derived
      void DetermineEnergy() override;

      //! derived
      void DetermineOptionalQuantity() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void ResetStepState() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! [derived]
      void PostOutput() override;

     private:
      //! all spring dashpot instances
      std::vector<Teuchos::RCP<UTILS::SpringDashpot>> springs_;

      //! number of spring dashpot conditions
      int n_conds_;

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! structural velocity at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> velnp_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<LINALG::SparseMatrix> stiff_spring_ptr_;

      //! spring forces at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fspring_np_ptr_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif /* STR_MODEL_EVALUATOR_SPRINGDASHPOT_H */
