/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_cardiovascular0d.H

\brief Evaluation and assembly of all 0D cardiovascular model terms

\maintainer Amadeus Gebauer

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef STR_MODEL_EVALUATOR_CARDIOVASCULAR0D_H
#define STR_MODEL_EVALUATOR_CARDIOVASCULAR0D_H

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

#include "../drt_cardiovascular0d/cardiovascular0d_manager.H"

// forward declaration ...
namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Cardiovascular0D : public Generic
    {
     public:
      //! constructor
      Cardiovascular0D();

      //! destructor
      virtual ~Cardiovascular0D(){};

      virtual void Setup();

      //! derived
      INPAR::STR::ModelType Type() const { return INPAR::STR::model_cardiovascular0d; }

      //! derived
      virtual void Reset(const Epetra_Vector& x);

      //! derived
      virtual bool EvaluateForce();

      //! derived
      virtual bool EvaluateStiff();

      //! derived
      virtual bool EvaluateForceStiff();

      //! derived
      virtual void PostEvaluate(){/* currently unused */};

      //! derived
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! derived
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! derived
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! derived
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! derived
      virtual void RunPreComputeX(
          const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
      {
        return;
      };

      //! derived
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew);

      //! derived
      virtual void RunPostIterate(const NOX::Solver::Generic& solver) { return; };

      //! derived
      virtual void UpdateStepState(const double& timefac_n);

      //! derived
      virtual void UpdateStepElement();

      //! derived
      virtual void DetermineStressStrain();

      //! derived
      virtual void DetermineEnergy();

      //! derived
      virtual void DetermineOptionalQuantity() override;

      //! derived
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

      //! derived
      virtual void ResetStepState();

      //! [derived]
      virtual void PostOutput();

      //! derived
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

      //! derived
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

      //! derived
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

     private:
      Teuchos::RCP<UTILS::Cardiovascular0DManager> cardvasc0dman_;  //!< Cardiovascular0D manager

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! cardio contributions to the structural stiffness matrix
      Teuchos::RCP<LINALG::SparseMatrix> stiff_cardio_ptr_;

      //! structural rhs contributions of the cardio model at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fstructcardio_np_ptr_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* STR_MODEL_EVALUATOR_CARDIOVASCULAR0D_H */
