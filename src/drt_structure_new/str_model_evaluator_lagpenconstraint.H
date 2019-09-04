/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all constraint terms

\maintainer Amadeus Gebauer

\level 3
*/
/*---------------------------------------------------------------------*/

#ifndef STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H
#define STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

#include "../drt_constraint/constraint_manager.H"

// forward declarations
namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace LAGPENCONSTRAINT
{
  class NoxInterface;
  class NoxInterfacePrec;
}  // namespace LAGPENCONSTRAINT

namespace STR
{
  namespace MODELEVALUATOR
  {
    class LagPenConstraint : public Generic
    {
     public:
      //! constructor
      LagPenConstraint();

      //! destructor
      virtual ~LagPenConstraint() = default;

      void Setup() override;

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_lag_pen_constraint; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override;

      //! derived
      bool EvaluateForce() override;

      //! derived
      bool EvaluateStiff() override;

      //! derived
      bool EvaluateForceStiff() override;

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
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

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

      //! [derived]
      void PostOutput() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      const Teuchos::RCP<UTILS::ConstrManager>& StrategyPtr();

      //! Return the NOX::NLN::CONSTRAINT::Interface::Required member object
      const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface>& NoxInterfacePtr();

      //! Return the NOX::NLN::CONSTRAINT::Interface::Preconditioner member object
      const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec>& NoxInterfacePrecPtr();

     protected:
      //! Returns the underlying contact strategy object
      UTILS::ConstrManager& Strategy();
      const UTILS::ConstrManager& Strategy() const;

     private:
      //! all constraint instances
      Teuchos::RCP<UTILS::ConstrManager> constrman_;  //!< Constraint manager

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<LINALG::SparseMatrix> stiff_constr_ptr_;

      //! constraint contributions to the structural rhs at \f%t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fstrconstr_np_ptr_;

      //! pointer to the NOX::NLN::CONSTRAINT::Interface::Required object
      Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface> noxinterface_ptr_;

      //! pointer to the NOX::NLN::CONSTRAINT::Interface::Preconditioner object
      Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec> noxinterface_prec_ptr_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H */
