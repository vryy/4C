/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all contact terms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_MODEL_EVALUATOR_CONTACT_H
#define BACI_STRUCTURE_NEW_MODEL_EVALUATOR_CONTACT_H

#include "baci_lib_utils_matrix_vector_enums.H"
#include "baci_structure_new_model_evaluator_generic.H"

// forward declarations
namespace CONTACT
{
  class Manager;
  class CoAbstractStrategy;
}  // namespace CONTACT

namespace MORTAR
{
  class StrategyBase;
}  // namespace MORTAR

namespace STR
{
  namespace MODELEVALUATOR
  {
    class ContactData;

    class Contact : public Generic
    {
     public:
      //! setup class variables [derived]
      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //!@{

      //! [derived]
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_contact; }

      //! reset class variables (without jacobian) [derived]
      void Reset(const Epetra_Vector& x) override;

      //! [derived]
      bool EvaluateForce() override;

      //! [derived]
      bool EvaluateStiff() override;

      //! [derived]
      bool EvaluateForceStiff() override;

      //! [derived]
      void PreEvaluate() override;

      //! [derived]
      void PostEvaluate() override;

      //! [derived]
      void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) override;

      //! [derived]
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool AssembleJacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! Perform a correction of adaptive parameters
      bool CorrectParameters(NOX::NLN::CorrectionType type) override;

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! [derived]
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! recover condensed Lagrange multipliers
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! [derived]
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;

      //! [derived]
      void RunPostIterate(const NOX::Solver::Generic& solver) override;

      /// [derived]
      void RunPreSolve(const NOX::Solver::Generic& solver) override;

      //! [derived]
      void RunPostApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void RunPreApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void UpdateStepState(const double& timefac_n) override;

      //! [derived]
      void UpdateStepElement() override;

      //! [derived]
      void DetermineStressStrain() override;

      //! [derived]
      void DetermineEnergy() override;

      //! [derived]
      void DetermineOptionalQuantity() override;

      //! [derived]
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! [derived]
      void ResetStepState() override;

      //! [derived]
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! [derived]
      void PostOutput() override;

      //! [derived]
      bool EvaluateCheapSOCRhs() override;

      //! [derived]
      bool AssembleCheapSOCRhs(Epetra_Vector& f, const double& timefac_np) const override;

      //! @}

      //! @name Call-back routines
      //!@{

      Teuchos::RCP<const CORE::LINALG::SparseMatrix> GetJacobianBlock(
          const DRT::UTILS::MatBlockType bt) const;

      /** \brief Assemble the structural right-hand side vector
       *
       *  \param[in] without_these_models  Exclude all models defined in this vector
       *                                   during the assembly
       *  \param[in] apply_dbc             Apply Dirichlet boundary conditions
       *
       *  \author hiermeier \date 08/17 */
      Teuchos::RCP<Epetra_Vector> AssembleForceOfModels(
          const std::vector<INPAR::STR::ModelType>* without_these_models = nullptr,
          const bool apply_dbc = false) const;

      virtual Teuchos::RCP<CORE::LINALG::SparseOperator> GetAuxDisplJacobian() const;

      void EvaluateWeightedGapGradientError();

      //!@}

      //! @name Accessors
      //!@{

      //! Returns a pointer to the underlying contact strategy object
      const Teuchos::RCP<CONTACT::CoAbstractStrategy>& StrategyPtr();

      //! Returns the underlying contact strategy object
      CONTACT::CoAbstractStrategy& Strategy();
      const CONTACT::CoAbstractStrategy& Strategy() const;

      //!@}

     protected:
      STR::MODELEVALUATOR::ContactData& EvalContact();
      const STR::MODELEVALUATOR::ContactData& EvalContact() const;

      virtual void CheckPseudo2D() const;

     private:
      void PostSetup(Teuchos::ParameterList& cparams);

      /// Set the correct time integration parameters within the contact strategy
      void SetTimeIntegrationInfo(CONTACT::CoAbstractStrategy& strategy) const;

      void PostUpdateStepState();

      void ExtendLagrangeMultiplierDomain(Teuchos::RCP<Epetra_Vector>& lm_vec) const;

      //! contact evaluation data container
      Teuchos::RCP<STR::MODELEVALUATOR::ContactData> eval_contact_ptr_;

      //! contact strategy
      Teuchos::RCP<CONTACT::CoAbstractStrategy> strategy_ptr_;

    };  // class Contact

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif  // STRUCTURE_NEW_MODEL_EVALUATOR_CONTACT_H
