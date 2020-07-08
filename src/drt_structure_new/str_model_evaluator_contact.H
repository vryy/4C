/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all contact terms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_

#include "str_model_evaluator_generic.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"

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
      //! constructor
      Contact();

      //! destructor
      virtual ~Contact() = default;

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
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! Perform a correction of adaptive parameters
      virtual bool CorrectParameters(NOX::NLN::CorrectionType type);

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! recover condensed Lagrange multipliers
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew);

      //! [derived]
      virtual void RunPreComputeX(
          const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp);

      //! [derived]
      virtual void RunPostIterate(const NOX::Solver::Generic& solver);

      /// [derived]
      void RunPreSolve(const NOX::Solver::Generic& solver) override;

      //! [derived]
      virtual void RunPostApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp);

      //! [derived]
      virtual void RunPreApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp);

      //! [derived]
      virtual void UpdateStepState(const double& timefac_n);

      //! [derived]
      virtual void UpdateStepElement();

      //! [derived]
      virtual void DetermineStressStrain();

      //! [derived]
      virtual void DetermineEnergy();

      //! [derived]
      virtual void DetermineOptionalQuantity() override;

      //! [derived]
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

      //! [derived]
      virtual void ResetStepState();

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

      //! [derived]
      virtual void PostOutput();

      //! [derived]
      virtual bool EvaluateCheapSOCRhs();

      //! [derived]
      virtual bool AssembleCheapSOCRhs(Epetra_Vector& f, const double& timefac_np) const;

      //! @}

      //! @name Call-back routines
      //!@{

      Teuchos::RCP<const LINALG::SparseMatrix> GetJacobianBlock(
          const DRT::UTILS::MatBlockType bt) const;

      /** \brief Assemble the structural right-hand side vector
       *
       *  \param[in] without_these_models  Exclude all models defined in this vector
       *                                   during the assembly
       *  \param[in] apply_dbc             Apply Dirichlet boundary conditions
       *
       *  \author hiermeier \date 08/17 */
      Teuchos::RCP<Epetra_Vector> AssembleForceOfModels(
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL,
          const bool apply_dbc = false) const;

      virtual Teuchos::RCP<LINALG::SparseOperator> GetAuxDisplJacobian() const;

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

#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_ */
