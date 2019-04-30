/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_meshtying.H

\brief Evaluation and assembly of all meshtying terms

\maintainer Matthias Mayr

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_

#include "str_model_evaluator_generic.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"

// forward declarations
namespace CONTACT
{
  class Manager;
  class MtAbstractStrategy;
}  // namespace CONTACT

namespace MORTAR
{
  class StrategyBase;
}  // namespace MORTAR

namespace STR
{
  namespace MODELEVALUATOR
  {
    class MeshtyingData;

    /*! \brief Model evaluator for meshtying problems
     *
     */
    class Meshtying : public Generic
    {
     public:
      //! constructor
      Meshtying();

      //! destructor
      virtual ~Meshtying() = default;

      /*! \brief Initialize class variables [derived]
       *
       * @param eval_data_ptr
       * @param gstate_ptr
       * @param gio_ptr
       * @param int_ptr
       * @param timint_ptr
       * @param dof_offset
       */
      virtual void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset);

      //! setup class variables [derived]
      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //!@{

      //! [derived]
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_meshtying; }

      //! [derived]
      void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) override;

      //! [derived]
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! [derived]
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! [derived]
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! [derived]
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! [derived]
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! [derived]
      void RunPostIterate(const NOX::Solver::Generic& solver) override { return; };

      //! [derived]
      void RunPostApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void RunPreApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp) override;

      //! [derived]
      void UpdateStepState(const double& timefac_n) override { return; };

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

      //! [derived]
      void ResetStepState() override { return; };

      //! [derived]
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! [derived]
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! [derived]
      void PostOutput() override { return; };

      /*! \brief Reset model specific variables (without jacobian) [derived]
       *
       * Nothing to do in case of meshtying.
       *
       * \param[in] x Current full state vector
       */
      void Reset(const Epetra_Vector& x) override { return; };

      //! \brief Perform actions just before the Evaluate() call [derived]
      void PreEvaluate() override { return; };

      //! \brief Perform actions right after the Evaluate() call [derived]
      void PostEvaluate() override { return; };

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

      virtual Teuchos::RCP<LINALG::SparseOperator> GetAuxDisplJacobian() const
      {
        return Teuchos::null;
      };

      void EvaluateWeightedGapGradientError();

      //! [derived]
      bool EvaluateForce() override;

      //! [derived]
      bool EvaluateStiff() override;

      //! [derived]
      bool EvaluateForceStiff() override;

      /*! Apply mesh relocation for conservation of angular momentum
       *
       * @param Xslavemod
       */
      void ApplyMeshInitialization(Teuchos::RCP<Epetra_Vector> Xslavemod);

      //!@}

      //! @name Accessors
      //!@{

      //! Returns a pointer to the underlying meshtying strategy object
      const Teuchos::RCP<CONTACT::MtAbstractStrategy>& StrategyPtr();

      //! Returns the underlying meshtying strategy object
      CONTACT::MtAbstractStrategy& Strategy();
      const CONTACT::MtAbstractStrategy& Strategy() const;

      //!@}

     protected:
     private:
      //! meshtying strategy
      Teuchos::RCP<CONTACT::MtAbstractStrategy> strategy_ptr_;

      //! Mesh relocation for conservation of angular momentum
      Teuchos::RCP<Epetra_Vector> mesh_relocation_;
    };  // namespace MODELEVALUATOR

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_ */
