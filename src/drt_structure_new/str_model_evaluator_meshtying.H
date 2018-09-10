/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_meshtying.H

\brief Evaluation and assembly of all meshtying terms

\maintainer Alexander Seitz

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_

#include "str_model_evaluator_generic.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"

// forward declarations...
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

    class Meshtying : public Generic
    {
     public:
      //! constructor
      Meshtying();

      //! destructor
      virtual ~Meshtying(){};

      //! init class variables [derived]
      virtual void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset);

      //! setup class variables [derived]
      virtual void Setup();

      //! @name Functions which are derived from the base generic class
      //! @{
      //! [derived]
      INPAR::STR::ModelType Type() const { return INPAR::STR::model_meshtying; }

      //    //! [derived]
      //    virtual bool EvaluateForce();
      //
      //    //! [derived]
      //    virtual bool EvaluateStiff();
      //
      //    //! [derived]
      //    virtual bool EvaluateForceStiff();
      //
      //    //! [derived]
      //    void PreEvaluate() override;
      //
      //    //! [derived]
      //    virtual void PostEvaluate();

      //! [derived]
      void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) override;

      //! [derived]
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! [derived]
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew);

      //! [derived]
      virtual void RunPreComputeX(
          const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp){};

      //! [derived]
      virtual void RunPostIterate(const NOX::Solver::Generic& solver){};

      //! [derived]
      virtual void RunPostApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp);

      //! [derived]
      virtual void RunPreApplyJacobianInverse(const Epetra_Vector& rhs, Epetra_Vector& result,
          const Epetra_Vector& xold, const NOX::NLN::Group& grp);

      //! [derived]
      virtual void UpdateStepState(const double& timefac_n){};

      //! [derived]
      virtual void UpdateStepElement(){};

      //! [derived]
      virtual void DetermineStressStrain(){};

      //! [derived]
      virtual void DetermineEnergy(){};

      //! [derived]
      virtual void DetermineOptionalQuantity() override{};

      //! [derived]
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const {};

      //! [derived]
      virtual void ResetStepState(){};

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

      //! [derived]
      virtual void PostOutput(){};

      //! @}

      /// @name call-back routines
      /// @{

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

      inline const DRT::DiscretizationInterface& StrDiscret() const { return Discret(); }

      //! [derived]
      virtual bool EvaluateForce();

      //! [derived]
      virtual bool EvaluateStiff();

      //! [derived]
      virtual bool EvaluateForceStiff();

      //! mesh relocation for conservation of angular momentum
      void ApplyMeshInitialization(Teuchos::RCP<Epetra_Vector> Xslavemod);

      /// @}

      //! @name Accessors
      //! @{

      //! Returns a pointer to the underlying meshtying strategy object
      const Teuchos::RCP<CONTACT::MtAbstractStrategy>& StrategyPtr();

      //! Returns the underlying meshtying strategy object
      CONTACT::MtAbstractStrategy& Strategy();
      const CONTACT::MtAbstractStrategy& Strategy() const;

      //! @}

      // things we don't need
      virtual void Reset(const Epetra_Vector& x) { return; };
      virtual void PostEvaluate(){};

     protected:
      STR::MODELEVALUATOR::MeshtyingData& EvalMeshtying();
      const STR::MODELEVALUATOR::MeshtyingData& EvalMeshtying() const;

     private:
      //! meshtying strategy
      Teuchos::RCP<CONTACT::MtAbstractStrategy> strategy_ptr_;
      //! mesh relocation of conservation of angular momentum
      Teuchos::RCP<Epetra_Vector> mesh_relocation_;
    };  // class Meshtying

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_MESHTYING_H_ */
