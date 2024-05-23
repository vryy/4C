/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of all beam interaction terms


\level 3
*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BEAMCONTACT_STR_MODEL_EVALUATOR_BEAMINTERACTION_OLD_HPP
#define FOUR_C_BEAMCONTACT_STR_MODEL_EVALUATOR_BEAMINTERACTION_OLD_HPP

#include "4C_config.hpp"

#include "4C_structure_new_model_evaluator_generic.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  class Beam3cmanager;
}  // namespace CONTACT

namespace STR
{
  namespace MODELEVALUATOR
  {
    class BeamInteractionOld : public Generic
    {
     public:
      //! constructor
      BeamInteractionOld();

      void Setup() override;

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_beam_interaction_old; }

      //! derived
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
      void PostEvaluate() override{/* currently unused */};

      //! derived
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      bool AssembleJacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

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
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void DetermineEnergy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void ResetStepState() override;

      //! [derived]
      void PostOutput() override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

     private:
      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! stiffness contributions from beam interactions
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_beaminteract_ptr_;

      //! force contributions from beam interaction at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> f_beaminteract_np_ptr_;

      //! beam contact manager
      Teuchos::RCP<CONTACT::Beam3cmanager> beamcman_;
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
