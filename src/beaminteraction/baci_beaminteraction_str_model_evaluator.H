/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of all beam interaction terms


\level 3
*/
/*-----------------------------------------------------------*/


#ifndef BACI_BEAMINTERACTION_STR_MODEL_EVALUATOR_H
#define BACI_BEAMINTERACTION_STR_MODEL_EVALUATOR_H

#include "baci_coupling_adapter.H"
#include "baci_inpar_beaminteraction.H"
#include "baci_linalg_mapextractor.H"
#include "baci_structure_new_enum_lists.H"
#include "baci_structure_new_model_evaluator_generic.H"  // base class


// forward declaration ...
namespace ADAPTER
{
  class Coupling;
}

namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class SparseMatrix;
  class MatrixRowTransform;
}  // namespace CORE::LINALG

namespace BINSTRATEGY
{
  class BinningStrategy;
}

namespace BEAMINTERACTION
{
  class BeamInteractionParams;

  class BeamCrosslinkerHandler;

  namespace SUBMODELEVALUATOR
  {
    class Generic;
  }
}  // namespace BEAMINTERACTION

namespace STR
{
  namespace MODELEVALUATOR
  {
    // forward declaration
    class BeamInteractionDataState;

    class BeamInteraction : public Generic
    {
     public:
      typedef std::map<enum INPAR::BEAMINTERACTION::SubModelType,
          Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic>>
          Map;
      typedef std::vector<Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic>> Vector;

      //! constructor
      BeamInteraction();

      void Setup() override;

      virtual void PostSetup();

      /// print welcome to biopolymer network simulation
      virtual void Logo() const;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{
      //! derived

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_beaminteraction; }

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

      //! derived
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
      void RunPostIterate(const NOX::Solver::Generic& solver) override;

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
      void RuntimeOutputStepState() const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! derived
      void PostOutput() override;

      //! derived
      void ResetStepState() override;
      //! @}

     protected:
      //! derived
      void Reset(const Epetra_Vector& x) override;

      //!@name routines for submodel management
      //! @{

     public:
      /// check if the given model type is active.
      bool HaveSubModelType(INPAR::BEAMINTERACTION::SubModelType const& submodeltype) const;

     private:
      void PartitionProblem();

      bool PostPartitionProblem();

      //! set beaminteraction sub models
      void SetSubModelTypes();


      //! build, init and setup submodel evaluator
      void InitAndSetupSubModelEvaluators();

      //! give submodels a certain order in which they are evaluated
      virtual Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> TransformToVector(
          STR::MODELEVALUATOR::BeamInteraction::Map submodel_map,
          std::vector<INPAR::BEAMINTERACTION::SubModelType>& sorted_submodel_types) const;

      //! @}

      //!@name routines that manage to discretizations with distinct parallel distribution
      //! @{

      /// check if interaction discretization needs to be redistributed completely
      bool CheckIfBeamDiscretRedistributionNeedsToBeDone();

      /// update coupling adapter and matrix transformation object with new maps
      void UpdateCouplingAdapterAndMatrixTransformation();

      /// transform force vector from ia_discret_ to Discret()
      virtual void TransformForce();

      /// transform stiffness matrix from ia_discret_ to Discret()
      virtual void TransformStiff();

      /// transform force vector and stiffness matrix from ia_discret_ to Discret()
      virtual void TransformForceStiff();

      /// update states based on bindis after its redistribution
      virtual void UpdateMaps();

      //! @}

      //!@name routines that manage binning strategy
      //! @{

      /// change parallel distribution of bindis and ia_discret and assign (beam) eles to bins
      virtual void ExtendGhosting();

      /// build ele to bin map
      virtual void BuildRowEleToBinMap();

      /// print some information about binning
      virtual void PrintBinningInfoToScreen() const;

      //! @}

     private:
      //! pointer to the problem discretization (cast of base class member)
      Teuchos::RCP<::DRT::Discretization> discret_ptr_;

      //! data container holding all beaminteraction related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamInteractionParams> beaminteraction_params_ptr_;

      //!@name data for submodel management
      //! @{
      /// current active model types for the model evaluator
      Teuchos::RCP<std::set<enum INPAR::BEAMINTERACTION::SubModelType>> submodeltypes_;

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> me_map_ptr_;

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> me_vec_ptr_;
      //! @}

      //!@name data for handling two distinct parallel distributed discretizations
      //! @{
      //! myrank
      int myrank_;

      //! coupling adapter to transfer vectors and matrices between Discret() and intactids_
      Teuchos::RCP<CORE::ADAPTER::Coupling> coupsia_;

      //! transform object for structure stiffness matrix
      Teuchos::RCP<CORE::LINALG::MatrixRowTransform> siatransform_;
      //! @}


      //!@name data for beaminteraction with binning strategy
      //! @{
      //! interaction discretization handling all interactions (e.g. crosslinker to beam,
      //! beam to beam, potential ...)
      Teuchos::RCP<::DRT::Discretization> ia_discret_;

      /// map extractor for split of different element types
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor> eletypeextractor_;

      //! pointer to the global state data container
      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> ia_state_ptr_;

      //! force based on ia_discret at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> ia_force_beaminteraction_;

      //! global force based on Discret() at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> force_beaminteraction_;

      //! structural stiffness matrix based on Discret()
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_beaminteraction_;

      //! beam crosslinker handler
      Teuchos::RCP<BEAMINTERACTION::BeamCrosslinkerHandler> beam_crosslinker_handler_;

      //! binning strategy
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

      //! crosslinker and bin discretization
      Teuchos::RCP<::DRT::Discretization> bindis_;

      //! elerowmap of bindis
      Teuchos::RCP<Epetra_Map> rowbins_;

      //! displacement of nodes since last redistribution
      Teuchos::RCP<Epetra_Vector> dis_at_last_redistr_;

      //! depending on which submodels are active this variable has different values
      //! and determines how often a redistribution needs to be done
      double half_interaction_distance_;

      //! @}
    };
  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif
