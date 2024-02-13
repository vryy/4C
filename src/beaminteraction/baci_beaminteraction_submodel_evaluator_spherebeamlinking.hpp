/*-----------------------------------------------------------*/
/*! \file

\brief class for managing rigid sphere to beam crosslinking


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_HPP
#define BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_HPP

#include "baci_config.hpp"

#include "baci_beaminteraction_submodel_evaluator_generic.hpp"
#include "baci_binstrategy_utils.hpp"
#include "baci_inpar_beaminteraction.hpp"
#include "baci_io_visualization_manager.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

#include <Epetra_MpiComm.h>

#include <unordered_set>

BACI_NAMESPACE_OPEN


namespace DRT
{
  class Element;
  class Node;
}  // namespace DRT
namespace BEAMINTERACTION
{
  class SphereBeamLinkingParams;
  class BeamLinkPinJointed;

  namespace SUBMODELEVALUATOR
  {
    class SphereBeamLinking : public Generic
    {
     public:
      //! constructor
      SphereBeamLinking();

      //! setup class variables
      void Setup() override;

      //! derived
      void PostSetup() override;

      //! Returns the type of the current model evaluator
      INPAR::BEAMINTERACTION::SubModelType Type() const override
      {
        return INPAR::BEAMINTERACTION::submodel_spherebeamlink;
      }

      //! derived
      void Reset() override;

      //! derived
      bool EvaluateForce() override;

      //! derived
      bool EvaluateStiff() override;

      //! derived
      bool EvaluateForceStiff() override;

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      bool PreUpdateStepElement(bool beam_redist) override;

      //! derived
      void UpdateStepElement(bool repartition_was_done) override;

      //! derived
      void PostUpdateStepElement() override;

      //! derived
      std::map<STR::EnergyType, double> GetEnergy() const override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void RuntimeOutputStepState() const override;

      //! derived
      void ResetStepState() override;

      //! derived
      void WriteRestart(
          IO::DiscretizationWriter& ia_writer, IO::DiscretizationWriter& bin_writer) const override;

      //! derived
      void PreReadRestart() override;

      //! derived
      void ReadRestart(
          IO::DiscretizationReader& ia_reader, IO::DiscretizationReader& bin_reader) override;

      //! derived
      void PostReadRestart() override;

      //! derived
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override{/*empty*/};

      //! derived
      void InitSubmodelDependencies(
          Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap) override;

      //! derived
      void AddBinsToBinColMap(std::set<int>& colbins) override;

      //! derived
      void AddBinsWithRelevantContentForIaDiscretColMap(std::set<int>& colbins) const override;

      //! derived
      void GetHalfInteractionDistance(double& half_interaction_distance) override;

      //! @}

     private:
      //! @name submodel specific, not derived member functions
      //! @{

      //! init output in VTP format
      void InitOutputRuntime();

      //! writes output for discretization structure in VTP format
      void WriteOutputRuntime() const;

      /// get neighboring eles in discret
      virtual void FindAndStoreNeighboringElements(
          std::map<int, std::vector<std::pair<int, int>>>& newlinks);

      /// get neighboring eles in discret
      virtual void CheckFeasibilityOfNewLink(DRT::Element const* currele,
          std::vector<DRT::Element const*> const& neighbors, std::unordered_set<int>& tobebonded,
          std::map<int, std::vector<std::pair<int, int>>>& newlinks) const;

      /// create new beam to sphere joint object
      virtual void CreateBeamToSphereJoint(
          std::map<int, std::vector<std::pair<int, int>>> const& newlinks);

      /// update linker length to mimic contractive cells
      virtual void UpdateLinkerLength();

      /// check if bond needs to be dissolved
      virtual void UnbindSphereBeamBonds(int& num_dissolved);

      /// compute force dependent off rate for a catch-slip bond
      virtual void CalcForceDependentCatchSlipBondUnbindProbability(
          Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr, double& p_unbind);

      //! @}

     private:
      //! cell discretization
      Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking> sm_crosslinkink_ptr;

      //! data container holding all beam contact related parameters
      Teuchos::RCP<BEAMINTERACTION::SphereBeamLinkingParams> spherebeamlinking_params_ptr_;

      //! runtime output for cell beam crosslinks (integrins)
      Teuchos::RCP<IO::VisualizationManager> visualization_manager_ptr_;

      //! step number for random stuff concerning sphere beam linking
      int random_number_sphere_beam_linking_step_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
