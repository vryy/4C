/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_spherebeamlinking.H

\brief class for managing rigid sphere to beam crosslinking

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_H_

#include "../drt_inpar/inpar_beaminteraction.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include <Epetra_MpiComm.h>

#include "../drt_binstrategy/binning_strategy_utils.H"
#include "beaminteraction_submodel_evaluator_generic.H"
#include <unordered_set>

// forward declaration
class RuntimeVtpWriter;

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

      //! destructor
      virtual ~SphereBeamLinking(){};

      //! setup class variables
      virtual void Setup();

      //! derived
      virtual void PostSetup();

      //! Returns the type of the current model evaluator
      INPAR::BEAMINTERACTION::SubModelType Type() const
      {
        return INPAR::BEAMINTERACTION::submodel_spherebeamlink;
      }

      //! derived
      virtual void Reset();

      //! derived
      virtual bool EvaluateForce();

      //! derived
      virtual bool EvaluateStiff();

      //! derived
      virtual bool EvaluateForceStiff();

      //! derived
      virtual void UpdateStepState(const double& timefac_n);

      //! derived
      virtual bool PreUpdateStepElement(bool beam_redist);

      //! derived
      virtual void UpdateStepElement(bool repartition_was_done);

      //! derived
      virtual void PostUpdateStepElement();

      //! derived
      std::map<STR::EnergyType, double> GetEnergy() const override;

      //! derived
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

      //! derived
      virtual void RuntimeOutputStepState() const;

      //! derived
      virtual void ResetStepState();

      //! derived
      virtual void WriteRestart(
          IO::DiscretizationWriter& ia_writer, IO::DiscretizationWriter& bin_writer) const;

      //! derived
      virtual void PreReadRestart() override;

      //! derived
      virtual void ReadRestart(
          IO::DiscretizationReader& ia_reader, IO::DiscretizationReader& bin_reader);

      //! derived
      virtual void PostReadRestart();

      //! derived
      virtual void RunPostIterate(const NOX::Solver::Generic& solver){/*empty*/};

      //! derived
      virtual void InitSubmodelDependencies(
          Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap);

      //! derived
      virtual void AddBinsToBinColMap(std::set<int>& colbins);

      //! derived
      virtual void AddBinsWithRelevantContentForIaDiscretColMap(std::set<int>& colbins) const;

      //! derived
      virtual void GetHalfInteractionDistance(double& half_interaction_distance);

      //! @}

     private:
      //! @name submodel specific, not derived member functions
      //! @{

      //! init output in VTP format
      void InitOutputRuntimeVtp();

      //! writes output for discretization structure in VTP format
      void WriteOutputRuntimeVtp() const;

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
      virtual void UnbindSphereBeamBonds(int& num_disolved);

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
      Teuchos::RCP<RuntimeVtpWriter> vtp_writer_ptr_;

      //! step number for random stuff concerning sphere beam linking
      int random_number_sphere_beam_linking_step_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

#endif
