/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief submodel for potential-based beam interactions


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_HPP
#define BACI_BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_HPP


#include "baci_config.hpp"

#include "baci_beaminteraction_submodel_evaluator_generic.hpp"
#include "baci_binstrategy_utils.hpp"
#include "baci_io_visualization_manager.hpp"

BACI_NAMESPACE_OPEN

// forward declaration ...
namespace DRT
{
  class Element;
  class Condition;
}  // namespace DRT
namespace BEAMINTERACTION
{
  class BeamPotentialParams;
  class BeamPotentialPair;

  namespace SUBMODELEVALUATOR
  {
    class BeamPotential : public Generic
    {
     public:
      //! constructor
      BeamPotential();

      //! setup class variables
      void Setup() override;

      //! derived
      void PostSetup() override;

      //! Returns the type of the current submodel evaluator
      INPAR::BEAMINTERACTION::SubModelType Type() const override
      {
        return INPAR::BEAMINTERACTION::submodel_potential;
      }

      //! @name Derived public BEAMINTERACTION::SUBMODELEVALUATOR::Generic methods
      //! @{
      //! \brief reset submodel specific variables
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
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override;

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

     protected:
      //!@name routines that are not derived and handle beam potential-based interactions
      //! @{
      /// print
      void PrintAllBeamPotentialElementPairs(std::ostream& out) const;

      /// print
      void PrintActiveBeamPotentialSet(std::ostream& out) const;

      //! @}

     private:
      inline BEAMINTERACTION::BeamPotentialParams const& BeamPotentialParams() const
      {
        CheckInit();
        return *beam_potential_params_ptr_;
      }

      inline BEAMINTERACTION::BeamPotentialParams& BeamPotentialParams()
      {
        CheckInit();
        return *beam_potential_params_ptr_;
      }

      inline Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> BeamPotentialParamsPtr() const
      {
        CheckInit();
        return beam_potential_params_ptr_;
      }

      //!@name routines that are not derived and handle beam potential-based interactions
      //! @{
      /// get neighbouring eles in discret
      void FindAndStoreNeighboringElements();

      /// exclude certain neighbors from interaction evaluation
      void SelectElesToBeConsideredForPotentialEvaluation(
          DRT::Element* currele, std::set<DRT::Element*>& neighbors) const;

      /// create instances of class BeamContactPair that will be evaluated
      //  to get force and stiffness contributions from beam interactions
      void CreateBeamPotentialElementPairs();

      void GetBeamPotentialConditionsAppliedToThisElementPair(
          BEAMINTERACTION::BeamPotentialPair const& elementpair,
          std::vector<DRT::Condition*>& conditions_element1,
          std::vector<DRT::Condition*>& conditions_element2) const;

      //! @}

      /** \brief print this beam potential-based element pair to screen
       *
       *  \author grill */
      void PrintConsoleWelcomeMessage(std::ostream& out) const;

      //!@name routines that handle visualization output for potential-based interactions
      //! @{

      //! init output for potential-based interactions in VTP format
      void InitOutputRuntimeBeamPotential();

      //! writes VTP output for potential-based interactions at the end of a time step
      void WriteTimeStepOutputRuntimeBeamPotential() const;

      //! writes VTP output for potential-based interactions at the end of a nonlinear iteration
      void WriteIterationOutputRuntimeBeamPotential(int iteration_number) const;

      //! writes VTP output for potential-based interactions
      void WriteOutputRuntimeBeamPotential(int timestep_number, double time) const;

      //! @}

     private:
      //! data container holding all beam contact related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> beam_potential_params_ptr_;

      //! type of eles in bins  // Todo kept line for future improvement
      //    BINSTRATEGY::UTILS::BinContentType bin_beamcontent_;

      //! interacting pairs of beam elements that might exert forces on each other
      std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>> beam_potential_element_pairs_;

      //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
      std::map<int, std::set<DRT::Element*>> nearby_elements_map_;

      //! runtime vtp writer for visualization of potential-based interactions
      Teuchos::RCP<IO::VisualizationManager> visualization_manager_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
