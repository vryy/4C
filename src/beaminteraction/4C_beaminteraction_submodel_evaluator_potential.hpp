/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief submodel for potential-based beam interactions


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_submodel_evaluator_generic.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_discretization_condition.hpp"
#include "4C_io_visualization_manager.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace Core::Elements
{
  class Element;
}

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
      void post_setup() override;

      //! Returns the type of the current submodel evaluator
      Inpar::BEAMINTERACTION::SubModelType Type() const override
      {
        return Inpar::BEAMINTERACTION::submodel_potential;
      }

      //! @name Derived public BEAMINTERACTION::SUBMODELEVALUATOR::Generic methods
      //! @{
      //! \brief reset submodel specific variables
      //! derived
      void Reset() override;

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      bool pre_update_step_element(bool beam_redist) override;

      //! derived
      void UpdateStepElement(bool repartition_was_done) override;

      //! derived
      void post_update_step_element() override;

      //! derived
      std::map<STR::EnergyType, double> get_energy() const override;

      //! derived
      void OutputStepState(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void runtime_output_step_state() const override;

      //! derived
      void ResetStepState() override;

      //! derived
      void write_restart(Core::IO::DiscretizationWriter& ia_writer,
          Core::IO::DiscretizationWriter& bin_writer) const override;

      //! derived
      void PreReadRestart() override;

      //! derived
      void read_restart(Core::IO::DiscretizationReader& ia_reader,
          Core::IO::DiscretizationReader& bin_reader) override;

      //! derived
      void PostReadRestart() override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override;

      //! derived
      void init_submodel_dependencies(
          Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap) override;

      //! derived
      void AddBinsToBinColMap(std::set<int>& colbins) override;

      //! derived
      void add_bins_with_relevant_content_for_ia_discret_col_map(
          std::set<int>& colbins) const override;

      //! derived
      void get_half_interaction_distance(double& half_interaction_distance) override;

      //! @}

     protected:
      //!@name routines that are not derived and handle beam potential-based interactions
      //! @{
      /// print
      void print_all_beam_potential_element_pairs(std::ostream& out) const;

      /// print
      void print_active_beam_potential_set(std::ostream& out) const;

      //! @}

     private:
      inline BEAMINTERACTION::BeamPotentialParams const& beam_potential_params() const
      {
        check_init();
        return *beam_potential_params_ptr_;
      }

      inline BEAMINTERACTION::BeamPotentialParams& beam_potential_params()
      {
        check_init();
        return *beam_potential_params_ptr_;
      }

      inline Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> beam_potential_params_ptr() const
      {
        check_init();
        return beam_potential_params_ptr_;
      }

      //!@name routines that are not derived and handle beam potential-based interactions
      //! @{
      /// get neighbouring eles in discret
      void find_and_store_neighboring_elements();

      /// exclude certain neighbors from interaction evaluation
      void select_eles_to_be_considered_for_potential_evaluation(
          Core::Elements::Element* currele, std::set<Core::Elements::Element*>& neighbors) const;

      /// create instances of class BeamContactPair that will be evaluated
      //  to get force and stiffness contributions from beam interactions
      void create_beam_potential_element_pairs();

      void get_beam_potential_conditions_applied_to_this_element_pair(
          BEAMINTERACTION::BeamPotentialPair const& elementpair,
          std::vector<Core::Conditions::Condition*>& conditions_element1,
          std::vector<Core::Conditions::Condition*>& conditions_element2) const;

      //! @}

      /** \brief print this beam potential-based element pair to screen
       *
       *  \author grill */
      void print_console_welcome_message(std::ostream& out) const;

      //!@name routines that handle visualization output for potential-based interactions
      //! @{

      //! init output for potential-based interactions in VTP format
      void init_output_runtime_beam_potential();

      //! writes VTP output for potential-based interactions at the end of a time step
      void write_time_step_output_runtime_beam_potential() const;

      //! writes VTP output for potential-based interactions at the end of a nonlinear iteration
      void write_iteration_output_runtime_beam_potential(int iteration_number) const;

      //! writes VTP output for potential-based interactions
      void write_output_runtime_beam_potential(int timestep_number, double time) const;

      //! @}

     private:
      //! data container holding all beam contact related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> beam_potential_params_ptr_;

      //! type of eles in bins  // Todo kept line for future improvement
      //    BINSTRATEGY::UTILS::BinContentType bin_beamcontent_;

      //! interacting pairs of beam elements that might exert forces on each other
      std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>> beam_potential_element_pairs_;

      //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
      std::map<int, std::set<Core::Elements::Element*>> nearby_elements_map_;

      //! runtime vtp writer for visualization of potential-based interactions
      Teuchos::RCP<Core::IO::VisualizationManager> visualization_manager_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
