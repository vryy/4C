/*-----------------------------------------------------------*/
/*! \file

\brief class for submodel beam contact


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_submodel_evaluator_generic.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_io_visualization_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}
namespace CORE::GEOMETRICSEARCH
{
  class geometric_search_params;
  class GeometricSearchVisualization;
}  // namespace CORE::GEOMETRICSEARCH
namespace BEAMINTERACTION
{
  class BeamContactParams;
  class BeamContactPair;
  class BeamInteractionConditions;
  class BeamToSolidVolumeMeshtyingVisualizationOutputWriter;
  class BeamToSolidSurfaceVisualizationOutputWriter;

  namespace SUBMODELEVALUATOR
  {
    class BeamContactAssemblyManager;

    class BeamContact : public Generic
    {
     public:
      //! constructor
      BeamContact();

      //! setup class variables
      void Setup() override;

      //! derived
      void PostSetup() override;

      //! Returns the type of the current submodel evaluator
      INPAR::BEAMINTERACTION::SubModelType Type() const override
      {
        return INPAR::BEAMINTERACTION::submodel_beamcontact;
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

      /**
       * \brief Call PreEvaluate on each pair.
       *
       * This has to be done outside of the the assembly managers, as PreEvaluate should only be
       * called once per pair, and it is technically possible that the same pair is in multiple
       * assembly managers.
       */
      void PreEvaluate();

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      bool pre_update_step_element(bool beam_redist) override;

      //! derived
      void UpdateStepElement(bool repartition_was_done) override;

      //! derived
      void post_update_step_element() override;

      //! derived
      std::map<STR::EnergyType, double> GetEnergy() const override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void runtime_output_step_state() const override;

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
      void init_submodel_dependencies(
          Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap) override;

      //! derived
      void AddBinsToBinColMap(std::set<int>& colbins) override;

      //! derived
      void add_bins_with_relevant_content_for_ia_discret_col_map(
          std::set<int>& colbins) const override;

      //! derived
      void get_half_interaction_distance(double& half_interaction_distance) override;

      /**
       * \brief Return the assembly managers in this submodel evaluator.
       */
      inline const std::vector<
          Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManager>>&
      GetAssemblyManagers() const
      {
        return assembly_managers_;
      }

      /**
       * \brief Return the geometry pairs in this submodel evaluator.
       */
      inline const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& GetContactPairs()
          const
      {
        return contact_elepairs_;
      }

      /**
       * \brief Return the conditions in this submodel evaluator.
       */
      inline Teuchos::RCP<const BEAMINTERACTION::BeamInteractionConditions> GetConditions() const
      {
        return beam_interaction_conditions_ptr_;
      }

      //! @}

     protected:
      //!@name routines that are not derived and handle beam to beam contact
      //! @{
      /// check which contact is considered
      bool HaveContactType(BINSTRATEGY::UTILS::BinContentType const& contacttype) const;

      /// print
      void print_all_beam_contact_element_pairs(std::ostream& out) const;

      /// print
      void print_active_beam_contact_set(std::ostream& out) const;

      //! @}

     private:
      inline BEAMINTERACTION::BeamContactParams const& BeamContactParams() const
      {
        CheckInit();
        return *beam_contact_params_ptr_;
      }

      inline BEAMINTERACTION::BeamContactParams& BeamContactParams()
      {
        CheckInit();
        return *beam_contact_params_ptr_;
      }

      inline Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_ptr() const
      {
        CheckInit();
        return beam_contact_params_ptr_;
      }

      //!@name routines that handle visualization output for contact forces
      //! @{

      //! init output for contact forces in VTP format
      void init_output_runtime_beam_contact();

      //! writes visualization output for contact forces at the end of a time step
      void write_time_step_output_runtime_beam_contact() const;

      //! writes visualization output for contact forces at the end of a nonlinear iteration
      void write_iteration_output_runtime_beam_contact(int iteration_number) const;

      //! writes visualization output for contact forces
      void write_output_runtime_beam_contact(int timestep_number, double time) const;


      //! @}

      //!@name routines that are not derived and handle beam to beam contact
      //! @{

      // init element types considered for beam to ? contact
      void init_element_types_considered_for_contact();

      /// get neighbouring eles in discret
      void find_and_store_neighboring_elements();

      /// exclude certain neighbors from interaction evaluation
      /*
       * Here we remove pairs based on very trivial rules, e.g.:
       * - an element can not be in contact with itself, thus is removed if same global indice
       *   is found
       * - evaluate beam contact pairs only once i.e., pair(ele1, ele2) and pair(ele2, ele1),
       *   thus only beam-to-beam pairs where the GID of the first element is smaller than the
       *   second are considered
       * - elements sharing nodes don't interact
       */
      void select_eles_to_be_considered_for_contact_evaluation(
          DRT::Element* currele, std::set<DRT::Element*>& neighbors) const;

      /// create instances of class BeamContactPair that will be evaluated
      //  to get force and stiffness contributions from beam interactions
      void create_beam_contact_element_pairs();

      /// Add the restart displacement to the pairs, if the coupling should be evaluated with
      /// respect to the restart state.
      void set_restart_displacement_in_pairs();

      //! @}

      //! data container holding all beam contact related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_ptr_;

      //! data container holding all beam interaction related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamInteractionParams> beam_interaction_params_ptr_;

      //! data container holding all beam interactions defined by conditions
      Teuchos::RCP<BEAMINTERACTION::BeamInteractionConditions> beam_interaction_conditions_ptr_;

      //! data container holding all geometric search related parameters
      Teuchos::RCP<CORE::GEOMETRICSEARCH::GeometricSearchParams> geometric_search_params_ptr_;

      //! element types considered for beam to ? contact
      std::vector<BINSTRATEGY::UTILS::BinContentType> contactelementtypes_;

      //! interacting pairs of beam elements that might exert forces on each other
      std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> contact_elepairs_;

      //! Objects to evaluate system contributions for stiffness and force terms.
      std::vector<Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManager>>
          assembly_managers_;

      //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
      std::map<int, std::set<DRT::Element*>> nearby_elements_map_;

      //! runtime visualization writer for visualization of contact forces
      Teuchos::RCP<IO::VisualizationManager> visualization_manager_ptr_;

      //! This object handles all beam to solid volume related visualization output.
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter>
          beam_to_solid_volume_meshtying_visualization_output_writer_ptr_;

      //! This object handles all beam to solid surface related visualization output.
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter>
          beam_to_solid_surface_visualization_output_writer_ptr_;

      //! This object handles all geometric search related visualization output.
      Teuchos::RCP<CORE::GEOMETRICSEARCH::GeometricSearchVisualization>
          geometric_search_visualization_ptr_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
