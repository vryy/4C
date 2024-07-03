/*-----------------------------------------------------------*/
/*! \file

\brief class for managing rigid sphere to beam crosslinking


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_HPP
#define FOUR_C_BEAMINTERACTION_SUBMODEL_EVALUATOR_SPHEREBEAMLINKING_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_submodel_evaluator_generic.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_MpiComm.h>

#include <unordered_set>

FOUR_C_NAMESPACE_OPEN


namespace Core::Nodes
{
  class Node;
}

namespace Core::Elements
{
  class Element;
}

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
      void setup() override;

      //! derived
      void post_setup() override;

      //! Returns the type of the current model evaluator
      Inpar::BEAMINTERACTION::SubModelType Type() const override
      {
        return Inpar::BEAMINTERACTION::submodel_spherebeamlink;
      }

      //! derived
      void reset() override;

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
      std::map<Solid::EnergyType, double> get_energy() const override;

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
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override{/*empty*/};

      //! derived
      void init_submodel_dependencies(
          Teuchos::RCP<Solid::MODELEVALUATOR::BeamInteraction::Map> const submodelmap) override;

      //! derived
      void AddBinsToBinColMap(std::set<int>& colbins) override;

      //! derived
      void add_bins_with_relevant_content_for_ia_discret_col_map(
          std::set<int>& colbins) const override;

      //! derived
      void get_half_interaction_distance(double& half_interaction_distance) override;

      //! @}

     private:
      //! @name submodel specific, not derived member functions
      //! @{

      //! init output in VTP format
      void init_output_runtime();

      //! writes output for discretization structure in VTP format
      void write_output_runtime() const;

      /// get neighboring eles in discret
      virtual void find_and_store_neighboring_elements(
          std::map<int, std::vector<std::pair<int, int>>>& newlinks);

      /// get neighboring eles in discret
      virtual void check_feasibility_of_new_link(Core::Elements::Element const* currele,
          std::vector<Core::Elements::Element const*> const& neighbors,
          std::unordered_set<int>& tobebonded,
          std::map<int, std::vector<std::pair<int, int>>>& newlinks) const;

      /// create new beam to sphere joint object
      virtual void create_beam_to_sphere_joint(
          std::map<int, std::vector<std::pair<int, int>>> const& newlinks);

      /// update linker length to mimic contractive cells
      virtual void update_linker_length();

      /// check if bond needs to be dissolved
      virtual void unbind_sphere_beam_bonds(int& num_dissolved);

      /// compute force dependent off rate for a catch-slip bond
      virtual void calc_force_dependent_catch_slip_bond_unbind_probability(
          Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr, double& p_unbind);

      //! @}

     private:
      //! cell discretization
      Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking> sm_crosslinkink_ptr_;

      //! data container holding all beam contact related parameters
      Teuchos::RCP<BEAMINTERACTION::SphereBeamLinkingParams> spherebeamlinking_params_ptr_;

      //! runtime output for cell beam crosslinks (integrins)
      Teuchos::RCP<Core::IO::VisualizationManager> visualization_manager_ptr_;

      //! step number for random stuff concerning sphere beam linking
      int random_number_sphere_beam_linking_step_;
    };

  }  // namespace SUBMODELEVALUATOR
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
