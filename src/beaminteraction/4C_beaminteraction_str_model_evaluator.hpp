/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of all beam interaction terms


\level 3
*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_STR_MODEL_EVALUATOR_HPP
#define FOUR_C_BEAMINTERACTION_STR_MODEL_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"  // base class

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace Adapter
{
  class Coupling;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class MatrixRowTransform;
}  // namespace Core::LinAlg

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
      typedef std::map<enum Inpar::BEAMINTERACTION::SubModelType,
          Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic>>
          Map;
      typedef std::vector<Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic>> Vector;

      //! constructor
      BeamInteraction();

      void setup() override;

      virtual void post_setup();

      /// print welcome to biopolymer network simulation
      virtual void Logo() const;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{
      //! derived

      //! derived
      Inpar::STR::ModelType Type() const override { return Inpar::STR::model_beaminteraction; }

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override { return; };

      //! derived
      void post_evaluate() override{/* currently unused */};

      //! derived
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! derived
      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const Inpar::STR::PredEnum& pred_type) override { return; };

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override;

      //! derived
      void update_step_state(const double& timefac_n) override;

      //! derived
      void update_step_element() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void determine_energy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void runtime_output_step_state() const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! derived
      void post_output() override;

      //! derived
      void reset_step_state() override;
      //! @}

     protected:
      //! derived
      void reset(const Epetra_Vector& x) override;

      //!@name routines for submodel management
      //! @{

     public:
      /// check if the given model type is active.
      bool HaveSubModelType(Inpar::BEAMINTERACTION::SubModelType const& submodeltype) const;

     private:
      void partition_problem();

      bool post_partition_problem();

      //! set beaminteraction sub models
      void set_sub_model_types();


      //! build, init and setup submodel evaluator
      void init_and_setup_sub_model_evaluators();

      //! give submodels a certain order in which they are evaluated
      virtual Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> transform_to_vector(
          STR::MODELEVALUATOR::BeamInteraction::Map submodel_map,
          std::vector<Inpar::BEAMINTERACTION::SubModelType>& sorted_submodel_types) const;

      //! @}

      //!@name routines that manage to discretizations with distinct parallel distribution
      //! @{

      /// check if interaction discretization needs to be redistributed completely
      bool check_if_beam_discret_redistribution_needs_to_be_done();

      /// update coupling adapter and matrix transformation object with new maps
      void update_coupling_adapter_and_matrix_transformation();

      /// transform force vector from ia_discret_ to discret()
      virtual void transform_force();

      /// transform stiffness matrix from ia_discret_ to discret()
      virtual void transform_stiff();

      /// transform force vector and stiffness matrix from ia_discret_ to discret()
      virtual void transform_force_stiff();

      /// update states based on bindis after its redistribution
      virtual void update_maps();

      //! @}

      //!@name routines that manage binning strategy
      //! @{

      /// change parallel distribution of bindis and ia_discret and assign (beam) eles to bins
      virtual void extend_ghosting();

      /// build ele to bin map
      virtual void build_row_ele_to_bin_map();

      /// print some information about binning
      virtual void print_binning_info_to_screen() const;

      //! @}

     private:
      //! pointer to the problem discretization (cast of base class member)
      Teuchos::RCP<Core::FE::Discretization> discret_ptr_;

      //! data container holding all beaminteraction related parameters
      Teuchos::RCP<BEAMINTERACTION::BeamInteractionParams> beaminteraction_params_ptr_;

      //!@name data for submodel management
      //! @{
      /// current active model types for the model evaluator
      Teuchos::RCP<std::set<enum Inpar::BEAMINTERACTION::SubModelType>> submodeltypes_;

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> me_map_ptr_;

      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> me_vec_ptr_;
      //! @}

      //!@name data for handling two distinct parallel distributed discretizations
      //! @{
      //! myrank
      int myrank_;

      //! coupling adapter to transfer vectors and matrices between discret() and intactids_
      Teuchos::RCP<Core::Adapter::Coupling> coupsia_;

      //! transform object for structure stiffness matrix
      Teuchos::RCP<Core::LinAlg::MatrixRowTransform> siatransform_;
      //! @}


      //!@name data for beaminteraction with binning strategy
      //! @{
      //! interaction discretization handling all interactions (e.g. crosslinker to beam,
      //! beam to beam, potential ...)
      Teuchos::RCP<Core::FE::Discretization> ia_discret_;

      /// map extractor for split of different element types
      Teuchos::RCP<Core::LinAlg::MultiMapExtractor> eletypeextractor_;

      //! pointer to the global state data container
      Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> ia_state_ptr_;

      //! force based on ia_discret at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> ia_force_beaminteraction_;

      //! global force based on discret() at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> force_beaminteraction_;

      //! structural stiffness matrix based on discret()
      Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_beaminteraction_;

      //! beam crosslinker handler
      Teuchos::RCP<BEAMINTERACTION::BeamCrosslinkerHandler> beam_crosslinker_handler_;

      //! binning strategy
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

      //! crosslinker and bin discretization
      Teuchos::RCP<Core::FE::Discretization> bindis_;

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


FOUR_C_NAMESPACE_CLOSE

#endif
