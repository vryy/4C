/*-----------------------------------------------------------*/
/*! \file

\brief class for submodel beam contact


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_submodel_evaluator_beamcontact.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_writer.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_writer.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_conditions.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_beaminteraction_contact_runtime_visualization_output_params.hpp"
#include "4C_beaminteraction_data.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.hpp"
#include "4C_binstrategy.hpp"
#include "4C_fem_geometric_search_bvh.hpp"
#include "4C_fem_geometric_search_params.hpp"
#include "4C_fem_geometric_search_visualization.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_inpar_geometry_pair.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_rigidsphere.hpp"
#include "4C_so3_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_FEVector.h>
#include <NOX_Solver_Generic.H>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::BeamContact()
    : beam_contact_params_ptr_(Teuchos::null),
      beam_interaction_params_ptr_(Teuchos::null),
      beam_interaction_conditions_ptr_(Teuchos::null),
      geometric_search_params_ptr_(Teuchos::null),
      contact_elepairs_(Teuchos::null),
      assembly_managers_(Teuchos::null),
      beam_to_solid_volume_meshtying_visualization_output_writer_ptr_(Teuchos::null),
      beam_to_solid_surface_visualization_output_writer_ptr_(Teuchos::null)
{
  // clear stl stuff
  nearby_elements_map_.clear();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::setup()
{
  check_init();

  // build a new data container to manage beam interaction parameters
  beam_interaction_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::BeamInteractionParams());
  beam_interaction_params_ptr_->init();
  beam_interaction_params_ptr_->setup();

  // build a new data container to manage geometric search parameters
  geometric_search_params_ptr_ = Teuchos::rcp(new Core::GeometricSearch::GeometricSearchParams(
      Global::Problem::Instance()->geometric_search_params(),
      Global::Problem::Instance()->IOParams()));
  if (beam_interaction_params_ptr_->GetSearchStrategy() ==
          Inpar::BEAMINTERACTION::SearchStrategy::bounding_volume_hierarchy &&
      geometric_search_params_ptr_->get_write_visualization_flag())
  {
    geometric_search_visualization_ptr_ =
        Teuchos::rcp(new Core::GeometricSearch::GeometricSearchVisualization(
            Core::IO::VisualizationParametersFactory(
                Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                *Global::Problem::Instance()->OutputControlFile(), GState().get_time_n()),
            DiscretPtr()->Comm(), "beam-interaction-geometric-search"));
  }

  // build a new data container to manage beam contact parameters
  beam_contact_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::BeamContactParams());

  // build runtime visualization writer if desired
  if ((bool)Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->beam_contact_params().sublist("RUNTIME VTK OUTPUT"),
          "VTK_OUTPUT_BEAM_CONTACT"))
  {
    beam_contact_params_ptr_->build_beam_contact_runtime_output_params(GState().get_time_n());

    init_output_runtime_beam_contact();
  }


  contactelementtypes_.clear();

  if (Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO BEAM CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::Beam);

    beam_contact_params_ptr_->build_beam_to_beam_contact_params();
  }

  // conditions for beam penalty point coupling
  std::vector<Core::Conditions::Condition*> beampenaltycouplingconditions(0);
  Discret().GetCondition("PenaltyPointCouplingCondition", beampenaltycouplingconditions);
  if (beampenaltycouplingconditions.size() > 0)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::Beam);
  }

  if (Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO SPHERE CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::RigidSphere);

    beam_contact_params_ptr_->build_beam_to_sphere_contact_params();
  }

  // Check if beam-to-solid volume mesh tying is present.
  const Teuchos::ParameterList& beam_to_solid_volume_parameters =
      Global::Problem::Instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID VOLUME MESHTYING");
  if (Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          beam_to_solid_volume_parameters, "CONTACT_DISCRETIZATION") !=
      Inpar::BeamToSolid::BeamToSolidContactDiscretization::none)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::Solid);

    beam_contact_params_ptr_->build_beam_to_solid_volume_meshtying_params();

    // Build the beam to solid volume meshtying output writer if desired.
    if (beam_contact_params_ptr_->beam_to_solid_volume_meshtying_params()
            ->get_visualization_output_params_ptr()
            ->GetOutputFlag())
    {
      beam_to_solid_volume_meshtying_visualization_output_writer_ptr_ =
          Teuchos::rcp<BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter>(
              new BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter(
                  Core::IO::VisualizationParametersFactory(
                      Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                      *Global::Problem::Instance()->OutputControlFile(), GState().get_time_n())));
      beam_to_solid_volume_meshtying_visualization_output_writer_ptr_->init();
      beam_to_solid_volume_meshtying_visualization_output_writer_ptr_->setup(
          GInOutput().get_runtime_output_params(),
          beam_contact_params_ptr_->beam_to_solid_volume_meshtying_params()
              ->get_visualization_output_params_ptr());
    }
  }

  // Check if beam-to-solid surface mesh tying is present.
  const Teuchos::ParameterList& beam_to_solid_surface_parameters =
      Global::Problem::Instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID SURFACE MESHTYING");
  if (Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          beam_to_solid_surface_parameters, "CONTACT_DISCRETIZATION") !=
      Inpar::BeamToSolid::BeamToSolidContactDiscretization::none)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::Solid);

    beam_contact_params_ptr_->build_beam_to_solid_surface_meshtying_params();

    // Build the beam to solid surface output writer if desired.
    if (beam_contact_params_ptr_->beam_to_solid_surface_meshtying_params()
            ->get_visualization_output_params_ptr()
            ->GetOutputFlag())
    {
      beam_to_solid_surface_visualization_output_writer_ptr_ =
          Teuchos::rcp<BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter>(
              new BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter(
                  Core::IO::VisualizationParametersFactory(
                      Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                      *Global::Problem::Instance()->OutputControlFile(), GState().get_time_n())));
      beam_to_solid_surface_visualization_output_writer_ptr_->init();
      beam_to_solid_surface_visualization_output_writer_ptr_->setup(
          GInOutput().get_runtime_output_params(),
          beam_contact_params_ptr_->beam_to_solid_surface_meshtying_params()
              ->get_visualization_output_params_ptr());
    }
  }

  // Check if beam-to-solid surface contact is present.
  const Teuchos::ParameterList& beam_to_solid_surface_contact_parameters =
      Global::Problem::Instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID SURFACE CONTACT");
  if (Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          beam_to_solid_surface_contact_parameters, "CONTACT_DISCRETIZATION") !=
      Inpar::BeamToSolid::BeamToSolidContactDiscretization::none)
  {
    contactelementtypes_.push_back(Core::Binstrategy::Utils::BinContentType::Solid);

    beam_contact_params_ptr_->build_beam_to_solid_surface_contact_params();
  }

  // Build the container to manage beam-to-solid conditions and get all coupling conditions.
  beam_interaction_conditions_ptr_ = Teuchos::rcp(new BEAMINTERACTION::BeamInteractionConditions());
  beam_interaction_conditions_ptr_->set_beam_interaction_conditions(
      DiscretPtr(), beam_contact_params_ptr_);

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::post_setup()
{
  check_init_setup();

  // Todo really needed here? maybe find better place
  // ensure that contact is evaluated correctly at beginning of first time step (initial overlap)
  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_contact_element_pairs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::init_submodel_dependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap)
{
  check_init_setup();
  // no active influence on other submodels
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::reset()
{
  check_init_setup();

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>::const_iterator iter;
  for (iter = contact_elepairs_.begin(); iter != contact_elepairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> elepairptr = *iter;

    std::vector<const Core::Elements::Element*> element_ptr(2);

    element_ptr[0] = elepairptr->Element1();
    element_ptr[1] = elepairptr->Element2();

    // element Dof values relevant for centerline interpolation
    std::vector<std::vector<double>> element_posdofvec_absolutevalues(2);

    for (unsigned int ielement = 0; ielement < 2; ++ielement)
    {
      // extract the Dof values of this element from displacement vector
      BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(Discret(), element_ptr[ielement],
          beam_interaction_data_state_ptr()->GetDisColNp(),
          element_posdofvec_absolutevalues[ielement]);
    }

    // update positional Dof values in the interaction element pair object
    elepairptr->ResetState(
        element_posdofvec_absolutevalues[0], element_posdofvec_absolutevalues[1]);

    // update rotational Dof values in the interaction pair object
    elepairptr->ResetRotationState(Discret(), beam_interaction_data_state_ptr()->GetDisColNp());
  }

  // Set restart displacements in the pairs.
  set_restart_displacement_in_pairs();

  // Update the geometry pair evaluation data.
  beam_interaction_conditions_ptr_->set_state(DiscretPtr(), beam_interaction_data_state_ptr());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::evaluate_force()
{
  check_init_setup();

  pre_evaluate();

  // Loop over the assembly manager and assemble contributions into the global force vector.
  for (auto& assembly_manager : assembly_managers_)
  {
    assembly_manager->evaluate_force_stiff(DiscretPtr(), beam_interaction_data_state_ptr(),
        beam_interaction_data_state_ptr()->GetForceNp(), Teuchos::null);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::evaluate_stiff()
{
  check_init_setup();

  pre_evaluate();

  // Loop over the assembly manager and assemble contributions into the global stiffness matrix.
  for (auto& assembly_manager : assembly_managers_)
  {
    assembly_manager->evaluate_force_stiff(DiscretPtr(), beam_interaction_data_state_ptr(),
        Teuchos::null, beam_interaction_data_state_ptr()->GetStiff());
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::evaluate_force_stiff()
{
  check_init_setup();

  pre_evaluate();

  // Loop over the assembly manager and assemble contributions into the global force vector and
  // stiffness matrix.
  for (auto& assembly_manager : assembly_managers_)
    assembly_manager->evaluate_force_stiff(DiscretPtr(), beam_interaction_data_state_ptr(),
        beam_interaction_data_state_ptr()->GetForceNp(),
        beam_interaction_data_state_ptr()->GetStiff());

  print_active_beam_contact_set(Core::IO::cout.os(Core::IO::verbose));

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::pre_evaluate()
{
  for (auto& elepairptr : contact_elepairs_) elepairptr->pre_evaluate();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::UpdateStepState(const double& timefac_n)
{
  check_init_setup();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::pre_update_step_element(bool beam_redist)
{
  check_init_setup();

  /* Fixme
   * writing vtk output needs to be done BEFORE updating (and thus clearing
   * element pairs)
   * move this to runtime_output_step_state as soon as we keep element pairs
   * from previous time step */
  /* Fixme
   * writing this output also must be done BEFORE re-distribution which
   * currently happens in STR::MODELEVALUATOR::BeamInteraction::update_step_element()
   * before calling update_step_element() on all submodels.
   * Hence, the only option currently is to call it from pre_update_step_element() */
  /* Note: another option would be to not use any data from state vectors or elements and only
   * write previously computed and (locally) stored data at this point. Like
   * this, it works in SUBMODELEVALUATOR::BeamPotential */
  if (visualization_manager_ptr_ != Teuchos::null and
      GState().get_step_np() % beam_contact_params()
                                   .beam_contact_runtime_visualization_output_params()
                                   ->output_interval_in_steps() ==
          0)
  {
    write_time_step_output_runtime_beam_contact();
  }
  if (beam_to_solid_volume_meshtying_visualization_output_writer_ptr_ != Teuchos::null)
    beam_to_solid_volume_meshtying_visualization_output_writer_ptr_->write_output_runtime(this);
  if (beam_to_solid_surface_visualization_output_writer_ptr_ != Teuchos::null)
    beam_to_solid_surface_visualization_output_writer_ptr_->write_output_runtime(this);

  // not repartition of binning discretization necessary
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::UpdateStepElement(bool repartition_was_done)
{
  check_init_setup();

  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_contact_element_pairs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::post_update_step_element()
{
  check_init_setup();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<STR::EnergyType, double> BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::get_energy()
    const
{
  check_init_setup();

  std::map<STR::EnergyType, double> contact_penalty_potential;

  for (auto& elepairptr : contact_elepairs_)
  {
    contact_penalty_potential[STR::beam_contact_penalty_potential] += elepairptr->get_energy();
  }

  for (const auto& assembly_manager : assembly_managers_)
  {
    contact_penalty_potential[STR::beam_contact_penalty_potential] +=
        assembly_manager->get_energy(GState().get_dis_np());
  }

  return contact_penalty_potential;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::OutputStepState(
    Core::IO::DiscretizationWriter& iowriter) const
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::runtime_output_step_state() const {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::init_output_runtime_beam_contact()
{
  check_init();

  visualization_manager_ptr_ = Teuchos::rcp(
      new Core::IO::VisualizationManager(beam_contact_params()
                                             .beam_contact_runtime_visualization_output_params()
                                             ->get_visualization_parameters(),
          Discret().Comm(), "beam-contact"));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::write_time_step_output_runtime_beam_contact()
    const
{
  check_init_setup();

  auto [output_time, output_step] = Core::IO::GetTimeAndTimeStepIndexForOutput(
      beam_contact_params()
          .beam_contact_runtime_visualization_output_params()
          ->get_visualization_parameters(),
      GState().get_time_n(), GState().get_step_n());
  write_output_runtime_beam_contact(output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::write_iteration_output_runtime_beam_contact(
    int iteration_number) const
{
  check_init_setup();

  auto [output_time, output_step] = Core::IO::GetTimeAndTimeStepIndexForOutput(
      beam_contact_params()
          .beam_contact_runtime_visualization_output_params()
          ->get_visualization_parameters(),
      GState().get_time_n(), GState().get_step_n(), iteration_number);
  write_output_runtime_beam_contact(output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::write_output_runtime_beam_contact(
    int timestep_number, double time) const
{
  check_init_setup();

  const unsigned int num_spatial_dimensions = 3;

  // number of active point contact point pairs * 2 = number of row points for writer object
  unsigned int num_row_points = 0;

  // loop over contact pairs and retrieve number of all active contact point pairs
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>::const_iterator pair_iter;
  for (pair_iter = contact_elepairs_.begin(); pair_iter != contact_elepairs_.end(); ++pair_iter)
  {
    num_row_points += 2 * (*pair_iter)->get_num_all_active_contact_point_pairs();
  }

  // get and prepare storage for point coordinate values
  std::vector<double>& point_coordinates =
      visualization_manager_ptr_->get_visualization_data().GetPointCoordinates();
  point_coordinates.clear();
  point_coordinates.reserve(num_spatial_dimensions * num_row_points);

  // contact force values: collect data and append to visualization results if desired
  std::vector<double> contact_force_vector(0);
  contact_force_vector.reserve(num_spatial_dimensions * num_row_points);

  // gap values: collect data and append to visualization results if desired
  std::vector<double> gaps(0);
  gaps.reserve(num_row_points);

  // loop over my points and collect the geometry/grid data, i.e. contact points
  std::vector<Core::LinAlg::Matrix<3, 1, double>> coordinates_ele1_this_pair;
  std::vector<Core::LinAlg::Matrix<3, 1, double>> coordinates_ele2_this_pair;

  std::vector<double> contact_forces_this_pair;
  std::vector<double> gaps_this_pair;

  // loop over contact pairs and retrieve all active contact point coordinates
  for (pair_iter = contact_elepairs_.begin(); pair_iter != contact_elepairs_.end(); ++pair_iter)
  {
    if ((*pair_iter)->GetContactFlag() == true)
    {
      // active contact points of element 1 and element 2
      (*pair_iter)->get_all_active_contact_point_coords_element1(coordinates_ele1_this_pair);
      (*pair_iter)->get_all_active_contact_point_coords_element2(coordinates_ele2_this_pair);
      (*pair_iter)->get_all_active_contact_forces(contact_forces_this_pair);
      (*pair_iter)->get_all_active_contact_gaps(gaps_this_pair);

      const unsigned int num_active_point_pairs = (unsigned int)coordinates_ele1_this_pair.size();

      FOUR_C_ASSERT(num_active_point_pairs == (unsigned int)coordinates_ele2_this_pair.size(),
          "number of active points on element 1 does not match number of active points "
          "on element 2!");

      FOUR_C_ASSERT(num_active_point_pairs == (unsigned int)contact_forces_this_pair.size(),
          "number of active points on element 1 does not match number of contact forces!");

      FOUR_C_ASSERT(num_active_point_pairs == (unsigned int)gaps_this_pair.size(),
          "number of active points on element 1 does not match number of gap values!");


      for (unsigned int ipointpair = 0; ipointpair < num_active_point_pairs; ++ipointpair)
      {
        Core::LinAlg::Matrix<3, 1, double> normal_vector;
        normal_vector.update(1.0, coordinates_ele1_this_pair[ipointpair], -1.0,
            coordinates_ele2_this_pair[ipointpair]);

        // contact point on first element
        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele1_this_pair[ipointpair](idim));

          contact_force_vector.push_back(
              contact_forces_this_pair[ipointpair] * normal_vector(idim));
        }
        gaps.push_back(gaps_this_pair[ipointpair]);

        // contact point on second element
        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele2_this_pair[ipointpair](idim));

          contact_force_vector.push_back(
              -1.0 * contact_forces_this_pair[ipointpair] * normal_vector(idim));
        }
        gaps.push_back(gaps_this_pair[ipointpair]);
      }
    }
  }


  // append all desired output data to the writer object's storage
  if (beam_contact_params()
          .beam_contact_runtime_visualization_output_params()
          ->is_write_contact_forces())
  {
    visualization_manager_ptr_->get_visualization_data().SetPointDataVector(
        "force", contact_force_vector, num_spatial_dimensions);
  }
  if (beam_contact_params().beam_contact_runtime_visualization_output_params()->IsWriteGaps())
  {
    visualization_manager_ptr_->get_visualization_data().SetPointDataVector("gap", gaps, 1);
  }

  // finalize everything and write all required vtk files to filesystem
  visualization_manager_ptr_->WriteToDisk(time, timestep_number);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::ResetStepState() { check_init_setup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::write_restart(
    Core::IO::DiscretizationWriter& ia_writer, Core::IO::DiscretizationWriter& bin_writer) const
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::PreReadRestart()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::read_restart(
    Core::IO::DiscretizationReader& ia_reader, Core::IO::DiscretizationReader& bin_reader)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::PostReadRestart()
{
  check_init_setup();
  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_contact_element_pairs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::run_post_iterate(
    const ::NOX::Solver::Generic& solver)
{
  check_init_setup();

  if (visualization_manager_ptr_ != Teuchos::null and
      beam_contact_params()
          .beam_contact_runtime_visualization_output_params()
          ->output_every_iteration())
  {
    write_iteration_output_runtime_beam_contact(solver.getNumIterations());
  }
  if (beam_to_solid_volume_meshtying_visualization_output_writer_ptr_ != Teuchos::null)
  {
    beam_to_solid_volume_meshtying_visualization_output_writer_ptr_->write_output_runtime_iteration(
        this, solver.getNumIterations());
  }
  if (beam_to_solid_surface_visualization_output_writer_ptr_ != Teuchos::null)
  {
    beam_to_solid_surface_visualization_output_writer_ptr_->write_output_runtime_iteration(
        this, solver.getNumIterations());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::AddBinsToBinColMap(std::set<int>& colbins)
{
  check_init_setup();
  // nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::
    add_bins_with_relevant_content_for_ia_discret_col_map(std::set<int>& colbins) const
{
  check_init_setup();
  // nothing to do
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::get_half_interaction_distance(
    double& half_interaction_distance)
{
  check_init_setup();

  // todo: choose meaningful safety factor
  double safe_fac = 1.5;

  // loop over all beams to get largest interaction radius
  double locmax_ia_distance = 0.0;
  double curr_ia_distance = 0.0;
  int const numroweles = ele_type_map_extractor_ptr()->beam_map()->NumMyElements();
  for (int rowele_i = 0; rowele_i < numroweles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->beam_map()->GID(rowele_i);
    Discret::ELEMENTS::Beam3Base* currele =
        dynamic_cast<Discret::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

    curr_ia_distance = currele->get_circular_cross_section_radius_for_interactions();

    if (curr_ia_distance > locmax_ia_distance) locmax_ia_distance = curr_ia_distance;
  }

  // get global maximum
  double globalmax_beam_ia_distance = 0.0;
  // build sum over all procs
  MPI_Allreduce(&locmax_ia_distance, &globalmax_beam_ia_distance, 1, MPI_DOUBLE, MPI_MAX,
      dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());

  // i) beam to beam contact
  if (have_contact_type(Core::Binstrategy::Utils::BinContentType::Beam))
  {
    // safety factor
    globalmax_beam_ia_distance *= safe_fac;

    half_interaction_distance = (globalmax_beam_ia_distance > half_interaction_distance)
                                    ? globalmax_beam_ia_distance
                                    : half_interaction_distance;

    // some screen output
    if (GState().get_my_rank() == 0)
      std::cout << " beam to beam contact half interaction distance " << globalmax_beam_ia_distance
                << std::endl;
  }

  // ii) beam to sphere contact
  if (have_contact_type(Core::Binstrategy::Utils::BinContentType::RigidSphere))
  {
    // loop over all spheres
    double curr_ia_dist = 0.0;
    double loc_max_ia_dist = 0.0;
    int unsigned const numrowsphereeles = EleTypeMapExtractor().sphere_map()->NumMyElements();
    for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
    {
      int const elegid = EleTypeMapExtractor().sphere_map()->GID(rowele_i);
      Discret::ELEMENTS::Rigidsphere* sphere =
          dynamic_cast<Discret::ELEMENTS::Rigidsphere*>(Discret().gElement(elegid));

      curr_ia_dist = sphere->Radius() + globalmax_beam_ia_distance;

      // update distance
      loc_max_ia_dist = (curr_ia_dist > loc_max_ia_dist) ? curr_ia_dist : loc_max_ia_dist;
    }

    // get global maximum
    double spherebeamlinking_half_interaction_distance_global = 0.0;
    // build sum over all procs
    MPI_Allreduce(&loc_max_ia_dist, &spherebeamlinking_half_interaction_distance_global, 1,
        MPI_DOUBLE, MPI_MAX, dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());

    half_interaction_distance =
        (spherebeamlinking_half_interaction_distance_global > half_interaction_distance)
            ? spherebeamlinking_half_interaction_distance_global
            : half_interaction_distance;

    // some screen output
    if (GState().get_my_rank() == 0)
      Core::IO::cout(Core::IO::verbose)
          << " sphere to beam contact half interaction distance "
          << spherebeamlinking_half_interaction_distance_global << Core::IO::endl;
  }

  // iii) beam to solid contact
  if (have_contact_type(Core::Binstrategy::Utils::BinContentType::Solid))
  {
    FOUR_C_THROW("Not yet implemented for beam to solid contact");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::have_contact_type(
    Core::Binstrategy::Utils::BinContentType const& contacttype) const
{
  check_init();
  return (std::find(contactelementtypes_.begin(), contactelementtypes_.end(), contacttype) !=
          contactelementtypes_.end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::find_and_store_neighboring_elements()
{
  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::find_and_store_neighboring_elements");

  check_init();

  // Build the ids of the elements for the beam-to-solid conditions.
  beam_interaction_conditions_ptr_->BuildIdSets(DiscretPtr());

  if (beam_interaction_params_ptr_->GetSearchStrategy() ==
      Inpar::BEAMINTERACTION::SearchStrategy::bruteforce_with_binning)
  {
    // loop over all row beam elements
    // note: like this we ensure that first element of pair is always a beam element, also only
    // beam to something contact considered
    int const numroweles = ele_type_map_extractor_ptr()->beam_map()->NumMyElements();
    for (int rowele_i = 0; rowele_i < numroweles; ++rowele_i)
    {
      int const elegid = ele_type_map_extractor_ptr()->beam_map()->GID(rowele_i);
      Core::Elements::Element* currele = DiscretPtr()->gElement(elegid);

      // (unique) set of neighboring bins for all col bins assigned to current element
      std::set<int> neighboring_binIds;

      // loop over all bins touched by currele
      std::set<int>::const_iterator biniter;
      for (biniter = beam_interaction_data_state_ptr()->GetRowEleToBinSet(elegid).begin();
           biniter != beam_interaction_data_state_ptr()->GetRowEleToBinSet(elegid).end(); ++biniter)
      {
        std::vector<int> loc_neighboring_binIds;
        // in three-dimensional space: 26 neighbouring-bins + 1 self
        loc_neighboring_binIds.reserve(27);

        // do not check on existence here -> shifted to GetBinContent
        BinStrategyPtr()->get_neighbor_and_own_bin_ids(*biniter, loc_neighboring_binIds);

        // build up comprehensive unique set of neighboring bins
        neighboring_binIds.insert(loc_neighboring_binIds.begin(), loc_neighboring_binIds.end());
      }
      // get set of elements that reside in neighboring bins
      std::vector<int> glob_neighboring_binIds(
          neighboring_binIds.begin(), neighboring_binIds.end());
      std::set<Core::Elements::Element*> neighboring_elements;
      BinStrategyPtr()->GetBinContent(
          neighboring_elements, contactelementtypes_, glob_neighboring_binIds);

      // sort out elements that should not be considered in contact evaluation
      select_eles_to_be_considered_for_contact_evaluation(currele, neighboring_elements);

      nearby_elements_map_[elegid] = neighboring_elements;
    }
  }
  else if (beam_interaction_params_ptr_->GetSearchStrategy() ==
           Inpar::BEAMINTERACTION::SearchStrategy::bounding_volume_hierarchy)
  {
    // Get vector of all beam element bounding boxes.
    int const numroweles = ele_type_map_extractor_ptr()->beam_map()->NumMyElements();
    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> beam_bounding_boxes;
    for (int rowele_i = 0; rowele_i < numroweles; ++rowele_i)
    {
      int const elegid = ele_type_map_extractor_ptr()->beam_map()->GID(rowele_i);
      Core::Elements::Element* currele = Discret().gElement(elegid);

      beam_bounding_boxes.emplace_back(std::make_pair(elegid,
          currele->GetBoundingVolume(Discret(), *beam_interaction_data_state_ptr()->GetDisColNp(),
              *geometric_search_params_ptr_)));
    }

    // Get vector of the bounding boxes of all possible interacting elements (also including beams
    // if beam-to-beam contact is activated).
    int const numcoleles = Discret().NumMyColElements();
    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> other_bounding_boxes;
    for (int colele_i = 0; colele_i < numcoleles; ++colele_i)
    {
      // Check if the current element is relevant for beam-to-xxx contact.
      Core::Elements::Element* currele = Discret().lColElement(colele_i);
      const Core::Binstrategy::Utils::BinContentType contact_type =
          BEAMINTERACTION::UTILS::ConvertElementToBinContentType(currele);
      if (std::find(contactelementtypes_.begin(), contactelementtypes_.end(), contact_type) !=
          contactelementtypes_.end())
      {
        other_bounding_boxes.emplace_back(std::make_pair(currele->Id(),
            currele->GetBoundingVolume(Discret(), *beam_interaction_data_state_ptr()->GetDisColNp(),
                *geometric_search_params_ptr_.getConst())));
      }
    }

    // Get colliding pairs.
    const auto& [indices, offsets] = CollisionSearch(other_bounding_boxes, beam_bounding_boxes,
        Discret().Comm(), geometric_search_params_ptr_->verbosity_);

    // Create the beam-to-xxx pair pointers according to the search.
    for (size_t i_beam = 0; i_beam < beam_bounding_boxes.size(); i_beam++)
    {
      const int beam_gid = beam_bounding_boxes[i_beam].first;
      Core::Elements::Element* currele = Discret().gElement(beam_gid);
      std::set<Core::Elements::Element*> neighboring_elements;
      for (int j = offsets[i_beam]; j < offsets[i_beam + 1]; j++)
      {
        neighboring_elements.insert(Discret().gElement(other_bounding_boxes[indices[j]].first));
      }

      // sort out elements that should not be considered in contact evaluation
      select_eles_to_be_considered_for_contact_evaluation(currele, neighboring_elements);

      nearby_elements_map_[beam_gid] = neighboring_elements;
    }

    // Check if the primitives and predicates should be output
    if (geometric_search_visualization_ptr_ != Teuchos::null)
    {
      // Output is desired, so create it right here, because we only search the pairs once per time
      // step anyways.
      geometric_search_visualization_ptr_->write_primitives_and_predicates_to_disk(
          GState().get_time_n(), GState().get_step_n(), other_bounding_boxes, beam_bounding_boxes);
    }
  }
  else
    FOUR_C_THROW("No appropriate search strategy for beam interaction chosen!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::
    select_eles_to_be_considered_for_contact_evaluation(
        Core::Elements::Element* currele, std::set<Core::Elements::Element*>& neighbors) const
{
  check_init();

  // sort out elements that should not be considered in contact evaluation
  std::set<Core::Elements::Element*>::iterator eiter;
  for (eiter = neighbors.begin(); eiter != neighbors.end();)
  {
    bool toerase = false;
    // 1) ensure that an element will not be in contact with it self
    if (currele->Id() == (*eiter)->Id())
    {
      toerase = true;
    }
    // 2) ensure each contact only evaluated once (keep in mind that we are
    //    using FEMatrices and FEvectors -> || (*eiter)->Owner() != myrank not necessary)
    // note: as we are only looping over beam elements, only beam to beam contact needs id check
    // here
    else if (dynamic_cast<Discret::ELEMENTS::Beam3Base*>(*eiter) != nullptr and
             not(currele->Id() < (*eiter)->Id()))
    {
      toerase = true;
    }
    // 3) ensure that two elements sharing the same node do not get into contact
    else
    {
      for (int i = 0; i < (*eiter)->num_node(); ++i)
        for (int j = 0; j < currele->num_node(); ++j)
          if ((*eiter)->NodeIds()[i] == currele->NodeIds()[j]) toerase = true;
    }

    if (toerase)
      neighbors.erase(eiter++);
    else
      ++eiter;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::create_beam_contact_element_pairs()
{
  // Todo maybe keep existing pairs and reuse them ?
  contact_elepairs_.clear();
  assembly_managers_.clear();

  // clear the geometry evaluation data
  beam_interaction_conditions_ptr_->clear();

  std::map<int, std::set<Core::Elements::Element*>>::const_iterator nearbyeleiter;

  for (nearbyeleiter = nearby_elements_map_.begin(); nearbyeleiter != nearby_elements_map_.end();
       ++nearbyeleiter)
  {
    const int elegid = nearbyeleiter->first;
    std::vector<Core::Elements::Element const*> ele_ptrs(2);
    ele_ptrs[0] = DiscretPtr()->gElement(elegid);

#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (dynamic_cast<Discret::ELEMENTS::Beam3Base const*>(ele_ptrs[0]) == nullptr)
      FOUR_C_THROW("first element of element pair must be a beam element");
#endif

    std::set<Core::Elements::Element*>::const_iterator secondeleiter;
    for (secondeleiter = nearbyeleiter->second.begin();
         secondeleiter != nearbyeleiter->second.end(); ++secondeleiter)
    {
      ele_ptrs[1] = *secondeleiter;

      // construct, init and setup contact pairs
      Teuchos::RCP<BEAMINTERACTION::BeamContactPair> newbeaminteractionpair =
          BEAMINTERACTION::BeamContactPair::Create(ele_ptrs, beam_interaction_conditions_ptr_);

      if (newbeaminteractionpair != Teuchos::null)
      {
        newbeaminteractionpair->init(beam_contact_params_ptr_, ele_ptrs);
        newbeaminteractionpair->setup();

        // add to list of current contact pairs
        contact_elepairs_.push_back(newbeaminteractionpair);
      }
    }
  }

  // Setup the geometry evaluation data.
  beam_interaction_conditions_ptr_->setup(DiscretPtr());

  // Get the pairs that can be assembled directly.
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_pairs_direct;
  for (auto& elepairptr : contact_elepairs_)
    if (elepairptr->IsAssemblyDirect()) assembly_pairs_direct.push_back(elepairptr);

  // Check if there are any processors that require a direct element assembly method.
  // We need to do this as in some assembly methods MPI communications are needed and the
  // simulation crashes if the assembly manager is not on all ranks.
  int my_direct_pairs = assembly_pairs_direct.size();
  int global_direct_pairs = 0;
  Discret().Comm().SumAll(&my_direct_pairs, &global_direct_pairs, 1);

  // Create the needed assembly manager.
  if (global_direct_pairs > 0)
    assembly_managers_.push_back(
        Teuchos::rcp<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect>(
            new BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect(
                assembly_pairs_direct)));

  // Each indirect assembly manager depends on a beam interaction.
  beam_interaction_conditions_ptr_->create_indirect_assembly_managers(
      DiscretPtr(), assembly_managers_);

  Core::IO::cout(Core::IO::standard)
      << "PID " << std::setw(2) << std::right << GState().get_my_rank() << " currently monitors "
      << std::setw(5) << std::right << contact_elepairs_.size() << " beam contact pairs"
      << Core::IO::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::set_restart_displacement_in_pairs()
{
  check_init_setup();

  if (beam_interaction_data_state_ptr()->get_restart_coupling_flag())
  {
    for (auto& pair : contact_elepairs_)
    {
      // Element Dof values at the restart state.
      std::vector<std::vector<double>> element_restart_dispalcement_(2);

      for (unsigned int i_element = 0; i_element < 2; ++i_element)
      {
        // Extract the Dof values of this element from the restart vector
        BEAMINTERACTION::UTILS::ExtractPosDofVecValues(Discret(), pair->GetElement(i_element),
            beam_interaction_data_state_ptr()->GetDisRestartCol(),
            element_restart_dispalcement_[i_element]);
      }

      // Set the displacement state in the pair.
      pair->set_restart_displacement(element_restart_dispalcement_);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::print_all_beam_contact_element_pairs(
    std::ostream& out) const
{
  out << "\n\nCurrent BeamContactElementPairs: ";
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>::const_iterator iter;
  for (iter = contact_elepairs_.begin(); iter != contact_elepairs_.end(); ++iter)
    (*iter)->print(out);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::print_active_beam_contact_set(
    std::ostream& out) const
{
  bool atleastoneactivepair = false;

  for (auto& elepairptr : contact_elepairs_)
    if (elepairptr->GetContactFlag() == true) atleastoneactivepair = true;


  if (atleastoneactivepair)
  {
    out << "\n    Active Beam-To-? Contact Set (PID " << GState().get_my_rank()
        << "):-----------------------------------------\n";
    out << "    ID1            ID2              T    xi       eta      angle    gap         "
           "force\n";


    for (auto& elepairptr : contact_elepairs_)
      elepairptr->print_summary_one_line_per_active_segment_pair(out);

    out << std::endl;
  }
}

FOUR_C_NAMESPACE_CLOSE
