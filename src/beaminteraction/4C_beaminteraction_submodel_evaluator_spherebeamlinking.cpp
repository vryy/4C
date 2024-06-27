/*-----------------------------------------------------------*/
/*! \file

\brief class for managing rigid sphere to beam crosslinking


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_beaminteraction_submodel_evaluator_spherebeamlinking.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_crosslinker_handler.hpp"
#include "4C_beaminteraction_crosslinker_node.hpp"
#include "4C_beaminteraction_link_beam3_reissner_line2_pinjointed.hpp"
#include "4C_beaminteraction_link_pinjointed.hpp"
#include "4C_beaminteraction_spherebeamlinking_params.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_beaminteraction_submodel_evaluator_crosslinking.hpp"
#include "4C_fem_geometry_periodic_boundingbox.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_rigidsphere.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtp_output.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::SphereBeamLinking()
    : sm_crosslinkink_ptr_(Teuchos::null),
      spherebeamlinking_params_ptr_(Teuchos::null),
      visualization_manager_ptr_(Teuchos::null),
      random_number_sphere_beam_linking_step_(-1)
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::setup()
{
  check_init();

  // construct, init and setup data container for crosslinking
  spherebeamlinking_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::SphereBeamLinkingParams());
  spherebeamlinking_params_ptr_->init(GState());
  spherebeamlinking_params_ptr_->setup();

  random_number_sphere_beam_linking_step_ = -1;

  // this includes temporary change in ghosting
  BEAMINTERACTION::UTILS::SetFilamentBindingSpotPositions(
      DiscretPtr(), spherebeamlinking_params_ptr_);

  // build runtime visualization output writer
  if (GInOutput().get_runtime_vtp_output_params() != Teuchos::null) init_output_runtime();

  // set flag
  issetup_ = true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::post_setup()
{
  check_init_setup();
  // nothing to do (yet)
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::init_submodel_dependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap)
{
  check_init_setup();

  // init pointer to crosslinker submodel
  STR::MODELEVALUATOR::BeamInteraction::Map::const_iterator miter;
  for (miter = (*submodelmap).begin(); miter != (*submodelmap).end(); ++miter)
    if (miter->first == Inpar::BEAMINTERACTION::submodel_crosslinking)
      sm_crosslinkink_ptr_ =
          Teuchos::rcp_dynamic_cast<BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking>(
              miter->second);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::reset()
{
  check_init_setup();

  // reset crosslinker pairs
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      // get elements
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (sphere != dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(
                        Discret().gElement(elepairptr->GetEleGid(0))))
        FOUR_C_THROW(" Rigid Sphere element has stored wrong linker. ");
#endif

      Discret::ELEMENTS::Beam3Base const* beamele =
          dynamic_cast<Discret::ELEMENTS::Beam3Base const*>(
              Discret().gElement(elepairptr->GetEleGid(1)));

      // init position of linker nodes
      std::vector<Core::LinAlg::Matrix<3, 1>> pos(2, Core::LinAlg::Matrix<3, 1>(true));

      // sphere current position
      std::vector<double> sphereeledisp;
      BEAMINTERACTION::UTILS::GetCurrentElementDis(
          Discret(), sphere, beam_interaction_data_state().GetDisColNp(), sphereeledisp);

      // note: sphere has just one node (with three translational dofs)
      for (unsigned int dim = 0; dim < 3; ++dim)
        pos[0](dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

      // beam bspot pos
      std::vector<double> beameledisp;
      BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele,
          beam_interaction_data_state().GetDisColNp(), PeriodicBoundingBox(), beameledisp);
      beamele->GetPosOfBindingSpot(pos[1], beameledisp,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(),
          elepairptr->GetLocBSpotNum(1), PeriodicBoundingBox());

      // unshift one of the positions if both are separated by a periodic boundary
      // condition, i.e. have been shifted before
      periodic_bounding_box_ptr()->UnShift3D(pos[1], pos[0]);

      // dummy triad
      std::vector<Core::LinAlg::Matrix<3, 3>> dummy_triad(2, Core::LinAlg::Matrix<3, 3>(true));

      // finally reset state
      elepairptr->ResetState(pos, dummy_triad);
    }
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::evaluate_force()
{
  check_init_setup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector<Core::LinAlg::SerialDenseVector> bspotforce(2, Core::LinAlg::SerialDenseVector(6));

  // resulting discrete element force vectors of the two parent elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      for (unsigned int i = 0; i < 2; ++i)
      {
        elegids[i] = elepairptr->GetEleGid(i);
        bspotforce[i].putScalar(0.0);
      }

      // evaluate beam linkage object to get forces of binding spots
      elepairptr->evaluate_force(bspotforce[0], bspotforce[1]);

      // apply forces on binding spots to parent elements
      // and get their discrete element force vectors
      BEAMINTERACTION::UTILS::ApplyBindingSpotForceToParentElements<Discret::ELEMENTS::Rigidsphere,
          Discret::ELEMENTS::Beam3Base>(Discret(), periodic_bounding_box_ptr(),
          beam_interaction_data_state_ptr()->GetDisColNp(), elepairptr, bspotforce, eleforce);

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BEAMINTERACTION::UTILS::fe_assemble_ele_force_stiff_into_system_vector_matrix(Discret(),
          elegids, eleforce, dummystiff, beam_interaction_data_state_ptr()->GetForceNp(),
          Teuchos::null);
    }
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::evaluate_stiff()
{
  check_init_setup();

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> bspotstiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2, Core::LinAlg::SerialDenseMatrix(6, 6)));

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  std::vector<Core::LinAlg::SerialDenseVector> dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      for (unsigned int i = 0; i < 2; ++i)
      {
        elegids[i] = elepairptr->GetEleGid(i);

        for (unsigned int j = 0; j < 2; ++j) bspotstiff[i][j].putScalar(0.0);
      }

      // evaluate beam linkage object to get linearizations of forces on binding spots
      elepairptr->evaluate_stiff(
          bspotstiff[0][0], bspotstiff[0][1], bspotstiff[1][0], bspotstiff[1][1]);

      // apply linearizations to parent elements and get their discrete element stiffness matrices
      BEAMINTERACTION::UTILS::ApplyBindingSpotStiffToParentElements<Discret::ELEMENTS::Rigidsphere,
          Discret::ELEMENTS::Beam3Base>(Discret(), periodic_bounding_box_ptr(),
          beam_interaction_data_state_ptr()->GetDisColNp(), elepairptr, bspotstiff, elestiff);

      // assemble the contributions into stiffness matrix class variable
      // stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
      BEAMINTERACTION::UTILS::fe_assemble_ele_force_stiff_into_system_vector_matrix(Discret(),
          elegids, dummyforce, elestiff, Teuchos::null,
          beam_interaction_data_state_ptr()->GetStiff());
    }
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::evaluate_force_stiff()
{
  check_init_setup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector<Core::LinAlg::SerialDenseVector> bspotforce(2, Core::LinAlg::SerialDenseVector(6));

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> bspotstiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2, Core::LinAlg::SerialDenseMatrix(6, 6)));

  // resulting discrete element force vectors of the two parent elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      for (unsigned int i = 0; i < 2; ++i)
      {
        elegids[i] = elepairptr->GetEleGid(i);
        bspotforce[i].putScalar(0.0);

        for (int j = 0; j < 2; ++j) bspotstiff[i][j].putScalar(0.0);
      }

      // evaluate beam linkage object to get forces on binding spots
      elepairptr->evaluate_force_stiff(bspotforce[0], bspotforce[1], bspotstiff[0][0],
          bspotstiff[0][1], bspotstiff[1][0], bspotstiff[1][1]);

      // apply forces on binding spots and corresponding linearizations to parent elements
      // and get their discrete element force vectors and stiffness matrices
      BEAMINTERACTION::UTILS::ApplyBindingSpotForceStiffToParentElements<
          Discret::ELEMENTS::Rigidsphere, Discret::ELEMENTS::Beam3Base>(Discret(),
          periodic_bounding_box_ptr(), beam_interaction_data_state_ptr()->GetDisColNp(), elepairptr,
          bspotforce, bspotstiff, eleforce, elestiff);

      // assemble the contributions into force and stiffness class variables
      // f_crosslink_np_ptr_, stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
      BEAMINTERACTION::UTILS::fe_assemble_ele_force_stiff_into_system_vector_matrix(Discret(),
          elegids, eleforce, elestiff, beam_interaction_data_state_ptr()->GetForceNp(),
          beam_interaction_data_state_ptr()->GetStiff());
    }
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateStepState(const double& timefac_n)
{
  check_init_setup();
}
/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::pre_update_step_element(
    bool beam_redist)
{
  check_init_setup();
  // not repartition of binning discretization necessary
  return false;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateStepElement(
    bool repartition_was_done)
{
  check_init_setup();

  // some screen output ( 0 = num_linker, 1 = num_newlinker, 2 = num_dissolved)
  std::vector<int> num_local(3, 0);
  std::vector<int> num_global(3, 0);

  // consider new bonds
  std::map<int, std::vector<std::pair<int, int>>> newlinks;
  newlinks.clear();
  if (spherebeamlinking_params_ptr_->MaxNumLinkerPerType()[0] > 0)
  {
    find_and_store_neighboring_elements(newlinks);
    create_beam_to_sphere_joint(newlinks);
  }

  for (auto const& iter : newlinks) num_local[1] += static_cast<int>(iter.second.size());

  // consider possible unbinding
  unbind_sphere_beam_bonds(num_local[2]);

  // sphere loop
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere*>(Discret().gElement(elegid));

    num_local[0] += sphere->GetNumberOfBonds();
  }

  // consider sphere linker contraction
  update_linker_length();

  // build sum over all procs
  MPI_Reduce(num_local.data(), num_global.data(), 3, MPI_INT, MPI_SUM, 0,
      dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());

  if (GState().get_my_rank() == 0)
  {
    Core::IO::cout(Core::IO::standard)
        << "\n************************************************" << Core::IO::endl;
    Core::IO::cout(Core::IO::standard) << "Sphere Beam Links: " << num_global[0];
    Core::IO::cout(Core::IO::standard) << " (New: " << num_global[1];
    Core::IO::cout(Core::IO::standard) << " Dissolved: " << num_global[2];
    Core::IO::cout(Core::IO::standard) << ")\n************************************************\n"
                                       << Core::IO::endl;
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::post_update_step_element()
{
  check_init_setup();

  //  std::set<int> spherebingids;
  //  int spheregid = -1;
  //  for( int i = 0; i < ele_type_map_extractor_ptr()->sphere_map()->NumMyElements(); ++i )
  //  {
  //    spheregid =  ele_type_map_extractor_ptr()->sphere_map()->GID(i);
  //    spherebingids.insert(
  //    beam_interaction_data_state_ptr()->GetRowEleToBinSet(spheregid).begin(),
  //                          beam_interaction_data_state_ptr()->GetRowEleToBinSet(spheregid).end()
  //                          );
  //  }
  //
  //  sm_crosslinkink_ptr->unbind_crosslinker_in_bins_and_neighborhood( spherebingids, false );

  //  int const updateevery = 200;
  //
  //  if( (global_state().get_step_n() + 1) % updateevery == 0 and global_state().get_step_n() >
  //  updateevery)
  //    sm_crosslinkink_ptr->double_bind_crosslinker_in_bins_and_neighborhood( spherebingids );
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
std::map<STR::EnergyType, double>
BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::get_energy() const
{
  check_init_setup();

  std::map<STR::EnergyType, double> sp_beam_link_energies;

  int unsigned const numrowsphereeles = EleTypeMapExtractor().sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = EleTypeMapExtractor().sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    for (auto const& bond_iter : sphere->GetBondMap())
    {
      sp_beam_link_energies[STR::beam_to_sphere_link_internal_energy] +=
          bond_iter.second->GetInternalEnergy();
      sp_beam_link_energies[STR::beam_to_sphere_link_kinetic_energy] +=
          bond_iter.second->GetKineticEnergy();
    }
  }

  return sp_beam_link_energies;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::OutputStepState(
    Core::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::runtime_output_step_state() const
{
  check_init_setup();

  if (visualization_manager_ptr_ != Teuchos::null) write_output_runtime();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::ResetStepState()
{
  check_init_setup();

  // in case time step is same as structure time step, update it
  spherebeamlinking_params_ptr_->ResetTimeStep((*GState().get_delta_time())[0]);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::write_restart(
    Core::IO::DiscretizationWriter& ia_writer, Core::IO::DiscretizationWriter& bin_writer) const
{
  check_init_setup();

  // as bonds are stored in rigid sphere element, nothing to do here
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PreReadRestart()
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::read_restart(
    Core::IO::DiscretizationReader& ia_reader, Core::IO::DiscretizationReader& bin_reader)
{
  check_init_setup();

  // as bonds are stored in rigid sphere element, nothing to do here
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PostReadRestart()
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::AddBinsToBinColMap(
    std::set<int>& colbins)
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::
    add_bins_with_relevant_content_for_ia_discret_col_map(std::set<int>& colbins) const
{
  check_init_setup();
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::get_half_interaction_distance(
    double& half_interaction_distance)
{
  double spherebeamlinking_half_interaction_distance = 0.0;

  // loop over all sphere elements (needed in case interaction distance should be
  // radius dependent in the future)
  double curr_ia_dist = 0.0;
  int unsigned const numrowsphereeles = EleTypeMapExtractor().sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    //    int const elegid = EleTypeMapExtractor().sphere_map()->GID(rowele_i);
    //    Discret::ELEMENTS::Rigidsphere * sphere =
    //        dynamic_cast< Discret::ELEMENTS::Rigidsphere * >( discret().gElement(elegid) );

    curr_ia_dist =
        0.5 * spherebeamlinking_params_ptr_->GetLinkerMaterial()->linking_length_tolerance();

    // update distance
    spherebeamlinking_half_interaction_distance =
        (curr_ia_dist > spherebeamlinking_half_interaction_distance)
            ? curr_ia_dist
            : spherebeamlinking_half_interaction_distance;
  }

  // get global maximum
  double spherebeamlinking_half_interaction_distance_global = 0.0;
  // build sum over all procs
  MPI_Allreduce(&spherebeamlinking_half_interaction_distance,
      &spherebeamlinking_half_interaction_distance_global, 1, MPI_DOUBLE, MPI_MAX,
      dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());

  // some screen output
  if (GState().get_my_rank() == 0)
    Core::IO::cout(Core::IO::verbose)
        << "\n spherebeamlinking half interaction distance "
        << spherebeamlinking_half_interaction_distance_global << Core::IO::endl;

  half_interaction_distance =
      (spherebeamlinking_half_interaction_distance_global > half_interaction_distance)
          ? spherebeamlinking_half_interaction_distance_global
          : half_interaction_distance;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::init_output_runtime()
{
  check_init();

  visualization_manager_ptr_ = Teuchos::rcp(new Core::IO::VisualizationManager(
      Core::IO::VisualizationParametersFactory(
          Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::Instance()->OutputControlFile(), GState().get_time_n()),
      BinDiscretPtr()->Comm(), "spherebeamlinker"));
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::write_output_runtime() const
{
  check_init_setup();

  // get number of linker on current proc
  unsigned int num_row_points = 0;
  int unsigned const numrowsphereeles = EleTypeMapExtractor().sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = EleTypeMapExtractor().sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    num_row_points += sphere->GetNumberOfBonds();
  }

  // set geometry manually
  const unsigned int num_spatial_dimensions = 3;

  // get and prepare storage for point coordinate values
  std::vector<double>& point_coordinates =
      visualization_manager_ptr_->get_visualization_data().GetPointCoordinates();
  point_coordinates.clear();
  point_coordinates.reserve(num_spatial_dimensions * num_row_points);

  // init output desired output vectors
  std::vector<double> currlength(num_row_points, 0.0);
  std::vector<double> orientation(num_row_points * num_spatial_dimensions, 0.0);
  std::vector<double> force(num_row_points * num_spatial_dimensions, 0.0);
  Core::LinAlg::SerialDenseVector bspotforce(num_spatial_dimensions);

  // set position of spherebeamlinks
  unsigned int bond_i = 0;
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = EleTypeMapExtractor().sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere const* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      Discret::ELEMENTS::Beam3Base const* beamele =
          dynamic_cast<Discret::ELEMENTS::Beam3Base const*>(
              Discret().gElement(elepairptr->GetEleGid(1)));

      // init position of linker nodes
      std::vector<Core::LinAlg::Matrix<3, 1>> pos(2, Core::LinAlg::Matrix<3, 1>(true));

      // sphere current position
      std::vector<double> sphereeledisp;
      BEAMINTERACTION::UTILS::GetCurrentElementDis(
          Discret(), sphere, beam_interaction_data_state().GetDisColNp(), sphereeledisp);

      // note: sphere has just one node (with three translational dofs)
      for (unsigned int dim = 0; dim < 3; ++dim)
        pos[0](dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

      // beam bspot pos
      std::vector<double> beameledisp;
      BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele,
          beam_interaction_data_state().GetDisColNp(), PeriodicBoundingBox(), beameledisp);
      beamele->GetPosOfBindingSpot(pos[1], beameledisp,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(),
          elepairptr->GetLocBSpotNum(1), PeriodicBoundingBox());

      // unshift one of the positions if both are separated by a periodic boundary
      // condition, i.e. have been shifted before
      PeriodicBoundingBox().UnShift3D(pos[1], pos[0]);

      elepairptr->GetBindingSpotForce(0, bspotforce);
      // set point coordinate value
      for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
      {
        point_coordinates.push_back(pos[0](idim) + 0.5 * (pos[1](idim) - pos[0](idim)));
        orientation[bond_i * num_spatial_dimensions + idim] = (pos[1](idim) - pos[0](idim));
        force[bond_i * num_spatial_dimensions + idim] = bspotforce(idim);
      }

      ++bond_i;
    }
  }

  // append all desired output data to the writer object's storage
  // i) number of bonds
  visualization_manager_ptr_->get_visualization_data().SetPointDataVector(
      "orientation", orientation, 3);
  // ii) linker force
  visualization_manager_ptr_->get_visualization_data().SetPointDataVector("force", force, 3);

  // finalize everything and write all required VTU files to filesystem
  visualization_manager_ptr_->WriteToDisk(GState().get_time_n(), GState().get_step_n());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::find_and_store_neighboring_elements(
    std::map<int, std::vector<std::pair<int, int>>>& newlinks)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::"
      "SphereBeamLinking::find_and_store_neighboring_elements");

  check_init_setup();

  std::unordered_set<int> tobebonded;

  // loop over all sphere elements
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  std::vector<int> rand_row_sphere = BEAMINTERACTION::UTILS::Permutation(numrowsphereeles);

  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rand_row_sphere[rowele_i]);
    Core::Elements::Element* currsphere = DiscretPtr()->gElement(elegid);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all bins touched by currele
    for (auto const& biniter : beam_interaction_data_state_ptr()->GetRowEleToBinSet(elegid))
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->get_neighbor_and_own_bin_ids(biniter, loc_neighboring_binIds);

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert(loc_neighboring_binIds.begin(), loc_neighboring_binIds.end());
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds(neighboring_binIds.begin(), neighboring_binIds.end());

    // set of beam elements that reside in neighboring bins
    std::set<Core::Elements::Element*> neighboring_elements;
    std::vector<Core::Binstrategy::Utils::BinContentType> bc(
        1, Core::Binstrategy::Utils::BinContentType::Beam);
    BinStrategyPtr()->GetBinContent(neighboring_elements, bc, glob_neighboring_binIds);

    // -------------------------------------------------------------------------
    // NOTE: This is crucial for reproducibility to ensure that computation does
    // not depend on pointer addresses (see also comment of class Less)
    // -------------------------------------------------------------------------
    std::vector<Core::Elements::Element const*> beamvec(
        neighboring_elements.begin(), neighboring_elements.end());
    std::sort(beamvec.begin(), beamvec.end(), BEAMINTERACTION::UTILS::Less());

    // sort out elements that do not meet bind event criteria
    check_feasibility_of_new_link(currsphere, beamvec, tobebonded, newlinks);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::check_feasibility_of_new_link(
    Core::Elements::Element const* currele,
    std::vector<Core::Elements::Element const*> const& neighbors,
    std::unordered_set<int>& tobebonded,
    std::map<int, std::vector<std::pair<int, int>>>& newlinks) const
{
  check_init_setup();

  int numnewbondsthisstep = 0;

  Discret::ELEMENTS::Rigidsphere const* sphere =
      dynamic_cast<Discret::ELEMENTS::Rigidsphere const*>(currele);

  // current sphere position
  std::vector<double> sphereeledisp;
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      Discret(), currele, beam_interaction_data_state().GetDisColNp(), sphereeledisp);

  // note: sphere has just one node (with three translational dofs)
  Core::LinAlg::Matrix<3, 1> spherepos(true);
  for (unsigned int dim = 0; dim < 3; ++dim)
    spherepos(dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

  // loop over all neighboring elements
  std::vector<int> rand_ele = BEAMINTERACTION::UTILS::Permutation(neighbors.size());
  for (auto const& eiter : rand_ele)
  {
    Discret::ELEMENTS::Beam3Base const* beamele =
        dynamic_cast<Discret::ELEMENTS::Beam3Base const*>(neighbors[eiter]);

#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (sphere == nullptr or beamele == nullptr)
      FOUR_C_THROW(" First element should be a sphere, second element a beam.");
#endif

    std::vector<double> beameledisp;
    BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele,
        beam_interaction_data_state().GetDisColNp(), PeriodicBoundingBox(), beameledisp);

    Core::LinAlg::Matrix<3, 1> bspotpos(true);
    Core::LinAlg::Matrix<3, 3> bspottriad(true);

    // loop over binding spots of neighboring element
    unsigned int numbspots = beamele->get_number_of_binding_spots(
        spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType());
    std::vector<int> rand_bsp = BEAMINTERACTION::UTILS::Permutation(numbspots);
    for (unsigned int bspot_i = 0; bspot_i < numbspots; ++bspot_i)
    {
      // build unique linker id from elegid and local binding spot id
      std::pair<int, int> bspotpair = std::make_pair(beamele->Id(), rand_bsp[bspot_i]);
      int bpspotgid = BEAMINTERACTION::UTILS::CantorPairing(bspotpair);

      // criterion: has sphere reached maximum number of admissible bonds?
      if ((sphere->GetNumberOfBonds() + numnewbondsthisstep) ==
          spherebeamlinking_params_ptr_->MaxNumLinkerPerType()[0])
        continue;

      // criterion: probability check for integrin collagen binding
      if (spherebeamlinking_params_ptr_->GetLinkerMaterial()->KOn() > -1e-08)
      {
        double plink = 1.0 - exp((-1.0) * spherebeamlinking_params_ptr_->DeltaTime() *
                                 spherebeamlinking_params_ptr_->GetLinkerMaterial()->KOn());
        if (Global::Problem::Instance()->Random()->Uni() > plink) continue;
      }

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // todo: do only in debug mode as soon tested enough
      // safety check
      if ((sphere->GetNumberOfBonds() + numnewbondsthisstep) >
          spherebeamlinking_params_ptr_->MaxNumLinkerPerType()[0])
        FOUR_C_THROW(" sphere has more bonds than allowed. Something went wrong.");
#endif

      // criterion: does identical bond already exist?
      if (sphere->DoesBondExist(bpspotgid)) continue;

      // criterion: is beam binding spot free (covered by criterion 1? Only in case of separate
      // linker for cell to beam and beam to beam fixme: maybe replace beam to beam in this case
      // search literature if binding spots are the same


      // criterion: distance
      // get current position at binding spot xi
      bspotpos.clear();
      beamele->GetPosOfBindingSpot(bspotpos, beameledisp,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(), rand_bsp[bspot_i],
          PeriodicBoundingBox());

      double linkdistmin =
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLength() -
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->linking_length_tolerance();

      // exclude links inside sphere
      linkdistmin = (linkdistmin > sphere->Radius()) ? linkdistmin : sphere->Radius();

      double linkdistmax =
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLength() +
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->linking_length_tolerance();


      if (BEAMINTERACTION::UTILS::IsDistanceOutOfRange(
              spherepos, bspotpos, linkdistmin, linkdistmax))
        continue;

      // criterion: orientation
      /* NOTE: this works slightly different to crosslinking of two beams: here we check the angle
       * between the beams first base vector and the direction vector from the spheres mid point to
       * the binding spot, i.e. LINKINGANGLE in the input file means something slightly different in
       * this case */

      // get current triad at binding spot xi
      bspottriad.clear();
      beamele->get_triad_of_binding_spot(bspottriad, beameledisp,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(), rand_bsp[bspot_i]);

      // note: we use first base vector instead of tangent vector here
      Core::LinAlg::Matrix<3, 1> curr_bindingspot_beam_tangent(true);
      for (unsigned int idim = 0; idim < 3; ++idim)
        curr_bindingspot_beam_tangent(idim) = bspottriad(idim, 0);

      Core::LinAlg::Matrix<3, 1> dist_vec(true);
      dist_vec.update(1.0, bspotpos, -1.0, spherepos);

      double const linkanglemin =
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngle() -
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->linking_angle_tolerance();
      double const linkanglemax =
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngle() +
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->linking_angle_tolerance();

      if (BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
              dist_vec, curr_bindingspot_beam_tangent, linkanglemin, linkanglemax))

        // criterion: check if bspot will already be bonded this step
        if (tobebonded.find(bpspotgid) != tobebonded.end()) continue;

      // update control variables
      tobebonded.insert(bpspotgid);
      ++numnewbondsthisstep;

      // add new link
      newlinks[sphere->Id()].push_back(bspotpair);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::create_beam_to_sphere_joint(
    std::map<int, std::vector<std::pair<int, int>>> const& newlinks)
{
  check_init_setup();


  // loop over cells
  for (auto const& newlinkiter : newlinks)
  {
    // init position of linker nodes
    std::vector<Core::LinAlg::Matrix<3, 1>> pos(2, Core::LinAlg::Matrix<3, 1>(true));

    int const spheregid = newlinkiter.first;
    // get elements
    Discret::ELEMENTS::Rigidsphere* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere*>(Discret().gElement(spheregid));

    std::vector<std::pair<int, int>> eleids(2);
    // todo: for now, sphere has one binding spot, so local binding spot id 0
    eleids[0] = std::make_pair(spheregid, 0);

    // sphere current position
    std::vector<double> sphereeledisp;
    BEAMINTERACTION::UTILS::GetCurrentElementDis(
        Discret(), sphere, beam_interaction_data_state().GetDisColNp(), sphereeledisp);

    // note: sphere has just one node (with three translational dofs)
    for (unsigned int dim = 0; dim < 3; ++dim)
      pos[0](dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

    // loop over all integrins that are about to be bonded
    for (auto const& bspotiter : newlinkiter.second)
    {
      int const beamgid = bspotiter.first;
      eleids[1] = std::make_pair(beamgid, bspotiter.second);

      // get neighboring element
      Discret::ELEMENTS::Beam3Base* beamele =
          dynamic_cast<Discret::ELEMENTS::Beam3Base*>(Discret().gElement(beamgid));

      // beam bspot pos
      std::vector<double> beameledisp;
      BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele,
          beam_interaction_data_state().GetDisColNp(), PeriodicBoundingBox(), beameledisp);
      beamele->GetPosOfBindingSpot(pos[1], beameledisp,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(), bspotiter.second,
          PeriodicBoundingBox());

      // create and initialize objects of beam-to-beam connections
      // Todo introduce enum for type of linkage (only linear Beam3r element possible so far)
      //      and introduce corresponding input parameter
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr =
          BEAMINTERACTION::BeamLinkPinJointed::Create(Inpar::BEAMINTERACTION::truss);

      // unique linker id is bspot elegid and locspot id paired
      int id = BEAMINTERACTION::UTILS::CantorPairing(eleids[1]);

      // dummy triad
      std::vector<Core::LinAlg::Matrix<3, 3>> dummy_triad(2, Core::LinAlg::Matrix<3, 3>(true));

      // finally initialize and setup object
      linkelepairptr->init(id, eleids, pos, dummy_triad,
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkerType(), GState().get_time_np());
      // material id
      linkelepairptr->setup(
          spherebeamlinking_params_ptr_->GetLinkerMaterial()->beam_elast_hyper_mat_num());

      // set on status on element level
      sphere->AddBond(id, linkelepairptr);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::unbind_sphere_beam_bonds(
    int& num_dissolved)
{
  check_init_setup();

  // safety check
  if (spherebeamlinking_params_ptr_->DeltaTime() < 1.0e-8)
    FOUR_C_THROW("You are about to divide by (almost) zero");

  // check if unbinding needs to be checked this problem time step
  int random_number_sphere_beam_linking_step =
      static_cast<int>((GState().get_time_np() - (*GState().get_delta_time())[0]) /
                           spherebeamlinking_params_ptr_->DeltaTime() +
                       1.0e-8);
  if (random_number_sphere_beam_linking_step == random_number_sphere_beam_linking_step_)
    return;
  else
    random_number_sphere_beam_linking_step_ = random_number_sphere_beam_linking_step;

  // in case off rate is zero
  if (std::abs(spherebeamlinking_params_ptr_->GetLinkerMaterial()->KOff()) < 1e-08) return;

  // init variables
  double p_unbind = 0.0;

  // sphere loop
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  std::vector<int> rand_row_sphere = BEAMINTERACTION::UTILS::Permutation(numrowsphereeles);
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rand_row_sphere[rowele_i]);
    Discret::ELEMENTS::Rigidsphere* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere*>(Discret().gElement(elegid));

    // loop over bonds of current sphere
    std::vector<int> to_dissolve;
    to_dissolve.reserve(sphere->GetBondMap().size());
    for (auto const& ipair : sphere->GetBondMap())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      // check if linker was set this time step
      if (elepairptr->GetTimeLinkWasSet() == GState().get_time_np()) continue;

      // consider catch-slip bond behavior of integrin linkers
      calc_force_dependent_catch_slip_bond_unbind_probability(elepairptr, p_unbind);

      // if probability criterion is not met, we are done here
      if (Global::Problem::Instance()->Random()->Uni() > p_unbind) continue;

      to_dissolve.push_back(ipair.first);
      ++num_dissolved;
    }

    // dissolve all bonds
    for (unsigned int i_td = 0; i_td < to_dissolve.size(); ++i_td)
      sphere->dissolve_bond(to_dissolve[i_td]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::
    calc_force_dependent_catch_slip_bond_unbind_probability(
        Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr, double& p_unbind)
{
  check_init_setup();

  // Note: this needs to come after a contraction was equilibrated by the network. Doing this
  // directly after changing the linker reference length does not make sense.

  // todo: maybe add these to linker material input line
  double const phi_FA_s = 7.78;
  double const phi_FA_c = 4.02;
  double const f_FA = 5.38;
  double const k_FA_0 = spherebeamlinking_params_ptr_->GetLinkerMaterial()->KOff();
  double dt = spherebeamlinking_params_ptr_->DeltaTime();

  // Fixme: is force 1 the correct one, do we need to take the force on the end that is connected to
  // the beam? note: as we only check unbinding in links that were set before the current time step,
  // we do not need to calculate the forces again.
  Core::LinAlg::SerialDenseVector bspotforce_one(6);
  linkelepairptr->GetBindingSpotForce(1, bspotforce_one);
  double f = Core::LinAlg::Norm2(bspotforce_one);


  // check if linker is stretched -> sgn+ or compressed -> sgn- by checking orientation of force
  // vector note: this works only if there are no other forces (like inertia, stochastic, damping)
  // acting on the linker
  Core::LinAlg::Matrix<3, 1> dist_vec(true);
  Core::LinAlg::Matrix<3, 1> bspotforceone(true);
  dist_vec.update(-1.0, linkelepairptr->GetBindSpotPos1(), 1.0, linkelepairptr->GetBindSpotPos2());
  for (unsigned int j = 0; j < 3; ++j) bspotforceone(j) = bspotforce_one(j);
  double sgn = (dist_vec.dot(bspotforceone) < 0.0) ? -1.0 : 1.0;

  /* alternative for linear centerline interpolation would be to compare
   * reference length and current length to see if element is stretched or compressed
   */

  // fixme: is that correct / does that make sense for compressive forces?
  f *= sgn;

  // calculate force dependent off rate for catch slip bond
  // ( see Wang et al. 2016 "Mechanosensitive subcellular rheostasis drives emergent
  //   single-cell mechanical homeostasis", nature materials. See supplementary information)
  double k_FA_off = k_FA_0 * (exp(f / f_FA - phi_FA_s) + exp(((-1.0) * f / f_FA) + phi_FA_c));

  // get respective force dependent unbind probability for each cl binding spot
  if (std::isfinite(k_FA_off))
  {
    p_unbind = 1.0 - exp((-1.0) * dt * k_FA_off);
  }
  else
  {
    p_unbind = 1.0;
    std::cout << "WARNING: You have very high forces (" << f
              << ") acting on your integrins. Are you "
                 "sure this is what you want? "
              << std::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::update_linker_length()
{
  check_init_setup();

  // adapt force/strain in linker
  // note: problem time step is used here
  double contraction_per_dt =
      spherebeamlinking_params_ptr_->ContractionRate(Inpar::BEAMINTERACTION::linkertype_integrin) *
      (*GState().get_delta_time())[0];
  double scalefac = 0.0;
  int unsigned const numrowsphereeles = ele_type_map_extractor_ptr()->sphere_map()->NumMyElements();
  for (unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->sphere_map()->GID(rowele_i);
    Discret::ELEMENTS::Rigidsphere* sphere =
        dynamic_cast<Discret::ELEMENTS::Rigidsphere*>(Discret().gElement(elegid));

    // todo: do we want this (note: this would be a compressible behavior in this case)
    // change radius of contracting cell
    //    if ( global_state().get_step_n() % 10 == 0 )
    //      sphere->ScaleRadius(0.95);

    // loop over bonds of current sphere
    for (auto const& ipair : sphere->GetBondMap())
    {
      // get pair object
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      // only contract if linker size > sphere radius * factor
      double factor = 1.01;
      if ((elepairptr->get_current_linker_length() <= sphere->Radius() * factor) or
          (GState().get_step_n() == 0))
        continue;

      // compute scaling factor for linker length
      scalefac = 1.0 - (contraction_per_dt / elepairptr->GetReferenceLength());

      // some safety checks
      if (scalefac < 0.0)
        FOUR_C_THROW(
            "Contraction %f of a linker of more than its reference length %f in one time step "
            "does not make sense.",
            contraction_per_dt, elepairptr->GetReferenceLength());
      if (contraction_per_dt > elepairptr->get_current_linker_length())
        FOUR_C_THROW("Contraction of %f for a linker with current length %f does not make sense.",
            contraction_per_dt, elepairptr->get_current_linker_length());

      // scale linker length / contract linker
      elepairptr->scale_linker_reference_length(scalefac);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
