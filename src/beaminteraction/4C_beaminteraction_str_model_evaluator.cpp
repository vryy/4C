/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation of all beam interaction terms, managing
       everything that has to do with parallelity


\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_str_model_evaluator.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_crosslinker_handler.hpp"
#include "4C_beaminteraction_crosslinker_node.hpp"
#include "4C_beaminteraction_data.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact.hpp"
#include "4C_beaminteraction_submodel_evaluator_crosslinking.hpp"
#include "4C_beaminteraction_submodel_evaluator_factory.hpp"
#include "4C_beaminteraction_submodel_evaluator_generic.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fem_geometry_periodic_boundingbox.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_beam_to_solid.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_rebalance_print.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteraction::BeamInteraction()
    : discret_ptr_(Teuchos::null),
      beaminteraction_params_ptr_(Teuchos::null),
      submodeltypes_(Teuchos::null),
      me_map_ptr_(Teuchos::null),
      me_vec_ptr_(Teuchos::null),
      myrank_(-1),
      coupsia_(Teuchos::null),
      siatransform_(Teuchos::null),
      ia_discret_(Teuchos::null),
      eletypeextractor_(Teuchos::null),
      ia_state_ptr_(Teuchos::null),
      ia_force_beaminteraction_(Teuchos::null),
      force_beaminteraction_(Teuchos::null),
      stiff_beaminteraction_(Teuchos::null),
      beam_crosslinker_handler_(Teuchos::null),
      binstrategy_(Teuchos::null),
      bindis_(Teuchos::null),
      rowbins_(Teuchos::null),
      dis_at_last_redistr_(Teuchos::null),
      half_interaction_distance_(-1.0)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::setup()
{
  check_init();

  // -------------------------------------------------------------------------
  // setup variables
  // -------------------------------------------------------------------------
  // discretization pointer
  discret_ptr_ = Teuchos::rcp_dynamic_cast<Core::FE::Discretization>(discret_ptr(), true);
  // stiff
  stiff_beaminteraction_ = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*global_state().dof_row_map_view(), 81, true, true));
  // force and displacement at last redistribution
  force_beaminteraction_ = Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map(), true));
  dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map(), true));
  // get myrank
  myrank_ = discret_ptr()->Comm().MyPID();

  beaminteraction_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::BeamInteractionParams());
  beaminteraction_params_ptr_->init();
  beaminteraction_params_ptr_->setup();

  // print logo
  Logo();

  // set submodel types
  set_sub_model_types();

  // -------------------------------------------------------------------------
  // clone problem discretization, the idea is simple: we redistribute only
  // the new discretization to enable all interactions (including the required
  // search), calculate the resulting force and stiffness contributions, export
  // them to our initial discretization where all evaluation, assembly and
  // solving is done. Therefore the maps of our initial discretization don't
  // change, i.e. there is no need to rebuild the global state.
  // -------------------------------------------------------------------------
  Teuchos::RCP<Core::FE::DiscretizationCreatorBase> discloner =
      Teuchos::rcp(new Core::FE::DiscretizationCreatorBase());
  ia_discret_ = discloner->create_matching_discretization(
      discret_ptr_, "ia_structure", true, true, false, true);
  // create discretization writer
  ia_discret_->SetWriter(Teuchos::rcp(new Core::IO::DiscretizationWriter(ia_discret_,
      Global::Problem::Instance()->OutputControlFile(),
      Global::Problem::Instance()->spatial_approximation_type())));

  // init data container
  ia_state_ptr_ = Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteractionDataState());
  ia_state_ptr_->init();
  ia_state_ptr_->setup(ia_discret_);

  ia_state_ptr_->GetDisNp() = Teuchos::rcp(new Epetra_Vector(*global_state_ptr()->get_dis_np()));
  BEAMINTERACTION::UTILS::PeriodicBoundaryConsistentDisVector(ia_state_ptr_->GetDisNp(),
      tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box(), ia_discret_);

  // -------------------------------------------------------------------------
  // initialize coupling adapter to transform matrices between the two discrets
  // (with distinct parallel distribution)
  // -------------------------------------------------------------------------
  coupsia_ = Teuchos::rcp(new Core::Adapter::Coupling());
  siatransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);

  // -------------------------------------------------------------------------
  // initialize and setup binning strategy and beam crosslinker handler
  // -------------------------------------------------------------------------
  // construct, init and setup binning strategy
  std::vector<Teuchos::RCP<Core::FE::Discretization>> discret_vec(1, ia_discret_);

  // We have to pass the displacement column vector to the initialization of the binning strategy.
  ia_state_ptr_->GetDisColNp() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  Core::LinAlg::Export(*ia_state_ptr_->GetDisNp(), *ia_state_ptr_->GetDisColNp());

  std::vector<Teuchos::RCP<const Epetra_Vector>> disp_vec(1, ia_state_ptr_->GetDisColNp());
  Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
  Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
      "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
      binning_params);
  binstrategy_ = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(binning_params,
      Global::Problem::Instance()->OutputControlFile(), ia_discret_->Comm(),
      ia_discret_->Comm().MyPID(), discret_vec, disp_vec));
  binstrategy_->set_deforming_binning_domain_handler(
      tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box());

  bindis_ = binstrategy_->BinDiscret();

  // construct, init and setup beam crosslinker handler and binning strategy
  // todo: move this and its single call during partition to crosslinker submodel
  if (HaveSubModelType(Inpar::BEAMINTERACTION::submodel_crosslinking))
  {
    beam_crosslinker_handler_ = Teuchos::rcp(new BEAMINTERACTION::BeamCrosslinkerHandler());
    beam_crosslinker_handler_->init(global_state().get_my_rank(), binstrategy_);
    beam_crosslinker_handler_->setup();
  }

  // some screen output for binning
  print_binning_info_to_screen();

  // extract map for each eletype that is in discretization
  eletypeextractor_ = Teuchos::rcp(new BEAMINTERACTION::UTILS::MapExtractor);
  BEAMINTERACTION::UTILS::SetupEleTypeMapExtractor(ia_discret_, eletypeextractor_);

  // initialize and setup submodel evaluators
  init_and_setup_sub_model_evaluators();

  // distribute problem according to bin distribution to procs ( in case of restart
  // partitioning is done during read_restart() )
  if (not Global::Problem::Instance()->restart()) partition_problem();

  // some actions need a partitioned system followed by a renewal of the partition
  if (not Global::Problem::Instance()->restart() and post_partition_problem()) partition_problem();

  post_setup();

  // some screen output
  Core::Rebalance::UTILS::print_parallel_distribution(*ia_discret_);
  Core::Rebalance::UTILS::print_parallel_distribution(*bindis_);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::post_setup()
{
  check_init();

  // post setup submodel loop
  Vector::iterator sme_iter;
  for (Vector::iterator sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->post_setup();

  if (beaminteraction_params_ptr_->get_repartition_strategy() ==
      Inpar::BEAMINTERACTION::repstr_adaptive)
  {
    // submodel loop to determine half interaction radius
    for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
      (*sme_iter)->get_half_interaction_distance(half_interaction_distance_);

    if (global_state().get_my_rank() == 0)
      std::cout << " half min bin size " << 0.5 * binstrategy_->GetMinBinSize() << std::endl;

    // safety checks
    if (half_interaction_distance_ > (0.5 * binstrategy_->GetMinBinSize()))
      FOUR_C_THROW(
          "Your half interaction distance %f is larger than half your smallest bin %f. You will "
          "not be\n"
          "able to track all interactions like this. Increase bin size?",
          half_interaction_distance_, 0.5 * binstrategy_->GetMinBinSize());
    if (half_interaction_distance_ < 0)
      FOUR_C_THROW("At least one model needs to define half interaction radius");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::set_sub_model_types()
{
  check_init();

  submodeltypes_ = Teuchos::rcp(new std::set<enum Inpar::BEAMINTERACTION::SubModelType>());

  // ---------------------------------------------------------------------------
  // check for crosslinking in biopolymer networks
  // ---------------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->beam_interaction_params().sublist("SPHERE BEAM LINK"),
          "SPHEREBEAMLINKING"))
    submodeltypes_->insert(Inpar::BEAMINTERACTION::submodel_spherebeamlink);

  // ---------------------------------------------------------------------------
  // check for crosslinking in biopolymer networks
  // ---------------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->beam_interaction_params().sublist("CROSSLINKING"),
          "CROSSLINKER"))
    submodeltypes_->insert(Inpar::BEAMINTERACTION::submodel_crosslinking);

  // ---------------------------------------------------------------------------
  // check for point to point penalty coupling conditions
  // ---------------------------------------------------------------------------

  // conditions for beam penalty point coupling
  std::vector<Core::Conditions::Condition*> beampenaltycouplingconditions(0);
  discret_ptr_->GetCondition("PenaltyPointCouplingCondition", beampenaltycouplingconditions);
  if (beampenaltycouplingconditions.size() > 0)
    submodeltypes_->insert(Inpar::BEAMINTERACTION::submodel_beamcontact);

  // ---------------------------------------------------------------------------
  // check for beam contact
  // ---------------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO BEAM CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO SPHERE CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID VOLUME MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE CONTACT"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none)
    submodeltypes_->insert(Inpar::BEAMINTERACTION::submodel_beamcontact);

  // ---------------------------------------------------------------------------
  // check for beam potential-based interactions
  // ---------------------------------------------------------------------------
  std::vector<Core::Conditions::Condition*> beampotconditions(0);
  discret().GetCondition("BeamPotentialLineCharge", beampotconditions);
  if (beampotconditions.size() > 0)
    submodeltypes_->insert(Inpar::BEAMINTERACTION::submodel_potential);

  // Check if all all combinations of submodel evaluators work
  if (Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO BEAM CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none and
      beampenaltycouplingconditions.size() > 0)
    FOUR_C_THROW(
        "It is not yet possible to use beam-to-beam contact in combination with beam-to-beam point "
        "coupling because every coupling point is also interpreted as a point of contact between 2 "
        "beams.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::init_and_setup_sub_model_evaluators()
{
  check_init();

  // model map
  me_map_ptr_ = BEAMINTERACTION::SUBMODELEVALUATOR::build_model_evaluators(*submodeltypes_);
  std::vector<enum Inpar::BEAMINTERACTION::SubModelType> sorted_submodeltypes(0);

  // build and sort submodel vector
  me_vec_ptr_ = transform_to_vector(*me_map_ptr_, sorted_submodeltypes);

  Vector::iterator sme_iter;
  for (sme_iter = (*me_vec_ptr_).begin(); sme_iter != (*me_vec_ptr_).end(); ++sme_iter)
  {
    (*sme_iter)->init(ia_discret_, bindis_, global_state_ptr(), global_in_output_ptr(),
        ia_state_ptr_, beam_crosslinker_handler_, binstrategy_,
        tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box(),
        Teuchos::rcp_dynamic_cast<BEAMINTERACTION::UTILS::MapExtractor>(eletypeextractor_, true));
    (*sme_iter)->setup();
  }

  // submodels build their pointer to other submodel objects to enable submodel dependencies
  // this is not particularly nice, at least the nicest way to handle such dependencies
  Vector::const_iterator iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->init_submodel_dependencies(me_map_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector>
STR::MODELEVALUATOR::BeamInteraction::transform_to_vector(
    STR::MODELEVALUATOR::BeamInteraction::Map submodel_map,
    std::vector<Inpar::BEAMINTERACTION::SubModelType>& sorted_submodel_types) const
{
  Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> me_vec_ptr =
      Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteraction::Vector(0));

  STR::MODELEVALUATOR::BeamInteraction::Map::iterator miter;

  // if there is a contractile cell submodel, put in first place
  miter = submodel_map.find(Inpar::BEAMINTERACTION::submodel_spherebeamlink);
  if (miter != submodel_map.end())
  {
    // put it in first place
    me_vec_ptr->push_back(miter->second);
    sorted_submodel_types.push_back(miter->first);
    submodel_map.erase(miter);
  }

  // insert the remaining submodel evaluators into the model vector
  for (miter = submodel_map.begin(); miter != submodel_map.end(); ++miter)
  {
    me_vec_ptr->push_back(miter->second);
    sorted_submodel_types.push_back(miter->first);
  }

  return me_vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::HaveSubModelType(
    Inpar::BEAMINTERACTION::SubModelType const& submodeltype) const
{
  check_init();
  return (submodeltypes_->find(submodeltype) != submodeltypes_->end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::partition_problem()
{
  check_init();

  // store structure discretization in vector
  std::vector<Teuchos::RCP<Core::FE::Discretization>> discret_vec(1, ia_discret_);

  // displacement vector according to periodic boundary conditions
  std::vector<Teuchos::RCP<Epetra_Vector>> mutabledisnp(
      1, Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap())));
  Core::LinAlg::Export(*ia_state_ptr_->GetDisNp(), *mutabledisnp[0]);

  std::vector<Teuchos::RCP<const Epetra_Vector>> disnp(
      1, Teuchos::rcp(new const Epetra_Vector(*mutabledisnp[0])));

  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int>>> nodesinbin(1);

  // weight for load balancing regarding the distribution of bins to procs
  // (this is experimental, choose what gives you best results)
  double const weight = 1.0;
  // get optimal row distribution of bins to procs
  rowbins_ =
      binstrategy_->weighted_distribution_of_bins_to_procs(discret_vec, disnp, nodesinbin, weight);

  // extract noderowmap because it will be called reset() after adding elements
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*bindis_->NodeRowMap()));
  // delete old bins ( in case you partition during your simulation or after a restart)
  bindis_->DeleteElements();
  binstrategy_->fill_bins_into_bin_discretization(rowbins_);

  // now node (=crosslinker) to bin (=element) relation needs to be
  // established in binning discretization. Therefore some nodes need to
  // change their owner according to the bins owner they reside in
  if (HaveSubModelType(Inpar::BEAMINTERACTION::submodel_crosslinking))
    beam_crosslinker_handler_->distribute_linker_to_bins(noderowmap);

  // determine boundary bins (physical boundary as well as boundary to other procs)
  binstrategy_->determine_boundary_row_bins();

  // determine one layer ghosting around boundary bins determined in previous step
  binstrategy_->determine_boundary_col_bins();

  // standard ghosting (if a proc owns a part of nodes (and therefore dofs) of
  // an element, the element and the rest of its nodes and dofs are ghosted
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmapdummy;
  binstrategy_->standard_discretization_ghosting(
      ia_discret_, rowbins_, ia_state_ptr_->GetDisNp(), stdelecolmap, stdnodecolmapdummy);

  // distribute elements that can be cut by the periodic boundary to bins
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  Core::LinAlg::Export(*ia_state_ptr_->GetDisNp(), *iadiscolnp);

  binstrategy_->distribute_row_elements_to_bins_using_ele_aabb(
      ia_discret_, ia_state_ptr_->GetBinToRowEleMap(), iadiscolnp);

  // build row elements to bin map
  build_row_ele_to_bin_map();

  // extend ghosting
  extend_ghosting();

  // assign Elements to bins
  binstrategy_->remove_all_eles_from_bins();
  binstrategy_->AssignElesToBins(ia_discret_, ia_state_ptr_->get_extended_bin_to_row_ele_map());

  // update maps of state vectors and matrices
  update_maps();

  // reset transformation
  update_coupling_adapter_and_matrix_transformation();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::post_partition_problem()
{
  check_init();

  bool repartition = false;

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    repartition = (*sme_iter)->post_partition_problem() ? true : repartition;

  return repartition;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::extend_ghosting()
{
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::BeamInteraction::extend_ghosting");

  ia_state_ptr_->get_extended_bin_to_row_ele_map().clear();

  check_init();

  std::set<int> colbins;

  // first, add default one layer ghosting
  std::set<int> const& boundcolbins = binstrategy_->BoundaryColBinsIds();
  colbins.insert(boundcolbins.begin(), boundcolbins.end());

  // if the bounding box of a row element of myrank touches a boundary col bin, we need
  // to ghost its neighborhood as well
  std::map<int, std::set<int>>::const_iterator it;
  std::vector<int> binvec(27);
  for (it = ia_state_ptr_->GetBinToRowEleMap().begin();
       it != ia_state_ptr_->GetBinToRowEleMap().end(); ++it)
  {
    // not doing the following if is only valid if you ensure that the largest element
    // in the discretization (in deformed state) is smaller than the smalles bin size
    // which is not necessarily needed e.g. for beam contact
    //    if( boundcolbins.find( it->first ) != boundcolbins.end() )
    {
      binstrategy_->get_neighbor_and_own_bin_ids(it->first, binvec);
      colbins.insert(binvec.begin(), binvec.end());
      binvec.clear();
    }
  }

  // enable submodel specific ghosting contributions to bin col map
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->AddBinsToBinColMap(colbins);

  // 1) extend ghosting of bin discretization
  // todo: think about if you really need to assign degrees of freedom for crosslinker
  // (now only needed in case you want to write output)
  binstrategy_->extend_ghosting_of_binning_discretization(rowbins_, colbins, true);

  // add submodel specific bins whose content should be ghosted in problem discret
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->add_bins_with_relevant_content_for_ia_discret_col_map(colbins);

  // build auxiliary bin col map
  std::vector<int> auxgids(colbins.begin(), colbins.end());
  Teuchos::RCP<Epetra_Map> auxmap = Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(auxgids.size()), auxgids.data(), 0, bindis_->Comm()));

  Teuchos::RCP<Epetra_Map> ia_elecolmap = binstrategy_->ExtendElementColMap(
      ia_state_ptr_->GetBinToRowEleMap(), ia_state_ptr_->GetBinToRowEleMap(),
      ia_state_ptr_->get_extended_bin_to_row_ele_map(), auxmap);

  // 2) extend ghosting of discretization
  BINSTRATEGY::UTILS::ExtendDiscretizationGhosting(ia_discret_, ia_elecolmap, true, false, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::reset(const Epetra_Vector& x)
{
  check_init_setup();

  // todo: somewhat illegal as of const correctness
  tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box()->ApplyDirichlet(
      global_state().get_time_n(), Global::Problem::Instance()->FunctionManager());

  // get current displacement state and export to interaction discretization dofmap
  BEAMINTERACTION::UTILS::UpdateDofMapOfVector(
      ia_discret_, ia_state_ptr_->GetDisNp(), global_state().get_dis_np());
  BEAMINTERACTION::UTILS::PeriodicBoundaryConsistentDisVector(ia_state_ptr_->GetDisNp(),
      tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box(), ia_discret_);

  // update column vector
  ia_state_ptr_->GetDisColNp() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  Core::LinAlg::Export(*ia_state_ptr_->GetDisNp(), *ia_state_ptr_->GetDisColNp());

  // update restart displacement vector
  if (ia_state_ptr_->get_restart_coupling_flag())
  {
    ia_state_ptr_->GetDisRestartCol() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
    Core::LinAlg::Export(*ia_state_ptr_->GetDisRestart(), *ia_state_ptr_->GetDisRestartCol());
  }

  // submodel loop
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->reset();

  // Zero out force and stiffness contributions
  force_beaminteraction_->PutScalar(0.0);
  ia_force_beaminteraction_->PutScalar(0.0);
  ia_state_ptr_->GetForceNp()->PutScalar(0.0);
  stiff_beaminteraction_->Zero();
  ia_state_ptr_->GetStiff()->Zero();

  // update gidmap_ and exporter in matrix transform object
  // Note: we need this in every evaluation call (i.e. every iteration) because a change
  //       in the active set of element pairs changes the entries of the used coarse
  //       system stiffness matrix (because we only assemble non-zero values).
  //       Therefore, the graph of the matrix changes and also the required gidmap
  //       (even in computation with one processor)
  // note: this is only necessary if active sets change in consecutive iteration steps
  // ( as crosslinker for example are only updated each time step, we only need to do this
  // every time step)
  if (HaveSubModelType(Inpar::BEAMINTERACTION::submodel_potential) ||
      HaveSubModelType(Inpar::BEAMINTERACTION::submodel_beamcontact))
    update_coupling_adapter_and_matrix_transformation();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::evaluate_force()
{
  check_init_setup();

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->evaluate_force();

  // do communication
  if (ia_state_ptr_->GetForceNp()->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble failed");
  // add to non fe vector
  if (ia_force_beaminteraction_->Update(1., *ia_state_ptr_->GetForceNp(), 1.))
    FOUR_C_THROW("update went wrong");

  // transformation from ia_discret to problem discret
  transform_force();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::evaluate_stiff()
{
  check_init_setup();

  ia_state_ptr_->GetStiff()->UnComplete();

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->evaluate_stiff();

  if (not ia_state_ptr_->GetStiff()->Filled()) ia_state_ptr_->GetStiff()->Complete();

  transform_stiff();

  if (not stiff_beaminteraction_->Filled()) stiff_beaminteraction_->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::evaluate_force_stiff()
{
  check_init_setup();

  ia_state_ptr_->GetStiff()->UnComplete();

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->evaluate_force_stiff();

  // do communication
  if (ia_state_ptr_->GetForceNp()->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble failed");

  // add to non fe vector
  if (ia_force_beaminteraction_->Update(1., *ia_state_ptr_->GetForceNp(), 1.))
    FOUR_C_THROW("update went wrong");
  if (not ia_state_ptr_->GetStiff()->Filled()) ia_state_ptr_->GetStiff()->Complete();

  transform_force_stiff();

  if (not stiff_beaminteraction_->Filled()) stiff_beaminteraction_->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  check_init_setup();

  Core::LinAlg::AssembleMyVector(1.0, f, timefac_np, *force_beaminteraction_);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  check_init_setup();

  Teuchos::RCP<Core::LinAlg::SparseMatrix> jac_dd_ptr = global_state().extract_displ_block(jac);
  jac_dd_ptr->Add(*stiff_beaminteraction_, false, timefac_np, 1.0);

  // no need to keep it
  stiff_beaminteraction_->Zero();
  ia_state_ptr_->GetStiff()->Zero();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();

  int const stepn = global_state().get_step_n();
  double const timen = global_state().get_time_n();
  Teuchos::RCP<Core::IO::DiscretizationWriter> ia_writer = ia_discret_->Writer();
  Teuchos::RCP<Core::IO::DiscretizationWriter> bin_writer = bindis_->Writer();

  // write restart of ia_discret
  ia_writer->write_mesh(stepn, timen);
  ia_writer->new_step(stepn, timen);

  // mesh is not written to disc, only maximum node id is important for output
  // fixme: can we just write mesh
  bin_writer->write_only_nodes_in_new_field_group_to_control_file(stepn, timen, true);
  bin_writer->new_step(stepn, timen);

  // as we know that our maps have changed every time we write output, we can empty
  // the map cache as we can't get any advantage saving the maps anyway
  ia_writer->clear_map_cache();
  bin_writer->clear_map_cache();

  // sub model loop
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->write_restart(*ia_writer, *bin_writer);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  check_init_setup();

  int const stepn = global_state().get_step_n();
  auto input_control_file = Global::Problem::Instance()->InputControlFile();

  // pre sub model loop
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->PreReadRestart();

  // read interaction discretization
  Core::IO::DiscretizationReader ia_reader(ia_discret_, input_control_file, stepn);
  // includes fill_complete()
  ia_reader.read_history_data(stepn);

  // rebuild bin discret correctly in case crosslinker were present
  // Fixme: do just read history data like with ia discret
  // read correct nodes
  Core::IO::DiscretizationReader bin_reader(bindis_, input_control_file, stepn);
  bin_reader.read_nodes_only(stepn);
  bindis_->fill_complete(false, false, false);

  // need to read step next (as it was written next, do safety check)
  if (stepn != ia_reader.read_int("step") or stepn != bin_reader.read_int("step"))
    FOUR_C_THROW("Restart step not consistent with read restart step. ");

  // rebuild binning
  partition_problem();

  // sub model loop
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->read_restart(ia_reader, bin_reader);

  // post sub model loop
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->PostReadRestart();

  // Check if we need to store the restart displacement in the data state container.
  const Teuchos::ParameterList& beam_interaction_params =
      Global::Problem::Instance()->beam_interaction_params();
  if (HaveSubModelType(Inpar::BEAMINTERACTION::submodel_beamcontact) &&
      (bool)Core::UTILS::IntegralValue<int>(
          beam_interaction_params.sublist("BEAM TO SOLID VOLUME MESHTYING"),
          "COUPLE_RESTART_STATE"))
  {
    ia_state_ptr_->set_restart_coupling_flag(true);
    ia_state_ptr_->GetDisRestart() = Teuchos::rcp(new Epetra_Vector(*ia_state_ptr_->GetDisNp()));
    ia_state_ptr_->GetDisRestartCol() = Teuchos::rcp(new Epetra_Vector(*ia_state_ptr_->GetDisNp()));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::run_post_iterate(const ::NOX::Solver::Generic& solver)
{
  check_init_setup();

  // submodel loop
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->run_post_iterate(solver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::update_step_state(const double& timefac_n)
{
  check_init_setup();

  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = global_state().get_fstructure_old();

  fstructold_ptr->Update(timefac_n, *force_beaminteraction_, 1.0);

  // submodel loop
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->UpdateStepState(timefac_n);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::update_step_element()
{
  check_init_setup();

  Vector::iterator sme_iter;

  /* the idea is the following: redistribution of elements is only necessary if
   * one node on any proc has moved "too far" compared to the time step of the
   * last redistribution. Therefore we only do the expensive redistribution,
   * change of ghosting and assigning of elements if necessary.
   * This holds for both the beam discretization as well as the linker/binning
   * discretization.
   */

  // repartition every time
  bool beam_redist = check_if_beam_discret_redistribution_needs_to_be_done();

  // submodel loop
  bool binning_redist = false;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    binning_redist = (*sme_iter)->pre_update_step_element(beam_redist) ? true : binning_redist;

  if (beam_redist)
  {
    binstrategy_->transfer_nodes_and_elements(
        ia_discret_, ia_state_ptr_->GetDisColNp(), ia_state_ptr_->GetBinToRowEleMap());

    build_row_ele_to_bin_map();

    // extend ghosting
    extend_ghosting();

    // assign Elements to bins
    binstrategy_->remove_all_eles_from_bins();
    binstrategy_->AssignElesToBins(ia_discret_, ia_state_ptr_->get_extended_bin_to_row_ele_map());

    // current displacement state gets new reference state
    dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*global_state().get_dis_n()));

    if (global_state().get_my_rank() == 0)
    {
      Core::IO::cout(Core::IO::verbose) << "\n************************************************\n"
                                        << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << "Complete redistribution was done " << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << "\n************************************************\n"
                                        << Core::IO::endl;
    }
  }
  else if (binning_redist)
  {
    extend_ghosting();
    binstrategy_->remove_all_eles_from_bins();
    binstrategy_->AssignElesToBins(ia_discret_, ia_state_ptr_->get_extended_bin_to_row_ele_map());

    if (global_state().get_my_rank() == 0)
    {
      Core::IO::cout(Core::IO::verbose) << "\n************************************************\n"
                                        << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << " binning redistribution was done " << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << "\n************************************************\n"
                                        << Core::IO::endl;
    }
  }

  // update maps of state vectors and matrices
  update_maps();

  // update coupling adapter, this should be done every time step as
  // interacting elements and therefore system matrix can change every time step
  update_coupling_adapter_and_matrix_transformation();

  // submodel loop update
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->UpdateStepElement(binning_redist || beam_redist);

  // submodel post update
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->post_update_step_element();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::check_if_beam_discret_redistribution_needs_to_be_done()
{
  if (beaminteraction_params_ptr_->get_repartition_strategy() !=
      Inpar::BEAMINTERACTION::repstr_adaptive)
    return true;

  Teuchos::RCP<Epetra_Vector> dis_increment =
      Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map(), true));
  int doflid[3];
  for (int i = 0; i < discret_ptr_->NumMyRowNodes(); ++i)
  {
    // get a pointer at i-th row node
    Core::Nodes::Node* node = discret_ptr_->lRowNode(i);

    /* Hermite Interpolation: Check whether node is a beam node which is NOT
     * used for centerline interpolation if so, we simply skip it because
     * it does not have position DoFs */
    if (BEAMINTERACTION::UTILS::IsBeamNode(*node) and
        not BEAMINTERACTION::UTILS::IsBeamCenterlineNode(*node))
      continue;

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret_ptr_->Dof(node);

    for (int dim = 0; dim < 3; ++dim)
    {
      // note: we are using the displacement vector of the problem discretization
      // ( this one also does not get shifted, therefore we do not need to worry
      // about a periodic boundary shift of a node between dis_at_last_redistr_ and the current
      // disp)
      doflid[dim] = dis_at_last_redistr_->Map().LID(dofnode[dim]);
      (*dis_increment)[doflid[dim]] =
          (*global_state().get_dis_np())[doflid[dim]] - (*dis_at_last_redistr_)[doflid[dim]];
    }
  }

  // get maximal displacement increment since last redistribution over all procs
  std::array<double, 2> extrema = {0.0, 0.0};
  dis_increment->MinValue(&extrema[0]);
  dis_increment->MaxValue(&extrema[1]);
  double gmaxdisincr = std::max(-extrema[0], extrema[1]);

  // some verbose screen output
  if (global_state().get_my_rank() == 0)
  {
    Core::IO::cout(Core::IO::debug)
        << " half interaction distance " << half_interaction_distance_ << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << " gmaxdisincr " << gmaxdisincr << Core::IO::endl;
    Core::IO::cout(Core::IO::debug)
        << " half min bin size " << 0.5 * binstrategy_->GetMinBinSize() << Core::IO::endl;
  }

  return ((half_interaction_distance_ + gmaxdisincr) > (0.5 * binstrategy_->GetMinBinSize()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::determine_stress_strain()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::determine_energy()
{
  check_init_setup();

  std::map<STR::EnergyType, double> energy_this_submodel;

  for (auto& submodel : (*me_vec_ptr_))
  {
    energy_this_submodel = submodel->get_energy();

    for (auto const& energy_type : energy_this_submodel)
      eval_data().add_contribution_to_energy_type(energy_type.second, energy_type.first);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::determine_optional_quantity()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->OutputStepState(iowriter);

  // visualize bins according to specification in input file ( MESHFREE -> WRITEBINS "" )
  binstrategy_->WriteBinOutput(global_state().get_step_n(), global_state().get_time_n());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::runtime_output_step_state() const
{
  check_init_setup();

  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->runtime_output_step_state();

  tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box()->runtime_output_step_state(
      global_state().get_time_n(), global_state().get_step_n());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::BeamInteraction::get_block_dof_row_map_ptr()
    const
{
  check_init_setup();
  return global_state().dof_row_map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::BeamInteraction::get_current_solution_ptr()
    const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::BeamInteraction::get_last_time_step_solution_ptr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::post_output()
{
  //  tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box()->ApplyDirichlet(
  //  global_state().get_time_n()
  //  );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::reset_step_state()
{
  Vector::iterator sme_iter;
  for (sme_iter = me_vec_ptr_->begin(); sme_iter != me_vec_ptr_->end(); ++sme_iter)
    (*sme_iter)->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::update_coupling_adapter_and_matrix_transformation()
{
  check_init();

  TEUCHOS_FUNC_TIME_MONITOR(
      "STR::MODELEVALUATOR::BeamInteraction::update_coupling_adapter_and_matrix_transformation");

  // reset transformation member variables (eg. exporter) by rebuilding
  // and provide new maps for coupling adapter
  siatransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  coupsia_->setup_coupling(*ia_discret_, *discret_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::build_row_ele_to_bin_map()
{
  check_init();

  // delete old map
  ia_state_ptr_->GetRowEleToBinMap().clear();
  // loop over bins
  std::map<int, std::set<int>>::const_iterator biniter;
  for (biniter = ia_state_ptr_->GetBinToRowEleMap().begin();
       biniter != ia_state_ptr_->GetBinToRowEleMap().end(); ++biniter)
  {
    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for (eleiter = biniter->second.begin(); eleiter != biniter->second.end(); ++eleiter)
    {
      int elegid = *eleiter;
      int bingid = biniter->first;
      // assign bins to elements
      ia_state_ptr_->GetRowEleToBinMap()[elegid].insert(bingid);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::update_maps()
{
  check_init();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export.
  // todo: check if update is necessary (->SameAs())

  // beam displacement
  BEAMINTERACTION::UTILS::UpdateDofMapOfVector(ia_discret_, ia_state_ptr_->GetDisNp());

  // get current displacement state and export to interaction discretization dofmap
  BEAMINTERACTION::UTILS::UpdateDofMapOfVector(
      ia_discret_, ia_state_ptr_->GetDisNp(), global_state().get_dis_np());
  BEAMINTERACTION::UTILS::PeriodicBoundaryConsistentDisVector(ia_state_ptr_->GetDisNp(),
      tim_int().get_data_sdyn_ptr()->get_periodic_bounding_box(), ia_discret_);

  // update column vector
  ia_state_ptr_->GetDisColNp() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  Core::LinAlg::Export(*ia_state_ptr_->GetDisNp(), *ia_state_ptr_->GetDisColNp());

  // update restart displacement vector
  if (ia_state_ptr_->get_restart_coupling_flag())
  {
    ia_state_ptr_->GetDisRestartCol() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
    Core::LinAlg::Export(*ia_state_ptr_->GetDisRestart(), *ia_state_ptr_->GetDisRestartCol());
  }

  // force
  ia_force_beaminteraction_ = Teuchos::rcp(new Epetra_Vector(*ia_discret_->dof_row_map(), true));
  ia_state_ptr_->GetForceNp() =
      Teuchos::rcp(new Epetra_FEVector(*ia_discret_->dof_row_map(), true));

  // stiff
  ia_state_ptr_->GetStiff() = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *ia_discret_->dof_row_map(), 81, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));

  BEAMINTERACTION::UTILS::SetupEleTypeMapExtractor(ia_discret_, eletypeextractor_);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::transform_force()
{
  check_init();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::BeamInteraction::transform_force");

  // transform force vector to problem discret layout/distribution
  force_beaminteraction_ = coupsia_->MasterToSlave(ia_force_beaminteraction_);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::transform_stiff()
{
  check_init();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::BeamInteraction::transform_stiff");

  stiff_beaminteraction_->UnComplete();
  // transform stiffness matrix to problem discret layout/distribution
  (*siatransform_)(*ia_state_ptr_->GetStiff(), 1.0,
      Core::Adapter::CouplingMasterConverter(*coupsia_), *stiff_beaminteraction_, false);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::transform_force_stiff()
{
  check_init();

  transform_force();
  transform_stiff();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::print_binning_info_to_screen() const
{
  std::vector<Teuchos::RCP<Core::FE::Discretization>> discret_vec(1, ia_discret_);
  std::vector<Teuchos::RCP<const Epetra_Vector>> disnp_vec(1, Teuchos::null);

  double bin_size_lower_bound =
      binstrategy_->compute_lower_bound_for_bin_size_as_max_edge_length_of_aabb_of_largest_ele(
          discret_vec, disnp_vec);
  Core::LinAlg::Matrix<3, 2> XAABB(true);
  binstrategy_->compute_min_binning_domain_containing_all_elements_of_multiple_discrets(
      discret_vec, disnp_vec, XAABB, false);
  if (global_state().get_my_rank() == 0)
  {
    Core::IO::cout(Core::IO::verbose)
        << " \n---------------------------------------------------------- " << Core::IO::endl;
    Core::IO::cout(Core::IO::verbose)
        << " chosen/computed cutoff radius                      : " << binstrategy_->GetMinBinSize()
        << Core::IO::endl;
    Core::IO::cout(Core::IO::verbose)
        << " largest edge length of largest element xaabb       : " << bin_size_lower_bound
        << Core::IO::endl;
    Core::IO::cout(Core::IO::verbose)
        << "DOMAINBOUNDINGBOX containing all elements of input discretization:\n " << XAABB(0, 0)
        << " " << XAABB(1, 0) << " " << XAABB(2, 0) << " " << XAABB(0, 1) << " " << XAABB(1, 1)
        << " " << XAABB(2, 1) << Core::IO::endl;
    Core::IO::cout(Core::IO::verbose)
        << " ----------------------------------------------------------\n " << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Logo() const
{
  check_init();

  if (myrank_ == 0)
  {
    Core::IO::cout << "\n****************************************************************"
                   << Core::IO::endl;
    Core::IO::cout << "*                                                              *"
                   << Core::IO::endl;
    Core::IO::cout << "*          Welcome to the Beam Interaction Model Evaluator     *"
                   << Core::IO::endl;
    Core::IO::cout << "*                                                              *"
                   << Core::IO::endl;
    Core::IO::cout << "****************************************************************"
                   << Core::IO::endl;
    Core::IO::cout << "                                                                  "
                   << Core::IO::endl;
    Core::IO::cout << "                                                                  "
                   << Core::IO::endl;
    Core::IO::cout << "                      0=========================0                 "
                   << Core::IO::endl;
    Core::IO::cout << "                    //|   \\            /       /||                "
                   << Core::IO::endl;
    Core::IO::cout << "                   // |    \\ |       |/       //||                "
                   << Core::IO::endl;
    Core::IO::cout << "                  //  |  /  \\|       /       // ||                "
                   << Core::IO::endl;
    Core::IO::cout << "                 //   |  \\   \\   /  /|\\     //  ||                "
                   << Core::IO::endl;
    Core::IO::cout << "                //    |  /   |\\ /  / | \\   //   ||                "
                   << Core::IO::endl;
    Core::IO::cout << "               //     |  \\   | \\     |  \\ //  / ||                "
                   << Core::IO::endl;
    Core::IO::cout << "              //  \\  /|  /   |/      |   //  /  ||                "
                   << Core::IO::endl;
    Core::IO::cout << "              0=========================0 \\ /   ||                "
                   << Core::IO::endl;
    Core::IO::cout << "             ||    /\\ |____          |  || \\    ||                "
                   << Core::IO::endl;
    Core::IO::cout << "             ||   /  \\|    \\   ------   ||/ \\   ||                "
                   << Core::IO::endl;
    Core::IO::cout << "             ||  /    |                 ||      ||                "
                   << Core::IO::endl;
    Core::IO::cout << "             || /     0----------/------||------0-                "
                   << Core::IO::endl;
    Core::IO::cout << "             ||      /   /       \\      ||     //                 "
                   << Core::IO::endl;
    Core::IO::cout << "             ||     /___/  \\     /    / ||    //                  "
                   << Core::IO::endl;
    Core::IO::cout << "             ||    /        \\    \\   /  ||   //                   "
                   << Core::IO::endl;
    Core::IO::cout << "             ||   /  \\/\\/\\/  \\   /  /   ||  //                    "
                   << Core::IO::endl;
    Core::IO::cout << "             ||  /      /     \\  \\ /    || //                     "
                   << Core::IO::endl;
    Core::IO::cout << "             || /      /         /      ||//                      "
                   << Core::IO::endl;
    Core::IO::cout << "             ||/                       /||/                       "
                   << Core::IO::endl;
    Core::IO::cout << "              0=========================0                         "
                   << Core::IO::endl;
    Core::IO::cout << "                                                                     "
                   << Core::IO::endl;
    Core::IO::cout << "                                                                     "
                   << Core::IO::endl;
  }
}

FOUR_C_NAMESPACE_CLOSE
