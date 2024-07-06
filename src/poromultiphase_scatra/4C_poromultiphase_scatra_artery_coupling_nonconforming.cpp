/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for non-conforming coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_porofluidmultiphase_ele_parameter.hpp"
#include "4C_porofluidmultiphase_utils.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_defines.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_pair.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
    PoroMultiPhaseScaTraArtCouplNonConforming(Teuchos::RCP<Core::FE::Discretization> arterydis,
        Teuchos::RCP<Core::FE::Discretization> contdis,
        const Teuchos::ParameterList& couplingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplBase(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname),
      couplingparams_(couplingparams),
      condname_(condname),
      porofluidmanagersset_(false),
      issetup_(false),
      porofluidprob_(false),
      has_varying_diam_(false),
      delete_free_hanging_eles_(Core::UTILS::IntegralValue<int>(
          Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
              "ARTERY COUPLING"),
          "DELETE_FREE_HANGING_ELES")),
      delete_free_hanging_eles_threshold_(Global::Problem::instance()
                                              ->poro_fluid_multi_phase_dynamic_params()
                                              .sublist("ARTERY COUPLING")
                                              .get<double>("DELETE_SMALL_FREE_HANGING_COMPS")),
      coupling_method_(Core::UTILS::IntegralValue<
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
          couplingparams, "ARTERY_COUPLING_METHOD")),
      timefacrhs_art_(0.0),
      timefacrhs_cont_(0.0),
      pp_(couplingparams_.get<double>("PENALTY"))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::init()
{
  // we do not have a moving mesh
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::porofluidmultiphase)
  {
    evaluate_in_ref_config_ = true;
    porofluidprob_ = true;
  }

  // fill the vectors
  fill_function_and_scale_vectors();

  // initialize phinp for continuous dis
  phinp_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->dof_row_map(), true));
  // initialize phin for continuous dis
  phin_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->dof_row_map(), true));
  // initialize phinp for artery dis
  phinp_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->dof_row_map(), true));

  // initialize phinp for continuous dis
  zeros_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->dof_row_map(), true));
  // initialize phinp for artery dis
  zeros_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->dof_row_map(), true));

  // -------------------------------------------------------------------
  // create empty D and M matrices (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  d_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(arterydis_->dof_row_map()), 27, false, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  m_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(arterydis_->dof_row_map()), 27, false, true, Core::LinAlg::SparseMatrix::FE_MATRIX));
  kappa_inv_ = Teuchos::rcp(new Epetra_FEVector(*arterydis_->dof_row_map(), true));

  // full map of continous and artery dofs
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(Teuchos::rcp(new Epetra_Map(*contdis_->dof_row_map())));
  maps.push_back(Teuchos::rcp(new Epetra_Map(*arterydis_->dof_row_map())));

  fullmap_ = Core::LinAlg::MultiMapExtractor::merge_maps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
  globalex_->setup(*fullmap_, maps);

  FEmat_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *fullmap_, 81, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));

  fe_rhs_ = Teuchos::rcp(new Epetra_FEVector(*fullmap_));

  // check global map extractor
  globalex_->check_for_valid_map_extractor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::setup()
{
  // get the coupling method
  auto arterycoupl =
      Core::UTILS::IntegralValue<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
          Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
              "ARTERY COUPLING"),
          "ARTERY_COUPLING_METHOD");

  // create the pairs
  if (arterycoupl == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
  {
    get_coupling_idsfrom_input();
    if (couplingnodes_ntp_.size() == 0)
      FOUR_C_THROW("No 1D Coupling Node Ids found for NTP Coupling");
    create_coupling_pairs_ntp();
  }
  else
  {
    create_coupling_pairs_line_surf_based();
  }


  // check if varying diameter is used
  if (contdis_->name() == "porofluid") set_varying_diam_flag();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::get_coupling_idsfrom_input()
{
  // get 1D coupling IDs from Input
  std::vector<Core::Conditions::Condition*> artCoupcond;

  arterydis_->get_condition(condname_, artCoupcond);

  couplingnodes_ntp_.resize(artCoupcond.size());

  for (unsigned iter = 0; iter < artCoupcond.size(); ++iter)
  {
    const std::vector<int>* ArteryNodeIds = (artCoupcond[iter])->get_nodes();
    for (auto couplingids : *ArteryNodeIds)
    {
      couplingnodes_ntp_[iter] = couplingids;
      if (myrank_ == 0)
      {
        std::cout << "Artery Coupling Node Id " << iter + 1
                  << " from Input = " << couplingnodes_ntp_[iter] << "\n";
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::evaluate(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  if (!issetup_) FOUR_C_THROW("setup() has not been called");

  if (!porofluidmanagersset_)
  {
    // set the right-hand side time factors (we assume constant time step size here)
    set_time_fac_rhs();
    for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
      coupl_elepairs_[i]->setup_fluid_managers_and_materials(
          contdis_->name(), timefacrhs_art_, timefacrhs_cont_);
    porofluidmanagersset_ = true;
  }

  // evaluate and assemble the pairs
  evaluate_coupling_pairs(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::setup_system(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_cont,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const Epetra_Map> dbcmap_art,
    Teuchos::RCP<const Epetra_Map> dbcmap_art_with_collapsed)
{
  // add normal part to rhs
  rhs->Update(1.0, *globalex_->insert_vector(rhs_cont, 0), 1.0);
  rhs->Update(1.0, *globalex_->insert_vector(rhs_art, 1), 1.0);

  // apply DBCs
  // 1) on vector
  Core::LinAlg::apply_dirichlet_to_system(*rhs, *zeros_cont_, *(dbcmap_cont->cond_map()));
  Core::LinAlg::apply_dirichlet_to_system(*rhs, *zeros_art_, *(dbcmap_art));
  // 2) on OD-matrices
  sysmat->matrix(0, 1).complete(sysmat_art->range_map(), sysmat_cont->range_map());
  sysmat->matrix(1, 0).complete(sysmat_cont->range_map(), sysmat_art->range_map());
  sysmat->matrix(0, 1).apply_dirichlet(*(dbcmap_cont->cond_map()), false);
  sysmat->matrix(1, 0).apply_dirichlet(*(dbcmap_art_with_collapsed), false);

  // 3) get also the main-diag terms into the global sysmat
  sysmat->matrix(0, 0).add(*sysmat_cont, false, 1.0, 1.0);
  sysmat->matrix(1, 1).add(*sysmat_art, false, 1.0, 1.0);
  sysmat->matrix(0, 0).complete();
  sysmat->matrix(1, 1).complete();
  // and apply DBC
  sysmat->matrix(0, 0).apply_dirichlet(*(dbcmap_cont->cond_map()), true);
  sysmat->matrix(1, 1).apply_dirichlet(*(dbcmap_art_with_collapsed), true);
  // Assign view to 3D system matrix (such that it now includes also contributions from coupling)
  // this is important! Monolithic algorithms use this matrix
  sysmat_cont->assign(Core::LinAlg::View, sysmat->matrix(0, 0));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
    create_coupling_pairs_line_surf_based()
{
  const Teuchos::ParameterList& fluidcouplingparams =
      Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
          "ARTERY COUPLING");
  // loop over pairs found by search
  std::map<int, std::set<int>>::const_iterator nearbyeleiter;
  int numactive_pairs = 0;
  for (nearbyeleiter = nearbyelepairs_.begin(); nearbyeleiter != nearbyelepairs_.end();
       ++nearbyeleiter)
    numactive_pairs += nearbyeleiter->second.size();

  coupl_elepairs_.resize(numactive_pairs);

  int mypair = 0;
  for (nearbyeleiter = nearbyelepairs_.begin(); nearbyeleiter != nearbyelepairs_.end();
       ++nearbyeleiter)
  {
    const int artelegid = nearbyeleiter->first;
    std::vector<Core::Elements::Element const*> ele_ptrs(2);
    ele_ptrs[0] = arterydis_->g_element(artelegid);

    std::set<int>::const_iterator secondeleiter;
    for (secondeleiter = nearbyeleiter->second.begin();
         secondeleiter != nearbyeleiter->second.end(); ++secondeleiter)
    {
      const int contelegid = *secondeleiter;
      ele_ptrs[1] = contdis_->g_element(contelegid);
      if (ele_ptrs[1]->owner() == myrank_)
      {
        // construct, init and setup coupling pairs
        Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPairBase> newpair =
            PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
                create_new_artery_coupling_pair(ele_ptrs);
        newpair->init(ele_ptrs, couplingparams_, fluidcouplingparams, coupleddofs_cont_,
            coupleddofs_art_, scale_vec_, funct_vec_, condname_,
            couplingparams_.get<double>("PENALTY"));

        // add to list of current contact pairs
        coupl_elepairs_[mypair] = newpair;
        mypair++;
      }
    }
  }
  coupl_elepairs_.resize(mypair);

  // output
  int total_numactive_pairs = 0;
  numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  get_comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);


  if (myrank_ == 0)
  {
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;
  }

  // not needed any more
  nearbyelepairs_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::create_coupling_pairs_ntp()
{
  const Teuchos::ParameterList& fluidcouplingparams =
      Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
          "ARTERY COUPLING");

  int numactive_pairs = std::accumulate(nearbyelepairs_.begin(), nearbyelepairs_.end(), 0,
      [](int a, auto b) { return a + (static_cast<int>(b.second.size())); });

  coupl_elepairs_.resize(numactive_pairs);

  // loop over pairs found by search
  int mypair = 0;
  for (const auto& nearbyeleiter : nearbyelepairs_)
  {
    // create vector of active coupling pairs
    std::vector<Core::Elements::Element const*> ele_ptrs(2);
    // assign artery element
    ele_ptrs[0] = arterydis_->g_element(nearbyeleiter.first);

    // get nodes of artery element
    const Core::Nodes::Node* const* artnodes = ele_ptrs[0]->nodes();

    // loop over nodes of artery element
    for (int i = 0; i < ele_ptrs[0]->num_node(); i++)
    {
      // loop over prescribed couplings nodes from input
      for (unsigned int j = 0; j < couplingnodes_ntp_.size(); j++)
      {
        // check if artery node is prescribed coupling node
        if (artnodes[i]->id() == couplingnodes_ntp_[j])
        {
          // get coupling type (ARTERY or AIRWAY ?)
          std::vector<Core::Conditions::Condition*> coupcond;
          arterydis_->get_condition(condname_, coupcond);
          std::string coupling_element_type_ =
              (coupcond[j])->parameters().get<std::string>("coupling_type");

          // recompute coupling dofs
          recompute_coupled_do_fs_for_ntp(coupcond, j);

          // get penalty parameter
          const auto penalty = coupcond[j]->parameters().get<double>("PENALTY");

          // get eta (parameter coordinate of corresponding node)
          const int eta_ntp = (i == 0) ? -1 : 1;

          // loop over assigned 2D/3D elements
          for (const auto continuouseleiter : nearbyeleiter.second)
          {
            // assign 2D/3D element
            ele_ptrs[1] = contdis_->g_element(continuouseleiter);

            // only those pairs, where the 3D element is owned by this proc actually evaluated by
            // this proc
            if (ele_ptrs[1]->owner() == myrank_)
            {
              // construct, init and setup coupling pairs
              Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPairBase>
                  newpair = PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
                      create_new_artery_coupling_pair(ele_ptrs);
              newpair->init(ele_ptrs, couplingparams_, fluidcouplingparams, coupleddofs_cont_,
                  coupleddofs_art_, scale_vec_, funct_vec_, condname_, penalty,
                  coupling_element_type_, eta_ntp);
              // add to list of current contact pairs
              coupl_elepairs_[mypair] = newpair;
              mypair++;
            }
          }
        }
      }
    }
  }
  coupl_elepairs_.resize(mypair);

  // output
  int total_numactive_pairs = 0;
  numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  get_comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);


  if (myrank_ == 0)
  {
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;
  }

  // not needed any more
  nearbyelepairs_.clear();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::set_varying_diam_flag()
{
  int has_varying_diam = 0;
  // check all column elements if one of them uses the diameter law by function
  for (int i = 0; i < arterydis_->num_my_col_elements(); ++i)
  {
    // pointer to current element
    Core::Elements::Element* actele = arterydis_->l_col_element(i);

    // get the artery-material
    Teuchos::RCP<Mat::Cnst1dArt> arterymat =
        Teuchos::rcp_dynamic_cast<Mat::Cnst1dArt>(actele->material());
    if (arterymat == Teuchos::null) FOUR_C_THROW("cast to artery material failed");

    if (arterymat->diameter_law() == Mat::PAR::ArteryDiameterLaw::diameterlaw_by_function)
    {
      has_varying_diam = 1;
      break;
    }
  }

  // sum over all procs.
  int sum_has_varying_diam = 0;
  get_comm().SumAll(&has_varying_diam, &sum_has_varying_diam, 1);
  // if one has a varying diameter set the flag to true
  if (sum_has_varying_diam > 0) has_varying_diam_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::evaluate_coupling_pairs(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  // reset
  if (coupling_method_ == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp)
  {
    d_->zero();
    m_->zero();
    kappa_inv_->PutScalar(0.0);
  }

  FEmat_->zero();
  fe_rhs_->PutScalar(0.0);

  // resulting discrete element force vectors of the two interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  // linearizations
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // element mortar coupling matrices
  Core::LinAlg::SerialDenseMatrix D_ele;
  Core::LinAlg::SerialDenseMatrix M_ele;
  Core::LinAlg::SerialDenseVector Kappa_ele;

  // set states
  if (contdis_->name() == "porofluid")
  {
    contdis_->set_state("phinp_fluid", phinp_cont_);
    contdis_->set_state("phin_fluid", phin_cont_);
    arterydis_->set_state("one_d_artery_pressure", phinp_art_);
    if (not evaluate_in_ref_config_ && not contdis_->has_state(1, "velocity field"))
      FOUR_C_THROW(
          "evaluation in current configuration wanted but solid phase velocity not available!");
    if (has_varying_diam_) reset_integrated_diam_to_zero();
  }
  else if (contdis_->name() == "scatra")
  {
    contdis_->set_state("phinp", phinp_cont_);
    arterydis_->set_state("one_d_artery_phinp", phinp_art_);
  }
  else
    FOUR_C_THROW(
        "Only porofluid and scatra-discretizations are supported for linebased-coupling so far");

  // evaluate all pairs
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // reset state on pairs
    coupl_elepairs_[i]->reset_state(contdis_, arterydis_);

    // get the segment lengths
    const std::vector<double> seglengths = get_ele_segment_lengths(coupl_elepairs_[i]->ele1_gid());

    // evaluate
    const double integrated_diam = coupl_elepairs_[i]->evaluate(&(eleforce[0]), &(eleforce[1]),
        &(elestiff[0][0]), &(elestiff[0][1]), &(elestiff[1][0]), &(elestiff[1][1]), &D_ele, &M_ele,
        &Kappa_ele, seglengths);

    // assemble
    fe_assemble_ele_force_stiff_into_system_vector_matrix(coupl_elepairs_[i]->ele1_gid(),
        coupl_elepairs_[i]->ele2_gid(), integrated_diam, eleforce, elestiff, sysmat, rhs);

    // in case of MP, assemble D, M and Kappa
    if (coupling_method_ == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp and
        num_coupled_dofs_ > 0)
      fe_assemble_dm_kappa(
          coupl_elepairs_[i]->ele1_gid(), coupl_elepairs_[i]->ele2_gid(), D_ele, M_ele, Kappa_ele);
  }

  // set artery diameter in material to be able to evalute the 1D elements with varying diameter
  // and evaluate additional linearization of (integrated) element diameters
  if (contdis_->name() == "porofluid" && has_varying_diam_)
  {
    set_artery_diam_in_material();
    evaluate_additional_linearizationof_integrated_diam();
  }

  if (fe_rhs_->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble of right hand side failed");
  rhs->Update(1.0, *fe_rhs_, 0.0);

  FEmat_->complete();
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockartery =
      FEmat_->split<Core::LinAlg::DefaultBlockMatrixStrategy>(*globalex_, *globalex_);

  blockartery->complete();
  sysmat->matrix(1, 0).add(blockartery->matrix(1, 0), false, 1.0, 0.0);
  sysmat->matrix(0, 1).add(blockartery->matrix(0, 1), false, 1.0, 0.0);
  sysmat->matrix(0, 0).add(blockartery->matrix(0, 0), false, 1.0, 0.0);
  sysmat->matrix(1, 1).add(blockartery->matrix(1, 1), false, 1.0, 0.0);

  // assemble D and M contributions into global force and stiffness
  if (coupling_method_ == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp and
      num_coupled_dofs_ > 0)
    sum_dm_into_global_force_stiff(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
    fe_assemble_ele_force_stiff_into_system_vector_matrix(const int& ele1gid, const int& ele2gid,
        const double& integrated_diam, std::vector<Core::LinAlg::SerialDenseVector> const& elevec,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  const Core::Elements::Element* ele1 = arterydis_->g_element(ele1gid);
  const Core::Elements::Element* ele2 = contdis_->g_element(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->location_vector(*arterydis_, lmrow1, lmrowowner1, lmstride);
  ele2->location_vector(*contdis_, lmrow2, lmrowowner2, lmstride);

  FEmat_->fe_assemble(elemat[0][0], lmrow1, lmrow1);
  FEmat_->fe_assemble(elemat[0][1], lmrow1, lmrow2);
  FEmat_->fe_assemble(elemat[1][0], lmrow2, lmrow1);
  FEmat_->fe_assemble(elemat[1][1], lmrow2, lmrow2);

  fe_rhs_->SumIntoGlobalValues(elevec[0].length(), lmrow1.data(), elevec[0].values());
  fe_rhs_->SumIntoGlobalValues(elevec[1].length(), lmrow2.data(), elevec[1].values());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::fe_assemble_dm_kappa(
    const int& ele1gid, const int& ele2gid, const Core::LinAlg::SerialDenseMatrix& D_ele,
    const Core::LinAlg::SerialDenseMatrix& M_ele, const Core::LinAlg::SerialDenseVector& Kappa_ele)
{
  const Core::Elements::Element* ele1 = arterydis_->g_element(ele1gid);
  const Core::Elements::Element* ele2 = contdis_->g_element(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->location_vector(*arterydis_, lmrow1, lmrowowner1, lmstride);
  ele2->location_vector(*contdis_, lmrow2, lmrowowner2, lmstride);

  d_->fe_assemble(D_ele, lmrow1, lmrow1);
  m_->fe_assemble(M_ele, lmrow1, lmrow2);
  kappa_inv_->SumIntoGlobalValues(Kappa_ele.length(), lmrow1.data(), Kappa_ele.values());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
    sum_dm_into_global_force_stiff(
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  // invert
  invert_kappa();

  // complete
  d_->complete();
  m_->complete(*contdis_->dof_row_map(), *arterydis_->dof_row_map());

  // get kappa matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kappaInvMat =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*new Epetra_Vector(Copy, *kappa_inv_, 0)));
  kappaInvMat->complete();

  // kappa^{-1}*M
  Teuchos::RCP<Core::LinAlg::SparseMatrix> km =
      Core::LinAlg::MLMultiply(*kappaInvMat, false, *m_, false, false, false, true);
  // kappa^{-1}*D
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kd =
      Core::LinAlg::MLMultiply(*kappaInvMat, false, *d_, false, false, false, true);

  // D^T*kappa^{-1}*D
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dtkd =
      Core::LinAlg::MLMultiply(*d_, true, *kd, false, false, false, true);
  // D^T*kappa^{-1}*M
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dtkm =
      Core::LinAlg::MLMultiply(*d_, true, *km, false, false, false, true);
  // M^T*kappa^{-1}*M
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mtkm =
      Core::LinAlg::MLMultiply(*m_, true, *km, false, false, false, true);

  // add matrices
  sysmat->matrix(0, 0).add(*mtkm, false, pp_ * timefacrhs_cont_, 1.0);
  sysmat->matrix(1, 1).add(*dtkd, false, pp_ * timefacrhs_art_, 1.0);
  sysmat->matrix(1, 0).add(*dtkm, false, -pp_ * timefacrhs_art_, 1.0);
  sysmat->matrix(0, 1).add(*dtkm, true, -pp_ * timefacrhs_cont_, 1.0);

  // add vector
  Teuchos::RCP<Epetra_Vector> art_contribution =
      Teuchos::rcp(new Epetra_Vector(*arterydis_->dof_row_map()));
  Teuchos::RCP<Epetra_Vector> cont_contribution =
      Teuchos::rcp(new Epetra_Vector(*contdis_->dof_row_map()));

  // Note: all terms are negative since rhs
  // pp*D^T*kappa^{-1}*D*phi_np^art
  dtkd->multiply(false, *phinp_art_, *art_contribution);
  rhs->Update(-pp_ * timefacrhs_art_, *globalex_->insert_vector(art_contribution, 1), 1.0);

  // -pp*D^T*kappa^{-1}*M*phi_np^cont
  dtkm->multiply(false, *phinp_cont_, *art_contribution);
  rhs->Update(pp_ * timefacrhs_art_, *globalex_->insert_vector(art_contribution, 1), 1.0);

  // pp*M^T*kappa^{-1}*M*phi_np^cont
  mtkm->multiply(false, *phinp_cont_, *cont_contribution);
  rhs->Update(-pp_ * timefacrhs_cont_, *globalex_->insert_vector(cont_contribution, 0), 1.0);

  // -pp*M^T*kappa^{-1}*D*phi_np^art = -pp*(D^T*kappa^{-1}*M)^T*phi_np^art
  dtkm->multiply(true, *phinp_art_, *cont_contribution);
  rhs->Update(pp_ * timefacrhs_cont_, *globalex_->insert_vector(cont_contribution, 0), 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::invert_kappa()
{
  // global assemble
  if (kappa_inv_->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble of kappaInv_ failed");

  // invert (pay attention to protruding elements)
  for (int i = 0; i < arterydis_->dof_row_map()->NumMyElements(); ++i)
  {
    const int artdofgid = arterydis_->dof_row_map()->GID(i);
    const double kappaVal = (*kappa_inv_)[0][kappa_inv_->Map().LID(artdofgid)];
    if (fabs(kappaVal) > KAPPAINVTOL)
      kappa_inv_->ReplaceGlobalValue(artdofgid, 0, 1.0 / kappaVal);
    else
      kappa_inv_->ReplaceGlobalValue(artdofgid, 0, 0.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPairBase>
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::create_new_artery_coupling_pair(
    std::vector<Core::Elements::Element const*> const& ele_ptrs)
{
  const Core::FE::CellType distypeart = ele_ptrs[0]->shape();
  switch (distypeart)
  {
    case Core::FE::CellType::line2:
    {
      const Core::FE::CellType distypecont = ele_ptrs[1]->shape();
      switch (distypecont)
      {
        case Core::FE::CellType::quad4:
        {
          switch (Global::Problem::instance()->n_dim())
          {
            case 1:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::quad4, 1>());
            case 2:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::quad4, 2>());
            case 3:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::quad4, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", Global::Problem::instance()->n_dim());
          }
        }
        case Core::FE::CellType::hex8:
        {
          switch (Global::Problem::instance()->n_dim())
          {
            case 1:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::hex8, 1>());
            case 2:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::hex8, 2>());
            case 3:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::hex8, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", Global::Problem::instance()->n_dim());
          }
        }
        case Core::FE::CellType::tet4:
        {
          switch (Global::Problem::instance()->n_dim())
          {
            case 1:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet4, 1>());
            case 2:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet4, 2>());
            case 3:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet4, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", Global::Problem::instance()->n_dim());
          }
        }
        case Core::FE::CellType::tet10:
        {
          switch (Global::Problem::instance()->n_dim())
          {
            case 1:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet10, 1>());
            case 2:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet10, 2>());
            case 3:
              return Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPair<
                  Core::FE::CellType::line2, Core::FE::CellType::tet10, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", Global::Problem::instance()->n_dim());
          }
        }
        default:
          FOUR_C_THROW(
              "only quad4, hex8, tet4 and tet10 elements supported for continuous elements so far");
      }
    }
    default:
      FOUR_C_THROW("only line 2 elements supported for artery elements so far");
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::setup_vector(
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
    Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // zero out
  vec->PutScalar(0.0);
  // set up global vector
  globalex_->insert_vector(*vec_cont, 0, *vec);
  globalex_->insert_vector(*vec_art, 1, *vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::extract_single_field_vectors(
    Teuchos::RCP<const Epetra_Vector> globalvec, Teuchos::RCP<const Epetra_Vector>& vec_cont,
    Teuchos::RCP<const Epetra_Vector>& vec_art)
{
  // process first field (continuous)
  vec_cont = globalex_->extract_vector(globalvec, 0);
  // process second field (artery)
  vec_art = globalex_->extract_vector(globalvec, 1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::artery_dof_row_map() const
{
  return globalex_->Map(1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::dof_row_map() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::set_solution_vectors(
    Teuchos::RCP<const Epetra_Vector> phinp_cont, Teuchos::RCP<const Epetra_Vector> phin_cont,
    Teuchos::RCP<const Epetra_Vector> phinp_art)
{
  phinp_cont_ = phinp_cont;
  if (phin_cont != Teuchos::null) phin_cont_ = phin_cont;
  phinp_art_ = phinp_art;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::blood_vessel_volume_fraction()
{
  FOUR_C_THROW("Not implemented in base class");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::print_out_coupling_method()
    const
{
  std::string name;
  if (coupling_method_ == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp)
    name = "Mortar Penalty";
  else if (coupling_method_ == Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::gpts)
    name = "Gauss-Point-To-Segment";
  else
    FOUR_C_THROW("unknown coupling method");

  std::cout << "<   Coupling-Method : " << std::left << std::setw(22) << name << "       >"
            << std::endl;
  std::cout << "<   Penalty         : " << std::left << std::setw(6) << pp_
            << "                       >" << std::endl;
  if (evaluate_in_ref_config_)
    std::cout << "<   Moving arteries : No                           >" << std::endl;
  else
    std::cout << "<   Moving arteries : Yes                          >" << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::
    fill_function_and_scale_vectors()
{
  scale_vec_.resize(2);
  funct_vec_.resize(2);

  // get the actual coupled DOFs  ----------------------------------------------------
  // 1) 1D artery discretization
  int word1;
  std::istringstream scale_art_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "SCALEREAC_ART"));
  while (scale_art_stream >> word1) scale_vec_[0].push_back((int)(word1));

  std::istringstream funct_art_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "REACFUNCT_ART"));
  while (funct_art_stream >> word1) funct_vec_[0].push_back((int)(word1 - 1));

  // 2) 2D, 3D continuous field discretization
  std::istringstream scale_cont_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "SCALEREAC_CONT"));
  while (scale_cont_stream >> word1) scale_vec_[1].push_back((int)(word1));

  std::istringstream funct_cont_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "REACFUNCT_CONT"));
  while (funct_cont_stream >> word1) funct_vec_[1].push_back((int)(word1 - 1));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::set_time_fac_rhs()
{
  // set the right hand side factor
  if (contdis_->name() == "porofluid")
  {
    Discret::ELEMENTS::PoroFluidMultiPhaseEleParameter* eleparams =
        Discret::ELEMENTS::PoroFluidMultiPhaseEleParameter::instance("porofluid");
    // artery
    timefacrhs_art_ = 1.0;
    // continuous
    timefacrhs_cont_ = eleparams->time_fac_rhs();
  }
  else if (contdis_->name() == "scatra")
  {
    Discret::ELEMENTS::ScaTraEleParameterTimInt* eleparams =
        Discret::ELEMENTS::ScaTraEleParameterTimInt::instance("scatra");
    // artery
    timefacrhs_art_ = eleparams->time_fac_rhs();
    // continuous
    timefacrhs_cont_ = eleparams->time_fac_rhs();
  }
  else
  {
    FOUR_C_THROW(
        "Only porofluid and scatra-discretizations are supported for non-conforming coupling so "
        "far");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNonConforming::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  nearbyelepairs_ = *nearbyelepairs;
}

FOUR_C_NAMESPACE_CLOSE
