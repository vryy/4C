/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for node-based coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_artery_coupling_nodebased.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_condition_selector.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_matrixtransform.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::PoroMultiPhaseScaTraArtCouplNodeBased(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& meshtyingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplBase(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname),
      condname_(condname)
{
  // user info
  if (myrank_ == 0)
  {
    std::cout << "<                                                  >" << std::endl;
    print_out_coupling_method();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::Init()
{
  // ARTERY COUPLING CONDITIONS
  std::vector<std::vector<int>> condIDs;
  std::vector<int> artIDs;
  std::vector<int> contfieldIDs;
  condIDs.push_back(artIDs);
  condIDs.push_back(contfieldIDs);

  // check if conditions are defined on both discretizations --------------------------
  // 1) 1D artery discretization
  std::vector<CORE::Conditions::Condition*> artCoupcond;
  arterydis_->GetCondition(condname_, artCoupcond);

  for (auto& iter : artCoupcond)
  {
    int myID = iter->parameters().Get<int>("coupling id");
    condIDs[0].push_back(myID);
  }

  // 2) 2D, 3D continuous field discretization
  std::vector<CORE::Conditions::Condition*> contfieldCoupcond;
  contdis_->GetCondition(condname_, contfieldCoupcond);

  for (auto& iter : contfieldCoupcond)
  {
    int myID = iter->parameters().Get<int>("coupling id");
    condIDs[1].push_back(myID);
  }

  if (condIDs[0].size() != condIDs[1].size())
    FOUR_C_THROW("Artery coupling conditions need to be defined on both discretizations");

  // -----------------------------------------------------------------------------------------------------------------
  // create map extractors needed for artery condition coupling --> continuous field part
  contfieldex_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
  setup_map_extractor(contfieldex_, contdis_, coupleddofs_cont_);
  check_dbc_on_coupled_dofs(contdis_, contfieldex_->Map(1));

  // -----------------------------------------------------------------------------------------------------------------
  // create map extractors needed for artery condition coupling --> artery part
  artex_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
  setup_map_extractor(artex_, arterydis_, coupleddofs_art_);
  check_dbc_on_coupled_dofs(arterydis_, artex_->Map(1));

  // setup coupling adapter
  artcontfieldcoup_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  artcontfieldcoup_->setup_condition_coupling(*contdis_, contfieldex_->Map(1), *arterydis_,
      artex_->Map(1), condname_, coupleddofs_cont_, coupleddofs_art_);

  // full map of continous field and uncoupled dofs of artery
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(contfieldex_->FullMap());
  maps.push_back(artex_->Map(0));

  fullmap_ = CORE::LINALG::MultiMapExtractor::MergeMaps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
  globalex_->Setup(*fullmap_, maps);

  // check global map extractor
  globalex_->check_for_valid_map_extractor();

  // needed for matrix transformations
  sbbtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::Setup()
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::setup_map_extractor(
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mapextractor,
    Teuchos::RCP<DRT::Discretization> dis, const std::vector<int>& coupleddofs)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> partialmaps_coupled;

  // build coupled maps for all coupled dofs
  for (int idof = 0; idof < num_coupled_dofs_; idof++)
  {
    CORE::Conditions::MultiConditionSelector mcs;
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> dummy =
        Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
    // selector for coupleddofs[idof]
    mcs.AddSelector(Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(
        *dis, condname_, coupleddofs[idof], coupleddofs[idof] + 1)));
    mcs.SetupExtractor(*dis, *dis->dof_row_map(), *dummy);

    partialmaps_coupled.push_back(dummy->Map(1));
  }
  // fullmap coupled -> all coupled dofs
  Teuchos::RCP<Epetra_Map> fullmap_coupled =
      CORE::LINALG::MultiMapExtractor::MergeMaps(partialmaps_coupled);

  // fullmap uncoupled -> all uncoupled dofs
  Teuchos::RCP<CORE::LINALG::MapExtractor> temp =
      Teuchos::rcp(new CORE::LINALG::MapExtractor(*dis->dof_row_map(), fullmap_coupled, false));
  Teuchos::RCP<Epetra_Map> fullmap_uncoupled = Teuchos::rcp(new Epetra_Map(*temp->CondMap()));

  // vector for setup of extractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fullmap_vector;
  fullmap_vector.push_back(fullmap_uncoupled);
  fullmap_vector.push_back(fullmap_coupled);

  mapextractor->Setup(*dis->dof_row_map(), fullmap_vector);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupSystem(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art)
{
  setup_rhs(rhs, rhs_cont, rhs_art);
  setup_matrix(sysmat, sysmat_cont, sysmat_art);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::setup_rhs(
    Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art)
{
  setup_vector(rhs, rhs_cont, rhs_art);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::setup_vector(
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
    Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // zero out
  vec->PutScalar(0.0);

  // inner (uncoupled) DOFs of artery
  Teuchos::RCP<Epetra_Vector> vec2_uncoupled = artex_->ExtractVector(vec_art, 0);

  // boundary (coupled) DOFs of artery
  Teuchos::RCP<Epetra_Vector> vec2_coupled = artex_->ExtractVector(vec_art, 1);

  // transform boundary DOFs to continuous dis
  Teuchos::RCP<Epetra_Vector> temp =
      contfieldex_->InsertVector(artcontfieldcoup_->SlaveToMaster(vec2_coupled), 1);

  // add to continous vec
  temp->Update(1.0, *vec_cont, 1.0);

  // set up global vector
  globalex_->InsertVector(*temp, 0, *vec);
  globalex_->InsertVector(*vec2_uncoupled, 1, *vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::setup_matrix(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art)
{
  // uncomplete
  sysmat_cont->UnComplete();

  // artery
  // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockartery =
      sysmat_art->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(*(artex_), *(artex_));
  blockartery->Complete();

  // inner artery dofs
  sysmat->Assign(1, 1, CORE::LINALG::View, blockartery->Matrix(0, 0));

  (*sibtransform_)(blockartery->FullRowMap(), blockartery->FullColMap(), blockartery->Matrix(0, 1),
      1.0, CORE::ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), sysmat->Matrix(1, 0));

  (*sbitransform_)(blockartery->Matrix(1, 0), 1.0,
      CORE::ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), sysmat->Matrix(0, 1));

  (*sbbtransform_)(blockartery->Matrix(1, 1), 1.0,
      CORE::ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_),
      CORE::ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), *sysmat_cont, true, true);

  // continuous field
  sysmat->Assign(0, 0, CORE::LINALG::View, *sysmat_cont);
  // complete
  sysmat->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::extract_single_field_vectors(
    Teuchos::RCP<const Epetra_Vector> globalvec, Teuchos::RCP<const Epetra_Vector>& vec_cont,
    Teuchos::RCP<const Epetra_Vector>& vec_art)
{
  // process second field (continuous)
  vec_cont = globalex_->ExtractVector(globalvec, 0);

  // process coupled (boundary) DOFs of the second field
  Teuchos::RCP<Epetra_Vector> boundary = contfieldex_->ExtractVector(vec_cont, 1);

  // process inner (uncoupled) and boundary (coupled) DOFs of artery
  Teuchos::RCP<const Epetra_Vector> artery_inner = globalex_->ExtractVector(globalvec, 1);
  Teuchos::RCP<Epetra_Vector> artery_boundary = artcontfieldcoup_->MasterToSlave(boundary);

  // build vector for artery
  // 1) inner DOFs
  Teuchos::RCP<Epetra_Vector> artery_temp = artex_->InsertVector(artery_inner, 0);
  // 2) boundary DOFs
  artex_->InsertVector(artery_boundary, 1, artery_temp);

  vec_art = artery_temp;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::check_dbc_on_coupled_dofs(
    Teuchos::RCP<DRT::Discretization> dis, const Teuchos::RCP<const Epetra_Map>& coupleddofmap)
{
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps = Teuchos::rcp(new CORE::LINALG::MapExtractor());
  {
    Teuchos::RCP<Epetra_Vector> zeros = CORE::LINALG::CreateVector(*dis->dof_row_map(), true);
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", 0.0);
    eleparams.set<const CORE::UTILS::FunctionManager*>(
        "function_manager", &GLOBAL::Problem::Instance()->FunctionManager());
    dis->evaluate_dirichlet(eleparams, zeros, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps);
  }
  // intersect DBC maps and coupled dof map to check if coupling and DBC are applied on same dofs
  std::vector<Teuchos::RCP<const Epetra_Map>> dummy;
  dummy.push_back(dbcmaps->CondMap());
  dummy.push_back(coupleddofmap);
  Teuchos::RCP<Epetra_Map> intersect_dbc_coupled =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(dummy);

  if (intersect_dbc_coupled->NumGlobalElements() > 0)
  {
    if (myrank_ == 0)
    {
      std::cout << "\n\n";
      std::cout << "You cannot define DBC and nodal coupling conditions on the same node\n"
                   "for discretization "
                << dis->Name()
                << "\n"
                   "The problematic DOFs are:"
                << std::endl;
    }
    intersect_dbc_coupled->Print(std::cout);
    FOUR_C_THROW("Re-think your Input file definition");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::CheckInitialFields(
    Teuchos::RCP<const Epetra_Vector> vec_cont, Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // boundary (coupled) DOFs of artery
  Teuchos::RCP<Epetra_Vector> vec2_coupled = artex_->ExtractVector(vec_art, 1);

  // transform boundary DOFs to continuous dis
  Teuchos::RCP<Epetra_Vector> temp = artcontfieldcoup_->SlaveToMaster(vec2_coupled);

  // process coupled (boundary) DOFs of the second field
  Teuchos::RCP<Epetra_Vector> boundary = contfieldex_->ExtractVector(vec_cont, 1);

  // subtract artery DOF values from continuous DOF values
  boundary->Update(-1.0, *temp, 1.0);

  // build L2 norm
  double diff(0.0);
  boundary->Norm2(&diff);

  if (diff > 1.0e-9)
  {
    FOUR_C_THROW("Your initial fields apparently are different with an L2 norm of %f", diff);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::ArteryDofRowMap() const
{
  return artex_->Map(0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::dof_row_map() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::ApplyMeshMovement()
{
  if (!evaluate_in_ref_config_)
    FOUR_C_THROW("Evaluation in current configuration not possible for node-based coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::blood_vessel_volume_fraction()
{
  FOUR_C_THROW("Output of vessel volume fraction not possible for node-based coupling");

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::print_out_coupling_method() const
{
  std::cout << "<   Coupling-Method : Nodebased                    >" << std::endl;
}

FOUR_C_NAMESPACE_CLOSE
