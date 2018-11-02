/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_artery_coupling_nodebased.cpp

 \brief base algorithm for node-based coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_nodebased.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_fsi/fsi_matrixtransform.H"


/*----------------------------------------------------------------------*
 | constructor                                         kremheller 04/18 |
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
    PrintOutCouplingMethod();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }

  return;
}

/*----------------------------------------------------------------------*
 | init the strategy                                 kremheller 07/18   |
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
  std::vector<DRT::Condition*> artCoupcond;
  arterydis_->GetCondition(condname_, artCoupcond);

  for (unsigned iter = 0; iter < artCoupcond.size(); ++iter)
  {
    int myID = (artCoupcond[iter])->GetInt("coupling id");
    condIDs[0].push_back(myID);
  }

  // 2) 2D, 3D continuous field discretization
  std::vector<DRT::Condition*> contfieldCoupcond;
  contdis_->GetCondition(condname_, contfieldCoupcond);

  for (unsigned iter = 0; iter < contfieldCoupcond.size(); ++iter)
  {
    int myID = (contfieldCoupcond[iter])->GetInt("coupling id");
    condIDs[1].push_back(myID);
  }

  if (condIDs[0].size() != condIDs[1].size())
    dserror("Artery coupling conditions need to be defined on both discretizations");

  // -----------------------------------------------------------------------------------------------------------------
  // create map extractors needed for artery condition coupling --> continuous field part
  contfieldex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  SetupMapExtractor(contfieldex_, contdis_, coupleddofs_cont_);
  CheckDbcOnCoupledDofs(contdis_, contfieldex_->Map(1));

  // -----------------------------------------------------------------------------------------------------------------
  // create map extractors needed for artery condition coupling --> artery part
  artex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  SetupMapExtractor(artex_, arterydis_, coupleddofs_art_);
  CheckDbcOnCoupledDofs(arterydis_, artex_->Map(1));

  // setup coupling adapter
  artcontfieldcoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  artcontfieldcoup_->SetupConditionCoupling(*contdis_, contfieldex_->Map(1), *arterydis_,
      artex_->Map(1), condname_, coupleddofs_cont_, coupleddofs_art_);

  // full map of continous field and uncoupled dofs of artery
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(contfieldex_->FullMap());
  maps.push_back(artex_->Map(0));

  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  globalex_->Setup(*fullmap_, maps);

  // check global map extractor
  globalex_->CheckForValidMapExtractor();

  // needed for matrix transformations
  sbbtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform());

  return;
}

/*----------------------------------------------------------------------*
 | setup map extractor for artery mesh tying         kremheller 04/18   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupMapExtractor(
    Teuchos::RCP<LINALG::MultiMapExtractor> mapextractor, Teuchos::RCP<DRT::Discretization> dis,
    const std::vector<int>& coupleddofs)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> partialmaps_coupled;

  // build coupled maps for all coupled dofs
  for (int idof = 0; idof < num_coupled_dofs_; idof++)
  {
    DRT::UTILS::MultiConditionSelector mcs;
    Teuchos::RCP<LINALG::MultiMapExtractor> dummy = Teuchos::rcp(new LINALG::MultiMapExtractor());
    // selector for coupleddofs[idof]
    mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(
        *dis, condname_, coupleddofs[idof], coupleddofs[idof] + 1)));
    mcs.SetupExtractor(*dis, *dis->DofRowMap(), *dummy);

    partialmaps_coupled.push_back(dummy->Map(1));
  }
  // fullmap coupled -> all coupled dofs
  Teuchos::RCP<Epetra_Map> fullmap_coupled =
      LINALG::MultiMapExtractor::MergeMaps(partialmaps_coupled);

  // fullmap uncoupled -> all uncoupled dofs
  Teuchos::RCP<LINALG::MapExtractor> temp =
      Teuchos::rcp(new LINALG::MapExtractor(*dis->DofRowMap(), fullmap_coupled, false));
  Teuchos::RCP<Epetra_Map> fullmap_uncoupled = Teuchos::rcp(new Epetra_Map(*temp->CondMap()));

  // vector for setup of extractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fullmap_vector;
  fullmap_vector.push_back(fullmap_uncoupled);
  fullmap_vector.push_back(fullmap_coupled);

  mapextractor->Setup(*dis->DofRowMap(), fullmap_vector);

  return;
}

/*----------------------------------------------------------------------*
 | setup the linear system of equations                kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art,
    Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art)
{
  SetupRHS(rhs, rhs_cont, rhs_art);
  SetupMatrix(sysmat, sysmat_cont, sysmat_art);

  return;
}


/*----------------------------------------------------------------------*
 | setup the global rhs                                kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupRHS(
    Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art)
{
  SetupVector(rhs, rhs_cont, rhs_art);
  return;
}

/*----------------------------------------------------------------------*
 | setup a global vector                               kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupVector(
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

  return;
}

/*----------------------------------------------------------------------*
 | setup system matrix                                 kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::SetupMatrix(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art)
{
  // uncomplete
  sysmat_cont->UnComplete();

  // artery
  // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockartery =
      sysmat_art->Split<LINALG::DefaultBlockMatrixStrategy>(*(artex_), *(artex_));
  blockartery->Complete();

  // inner artery dofs
  sysmat->Assign(1, 1, LINALG::View, blockartery->Matrix(0, 0));

  (*sibtransform_)(blockartery->FullRowMap(), blockartery->FullColMap(), blockartery->Matrix(0, 1),
      1.0, ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), sysmat->Matrix(1, 0));

  (*sbitransform_)(blockartery->Matrix(1, 0), 1.0,
      ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), sysmat->Matrix(0, 1));

  (*sbbtransform_)(blockartery->Matrix(1, 1), 1.0,
      ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_),
      ADAPTER::CouplingSlaveConverter(*artcontfieldcoup_), *sysmat_cont, true, true);

  // continuous field
  sysmat->Assign(0, 0, LINALG::View, *sysmat_cont);
  // complete
  sysmat->Complete();

  return;
}

/*----------------------------------------------------------------------*
 | extract single field vectors                        kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::ExtractSingleFieldVectors(
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

  return;
}

/*----------------------------------------------------------------------*
 | check if nodal coupling and DBC on same DOF         kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::CheckDbcOnCoupledDofs(
    Teuchos::RCP<DRT::Discretization> dis, const Teuchos::RCP<const Epetra_Map>& coupleddofmap)
{
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dis->DofRowMap(), true);
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", 0.0);
    dis->EvaluateDirichlet(eleparams, zeros, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps);
  }
  // intersect DBC maps and coupled dof map to check if coupling and DBC are applied on same dofs
  std::vector<Teuchos::RCP<const Epetra_Map>> dummy;
  dummy.push_back(dbcmaps->CondMap());
  dummy.push_back(coupleddofmap);
  Teuchos::RCP<Epetra_Map> intersect_dbc_coupled = LINALG::MultiMapExtractor::IntersectMaps(dummy);

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
    dserror("Re-think your Input file definition");
  }

  return;
}

/*----------------------------------------------------------------------*
 | check if initial fields on coupled DOFs match       kremheller 04/18 |
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
    dserror("Your initial fields apparently are different with an L2 norm of %f", diff);
  }


  return;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::ArteryDofRowMap() const
{
  return artex_->Map(0);
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::DofRowMap() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::ApplyMeshMovement()
{
  if (!evaluate_in_ref_config_)
    dserror("Evaluation in current configuration not possible for node-based coupling");

  return;
}

/*----------------------------------------------------------------------*
 | print out method                                    kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased::PrintOutCouplingMethod() const
{
  std::cout << "<   Coupling-Method : Nodebased                    >" << std::endl;

  return;
}
