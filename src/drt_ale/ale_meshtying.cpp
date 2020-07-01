/*--------------------------------------------------------------------------*/
/*! \file

\brief Mesh tying for ale problems

\level 2

*/
/*--------------------------------------------------------------------------*/


#define DIRECTMANIPULATION
#define ZEROSYSMAT

#include "ale_meshtying.H"
#include "ale_utils.H"
#include "ale_utils_mapextractor.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_krylov_projector.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>


ALE::Meshtying::Meshtying(Teuchos::RCP<DRT::Discretization> dis, LINALG::Solver& solver, int msht,
    int nsd, const UTILS::MapExtractor* surfacesplitter)
    : discret_(dis),
      solver_(solver),
      dofrowmap_(discret_->DofRowMap()),
      gsdofrowmap_(Teuchos::null),
      gmdofrowmap_(Teuchos::null),
      mergedmap_(Teuchos::null),
      //  msht_(msht),
      surfacesplitter_(surfacesplitter),
      problemrowmap_(Teuchos::null),
      gndofrowmap_(Teuchos::null),
      gsmdofrowmap_(Teuchos::null),
      valuesdc_(Teuchos::null),
      dconmaster_(false),
      firstnonliniter_(false),
      //  nsd_(nsd),
      is_multifield_(false)
{
  // get the processor ID from the communicator
  myrank_ = discret_->Comm().MyPID();
}

/*-------------------------------------------------------*/
/*  Setup mesh-tying problem                 wirtz 01/16 */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> ALE::Meshtying::Setup(
    std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  1)   Setup Meshtying");
  //  if(coupleddof[nsd_]==0)
  //    pcoupled_=false;

  AdapterMortar(coupleddof);

  if (myrank_ == 0) CompareNumDof();

  DofRowMaps();

  // merge dofrowmap for slave and master discretization
  gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_, *gsdofrowmap_, false);

  // dofrowmap for discretisation without slave and master dofrowmap
  gndofrowmap_ = LINALG::SplitMap(*dofrowmap_, *gsmdofrowmap_);

  // map for 2x2 (uncoupled dof's & master dof's)
  mergedmap_ = LINALG::MergeMap(*gndofrowmap_, *gmdofrowmap_, false);

  // std::cout << "number of n dof   " << gndofrowmap_->NumGlobalElements() << std::endl;
  // std::cout << "number of m dof   " << gmdofrowmap_->NumGlobalElements() << std::endl;
  // std::cout << "number of s dof   " << gsdofrowmap_->NumGlobalElements() << std::endl;

  // generate map for blockmatrix
  std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
  alemaps.push_back(gndofrowmap_);
  alemaps.push_back(gmdofrowmap_);
  alemaps.push_back(gsdofrowmap_);

  LINALG::MultiMapExtractor extractor;

  extractor.Setup(*dofrowmap_, alemaps);

  // check, if extractor maps are valid
  extractor.CheckForValidMapExtractor();

  // allocate 3x3 block sparse matrix with the interface split strategy
  // the interface split strategy speeds up the assembling process,
  // since the information, which nodes are part of the interface, is available
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------

  Teuchos::RCP<LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
  mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
      extractor, extractor, 108, false, true));
  // nodes on the interface
  Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->ConditionedElementMap(*discret_);

  mat->SetCondElements(condelements);

  // Important: right way to do it (Tobias W.)
  // allocate 2x2 solution matrix with the default block matrix strategy in order to solve the
  // reduced system memory is not allocated(1), since the matrix gets a Teuchos::RCP on the
  // respective blocks of the 3x3 block matrix
  // ---------------
  // | knn  | knm' |
  // | kmn' | kmm' |
  // ---------------

  LINALG::MapExtractor rowmapext(*mergedmap_, gmdofrowmap_, gndofrowmap_);
  LINALG::MapExtractor dommapext(*mergedmap_, gmdofrowmap_, gndofrowmap_);
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>> matsolve =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          dommapext, rowmapext, 1, false, true));
  sysmatsolve_ = matsolve;

  return mat;
}

/*-------------------------------------------------------*/
/*  Use the split of the ale mesh tying for the sysmat   */
/*                                           wirtz 01/16 */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> ALE::Meshtying::MshtSplit()
{
  // generate map for blockmatrix
  std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
  alemaps.push_back(gndofrowmap_);
  alemaps.push_back(gmdofrowmap_);
  alemaps.push_back(gsdofrowmap_);

  LINALG::MultiMapExtractor extractor;

  extractor.Setup(*dofrowmap_, alemaps);

  // check, if extractor maps are valid
  extractor.CheckForValidMapExtractor();

  // allocate 3x3 block sparse matrix with the interface split strategy
  // the interface split strategy speeds up the assembling process,
  // since the information, which nodes are part of the interface, is available
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------

  Teuchos::RCP<LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
  mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
      extractor, extractor, 108, false, true));
  // nodes on the interface
  Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->ConditionedElementMap(*discret_);
  mat->SetCondElements(condelements);

  return mat;
}

/*------------------------------------------------------------------------------*/
/*  Check if Dirichlet BC are defined on the master                 wirtz 01/16 */
/*------------------------------------------------------------------------------*/
void ALE::Meshtying::DirichletOnMaster(Teuchos::RCP<const Epetra_Map> bmaps)
{
  // This method checks if Dirichlet or Dirichlet-like boundary conditions are defined
  // on the master side of the internal interface.
  // In this case, the slave side has to be handled in a special way
  // strategies:
  // (a)  Apply DC on both master and slave side of the internal interface (->disabled)
  //      -> over-constraint system, but nevertheless, result is correct and no solver issues
  // (b)  DC are projected from the master to the slave side during PrepareTimeStep
  //      (in ProjectMasterToSlaveForOverlappingBC()) (-> disabled)
  //      -> DC also influence slave nodes which are not part of the inflow
  //
  //      if(msht_ != INPAR::ALE::no_meshtying)
  //        meshtying_->ProjectMasterToSlaveForOverlappingBC(dispnp_, dbcmaps_->CondMap());
  //
  // (c)  DC are included in the condensation process (-> actual strategy)

  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(bmaps);
  Teuchos::RCP<const Epetra_Map> gmdofrowmap = gmdofrowmap_;
  intersectionmaps.push_back(gmdofrowmap);
  Teuchos::RCP<Epetra_Map> intersectionmap =
      LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    dconmaster_ = true;
    if (myrank_ == 0)
    {
      std::cout
          << "Dirichlet or Dirichlet-like boundary condition defined on master side of the "
             "internal interface!\n "
          << "These conditions has to be also included at the slave side of the internal interface"
          << std::endl
          << std::endl;
    }
  }

  return;
}

/*---------------------------------------------------*/
/*  Prepare Meshtying system             wirtz 01/16 */
/*---------------------------------------------------*/
void ALE::Meshtying::PrepareMeshtyingSystem(Teuchos::RCP<LINALG::SparseOperator>& sysmat,
    Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation block matrix");

  CondensationOperationBlockMatrix(sysmat, residual, dispnp);
  return;
}

/*-------------------------------------------------------*/
/*  Split Vector                             wirtz 01/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::SplitVector(
    Teuchos::RCP<Epetra_Vector> vector, std::vector<Teuchos::RCP<Epetra_Vector>>& splitvector)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.2)   - Split Vector");

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/

  // do the vector splitting smn -> sm+n
  LINALG::SplitVector(*dofrowmap_, *vector, gsmdofrowmap_, fsm, gndofrowmap_, fn);

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  LINALG::SplitVector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

  // splitvector[ii]
  // fn [0]
  // fm [1]
  // fs [2]

  splitvector[0] = fn;
  splitvector[1] = fm;
  splitvector[2] = fs;

  return;
}

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
void ALE::Meshtying::SplitVectorBasedOn3x3(
    Teuchos::RCP<Epetra_Vector> orgvector, Teuchos::RCP<Epetra_Vector> vectorbasedon2x2)
{
  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvector(3);

  SplitVector(orgvector, splitvector);
  // build up the reduced residual
  LINALG::Export(*(splitvector[0]), *vectorbasedon2x2);
  LINALG::Export(*(splitvector[1]), *vectorbasedon2x2);

  return;
}



/*-------------------------------------------------------*/
/*  Set the flag for multifield problems     wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshtying::IsMultifield(
    const LINALG::MultiMapExtractor& interface,  ///< interface maps for split of ale matrix
    bool ismultifield                            ///< flag for multifield problems
)
{
  multifield_interface_ = interface;
  is_multifield_ = ismultifield;

  return;
}

/*-------------------------------------------------------*/
/*  Use the split of the ale mesh tying for the sysmat   */
/*                                           wirtz 01/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::MshtSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat)
{
  if (is_multifield_)
  {
    // generate map for blockmatrix
    std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
    alemaps.push_back(gndofrowmap_);
    alemaps.push_back(gmdofrowmap_);
    alemaps.push_back(gsdofrowmap_);

    LINALG::MultiMapExtractor extractor;

    extractor.Setup(*dofrowmap_, alemaps);

    // check, if extractor maps are valid
    extractor.CheckForValidMapExtractor();

    // allocate 3x3 block sparse matrix with the interface split strategy
    // the interface split strategy speeds up the assembling process,
    // since the information, which nodes are part of the interface, is available
    // -------------------
    // | knn | knm | kns |
    // | kmn | kmm | kms |
    // | ksn | ksm | kss |
    // -------------------

    Teuchos::RCP<LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
        extractor, extractor, 108, false, true));
    // nodes on the interface
    Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->ConditionedElementMap(*discret_);
    mat->SetCondElements(condelements);

    sysmat = mat;
  }
}

/*-------------------------------------------------------*/
/*  Use the split of the multifield problem for the      */
/*  sysmat                                   wirtz 01/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::MultifieldSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat)
{
  if (is_multifield_)
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew =
        Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

    Teuchos::RCP<Epetra_Vector> ones =
        Teuchos::rcp(new Epetra_Vector(sysmatnew->Matrix(2, 2).RowMap()));
    ones->PutScalar(1.0);

    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    sysmatnew->Matrix(0, 2).UnComplete();
    sysmatnew->Matrix(0, 2).Zero();

    sysmatnew->Matrix(1, 2).UnComplete();
    sysmatnew->Matrix(1, 2).Zero();

    sysmatnew->Matrix(2, 2).UnComplete();
    sysmatnew->Matrix(2, 2).Zero();
    sysmatnew->Matrix(2, 2).Add(*onesdiag, false, 1.0, 1.0);

    sysmatnew->Matrix(2, 0).UnComplete();
    sysmatnew->Matrix(2, 0).Zero();

    sysmatnew->Matrix(2, 1).UnComplete();
    sysmatnew->Matrix(2, 1).Zero();

    sysmatnew->Complete();

    Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = sysmatnew->Merge();

    Teuchos::RCP<LINALG::MapExtractor> extractor = Teuchos::rcp(
        new LINALG::MapExtractor(*multifield_interface_.FullMap(), multifield_interface_.Map(1)));

    Teuchos::RCP<LINALG::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat =
        mergedmatrix->Split<ALE::UTILS::InterfaceSplitStrategy>(*extractor, *extractor);

    mat->Complete();

    sysmat = mat;
  }
}

/*-------------------------------------------------------*/
/*  Call the constructor and the setup of the mortar     */
/*  coupling adapter                         wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::AdapterMortar(std::vector<int> coupleddof)
{
  adaptermeshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());

  // Setup of meshtying adapter
  adaptermeshtying_->Setup(
      discret_, discret_, Teuchos::null, coupleddof, "Mortar", discret_->Comm(), true);
}

/*-------------------------------------------------------*/
/*  Compare the size of the slave and master dof row map */
/*                                           wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::CompareNumDof()
{
  int numdofmaster = (adaptermeshtying_->MasterDofMap())->NumGlobalElements();
  int numdofslave = (adaptermeshtying_->SlaveDofMap())->NumGlobalElements();

  std::cout << std::endl << "number of master dof's:   " << numdofmaster << std::endl;
  std::cout << "number of slave dof's:   " << numdofslave << std::endl << std::endl;

  if (numdofmaster > numdofslave)
    std::cout << "The master side is discretized by more elements than the slave side" << std::endl;
  else
    std::cout << "The slave side is discretized by more elements than the master side" << std::endl;
}

/*-------------------------------------------------------*/
/*  Get function for the slave and master dof row map    */
/*                                           wirtz 02/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::DofRowMaps()
{
  // slave dof rowmap
  gsdofrowmap_ = adaptermeshtying_->SlaveDofMap();

  // master dof rowmap
  gmdofrowmap_ = adaptermeshtying_->MasterDofMap();
}

/*-------------------------------------------------------*/
/*  Get function for the P matrix            wirtz 02/16 */
/*                                                       */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ALE::Meshtying::GetMortarTrafo()
{
  return adaptermeshtying_->GetMortarTrafo();
}

/*-------------------------------------------------------*/
/*  Condensation operation block matrix      wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshtying::CondensationOperationBlockMatrix(Teuchos::RCP<LINALG::SparseOperator>& sysmat,
    Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  SplitVector(residual, splitres);

  /**********************************************************************/
  /* Condensate blockmatrix                                             */
  /**********************************************************************/

  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Condensation Operation");

  // cast Teuchos::RCP<LINALG::SparseOperator> to a Teuchos::RCP<LINALG::BlockSparseMatrixBase>
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

  /**********************************************************************/
  /* Build the final sysmat and residual                                */
  /**********************************************************************/

  // only the blocks nm, mn and mm are modified
  // the other blocks remain unchanged, since only the 2x2 block matrix system is solved
  // ---------------------        ------------------
  // | nn | nm | ns | 0  |        | nn  | nm' | ns  |
  // | mn | mm | ms | D  |   ->   | mn' | mm' | ms  |
  // | sn | sm | ss | -M |        | sn  | sm  | ss  |
  // |  0 | DT |-MT | 0  |        ------------------
  // ---------------------
  // solved system (2x2 matrix)
  // -------------
  // | nn  | nm' |
  // | mn' | mm' |
  // -------------

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process

  Teuchos::RCP<Epetra_Vector> dcnm = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dcmm = Teuchos::null;
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdcmaster(3);

  if (dconmaster_ == true and firstnonliniter_ == true)
  {
    dcnm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
    dcmm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));

    SplitVector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<LINALG::SparseMatrix> P = GetMortarTrafo();

  /*--------------------------------------------------------------------*/
  // block nm
  /*--------------------------------------------------------------------*/
  // compute modification for block nm
  Teuchos::RCP<LINALG::SparseMatrix> knm_mod =
      MLMultiply(sysmatnew->Matrix(0, 2), false, *P, false, false, false, true);

  // Add transformation matrix to nm
  sysmatnew->Matrix(0, 1).UnComplete();
  sysmatnew->Matrix(0, 1).Add(*knm_mod, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    knm_mod->Multiply(false, *(splitdcmaster[1]), *dcnm);

  /*--------------------------------------------------------------------*/
  // block mn
  /*--------------------------------------------------------------------*/
  // compute modification for block kmn
  Teuchos::RCP<LINALG::SparseMatrix> kmn_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 0), false, false, false, true);

  // Add transformation matrix to mn
  sysmatnew->Matrix(1, 0).UnComplete();
  sysmatnew->Matrix(1, 0).Add(*kmn_mod, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // block mm
  /*--------------------------------------------------------------------*/
  // compute modification for block kmm
  Teuchos::RCP<LINALG::SparseMatrix> kss_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<LINALG::SparseMatrix> kmm_mod =
      MLMultiply(*kss_mod, false, *P, false, false, false, true);

  // Add transformation matrix to mm
  sysmatnew->Matrix(1, 1).UnComplete();
  sysmatnew->Matrix(1, 1).Add(*kmm_mod, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  // complete matrix
  sysmatnew->Complete();

  //*************************************************
  //  condensation operation for the residual
  //*************************************************

  // r_m: add P^T*r_s
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  P->Multiply(true, *(splitres[2]), *fm_mod);

  // r_m: insert Dirichlet boundary conditions
  if (dconmaster_ == true and firstnonliniter_ == true) fm_mod->Update(-1.0, *dcmm, 1.0);

  // export and add r_m subvector to residual
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod, *fm_modexp);
  residual->Update(1.0, *fm_modexp, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
  {
    Teuchos::RCP<Epetra_Vector> fn_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
    LINALG::Export(*dcnm, *fn_exp);
    residual->Update(-1.0, *fn_exp, 1.0);
  }

  // export r_s = zero to residual
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  LINALG::Export(*fs_mod, *residual);

  return;
}

/*-------------------------------------------------------*/
/*  Compute and update Slave DOF's           wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshtying::UpdateSlaveDOF(
    Teuchos::RCP<Epetra_Vector>& inc, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.4)   - Update slave DOF");

  // get dof row map
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // split incremental and displacement vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitinc(3);
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdisp(3);
  SplitVector(inc, splitinc);
  SplitVector(dispnp, splitdisp);

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process

  // split vector containing Dirichlet boundary conditions, if any
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdcmaster(3);
  if (dconmaster_ == true and firstnonliniter_ == true) SplitVector(valuesdc_, splitdcmaster);

  // get transformation matrix
  Teuchos::RCP<LINALG::SparseMatrix> P = GetMortarTrafo();

  // define new incremental vector
  Teuchos::RCP<Epetra_Vector> incnew = LINALG::CreateVector(*dofrowmap, true);

  // delta_vp^s: add P*delta_vp^m
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  P->Multiply(false, *(splitinc[1]), *fs_mod);

  // delta_vp^s: subtract vp_i^s
  fs_mod->Update(-1.0, *(splitdisp[2]), 1.0);

  // delta_vp^s: add P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fs_mod_m = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  P->Multiply(false, *(splitdisp[1]), *fs_mod_m);
  fs_mod->Update(1.0, *fs_mod_m, 1.0);

  // set Dirichlet boundary conditions, if any
  if (dconmaster_ == true and firstnonliniter_ == true)
  {
    Teuchos::RCP<Epetra_Vector> fsdc_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
    P->Multiply(false, *(splitdcmaster[1]), *fsdc_mod);
    fs_mod->Update(1.0, *fsdc_mod, 1.0);
  }

  // export interior degrees of freedom
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitinc[0]), *fnexp);
  incnew->Update(1.0, *fnexp, 1.0);

  // export master degrees of freedom
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitinc[1]), *fmexp);
  incnew->Update(1.0, *fmexp, 1.0);

  // export slave degrees of freedom
  Teuchos::RCP<Epetra_Vector> fs_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*fs_mod, *fs_modexp);
  incnew->Update(1.0, *fs_modexp, 1.0);

  // set iteration counter for Dirichlet boundary conditions, if any
  if (dconmaster_ == true and firstnonliniter_ == true) firstnonliniter_ = false;

  // define incremental vector to new incremental vector
  inc = incnew;

  return;
}

/*-------------------------------------------------------*/
/*  solve mesh-tying system                  wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
int ALE::Meshtying::SolveMeshtying(LINALG::Solver& solver,
    Teuchos::RCP<LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector>& disi,
    Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3)   Solve meshtying system");

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatsolve =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmatsolve_);

  Teuchos::RCP<Epetra_Vector> res = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dis = Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::null;

  res = LINALG::CreateVector(*mergedmap_, true);
  dis = LINALG::CreateVector(*mergedmap_, true);

  mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap_, 108, false, true));

  int errorcode = 0;

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");
    SplitVectorBasedOn3x3(residual, res);

    // assign blocks to the solution matrix
    sysmatsolve->Assign(0, 0, LINALG::View, sysmatnew->Matrix(0, 0));
    sysmatsolve->Assign(0, 1, LINALG::View, sysmatnew->Matrix(0, 1));
    sysmatsolve->Assign(1, 0, LINALG::View, sysmatnew->Matrix(1, 0));
    sysmatsolve->Assign(1, 1, LINALG::View, sysmatnew->Matrix(1, 1));
    sysmatsolve->Complete();

    mergedmatrix = sysmatsolve->Merge();
  }

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");

    errorcode = solver_.Solve(mergedmatrix->EpetraOperator(), dis, res, true);

    LINALG::Export(*dis, *disi);
    LINALG::Export(*res, *residual);
    // compute and update slave dof's
    UpdateSlaveDOF(disi, dispnp);
  }
  return errorcode;
}
