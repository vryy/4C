/*--------------------------------------------------------------------------*/
/*! \file

\brief Mesh tying for ale problems

\level 2

*/
/*--------------------------------------------------------------------------*/


#define DIRECTMANIPULATION
#define ZEROSYSMAT

#include "4C_ale_meshtying.hpp"

#include "4C_ale_utils.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_mortar_node.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


ALE::Meshtying::Meshtying(Teuchos::RCP<Discret::Discretization> dis, Core::LinAlg::Solver& solver,
    int msht, int nsd, const UTILS::MapExtractor* surfacesplitter)
    : discret_(dis),
      solver_(solver),
      dofrowmap_(discret_->dof_row_map()),
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
Teuchos::RCP<Core::LinAlg::SparseOperator> ALE::Meshtying::Setup(
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
  gsmdofrowmap_ = Core::LinAlg::MergeMap(*gmdofrowmap_, *gsdofrowmap_, false);

  // dofrowmap for discretisation without slave and master dofrowmap
  gndofrowmap_ = Core::LinAlg::SplitMap(*dofrowmap_, *gsmdofrowmap_);

  // map for 2x2 (uncoupled dof's & master dof's)
  mergedmap_ = Core::LinAlg::MergeMap(*gndofrowmap_, *gmdofrowmap_, false);

  // std::cout << "number of n dof   " << gndofrowmap_->NumGlobalElements() << std::endl;
  // std::cout << "number of m dof   " << gmdofrowmap_->NumGlobalElements() << std::endl;
  // std::cout << "number of s dof   " << gsdofrowmap_->NumGlobalElements() << std::endl;

  // generate map for blockmatrix
  std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
  alemaps.push_back(gndofrowmap_);
  alemaps.push_back(gmdofrowmap_);
  alemaps.push_back(gsdofrowmap_);

  Core::LinAlg::MultiMapExtractor extractor;

  extractor.Setup(*dofrowmap_, alemaps);

  // check, if extractor maps are valid
  extractor.check_for_valid_map_extractor();

  // allocate 3x3 block sparse matrix with the interface split strategy
  // the interface split strategy speeds up the assembling process,
  // since the information, which nodes are part of the interface, is available
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
  mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
      extractor, extractor, 108, false, true));
  // nodes on the interface
  Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->conditioned_element_map(*discret_);

  mat->SetCondElements(condelements);

  // Important: right way to do it (Tobias W.)
  // allocate 2x2 solution matrix with the default block matrix strategy in order to solve the
  // reduced system memory is not allocated(1), since the matrix gets a Teuchos::RCP on the
  // respective blocks of the 3x3 block matrix
  // ---------------
  // | knn  | knm' |
  // | kmn' | kmm' |
  // ---------------

  Core::LinAlg::MapExtractor rowmapext(*mergedmap_, gmdofrowmap_, gndofrowmap_);
  Core::LinAlg::MapExtractor dommapext(*mergedmap_, gmdofrowmap_, gndofrowmap_);
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> matsolve =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          dommapext, rowmapext, 1, false, true));
  sysmatsolve_ = matsolve;

  return mat;
}

/*-------------------------------------------------------*/
/*  Use the split of the ale mesh tying for the sysmat   */
/*                                           wirtz 01/16 */
/*-------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseOperator> ALE::Meshtying::MshtSplit()
{
  // generate map for blockmatrix
  std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
  alemaps.push_back(gndofrowmap_);
  alemaps.push_back(gmdofrowmap_);
  alemaps.push_back(gsdofrowmap_);

  Core::LinAlg::MultiMapExtractor extractor;

  extractor.Setup(*dofrowmap_, alemaps);

  // check, if extractor maps are valid
  extractor.check_for_valid_map_extractor();

  // allocate 3x3 block sparse matrix with the interface split strategy
  // the interface split strategy speeds up the assembling process,
  // since the information, which nodes are part of the interface, is available
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
  mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
      extractor, extractor, 108, false, true));
  // nodes on the interface
  Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->conditioned_element_map(*discret_);
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
  // (b)  DC are projected from the master to the slave side during prepare_time_step
  //      (in project_master_to_slave_for_overlapping_bc()) (-> disabled)
  //      -> DC also influence slave nodes which are not part of the inflow
  //
  //      if(msht_ != Inpar::ALE::no_meshtying)
  //        meshtying_->project_master_to_slave_for_overlapping_bc(dispnp_, dbcmaps_->CondMap());
  //
  // (c)  DC are included in the condensation process (-> actual strategy)

  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(bmaps);
  Teuchos::RCP<const Epetra_Map> gmdofrowmap = gmdofrowmap_;
  intersectionmaps.push_back(gmdofrowmap);
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(intersectionmaps);

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
void ALE::Meshtying::prepare_meshtying_system(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation block matrix");

  condensation_operation_block_matrix(sysmat, residual, dispnp);
  return;
}

/*-------------------------------------------------------*/
/*  Split Vector                             wirtz 01/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::split_vector(
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
  Core::LinAlg::split_vector(*dofrowmap_, *vector, gsmdofrowmap_, fsm, gndofrowmap_, fn);

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  Core::LinAlg::split_vector(*gsmdofrowmap_, *fsm, gsdofrowmap_, fs, gmdofrowmap_, fm);

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
void ALE::Meshtying::split_vector_based_on3x3(
    Teuchos::RCP<Epetra_Vector> orgvector, Teuchos::RCP<Epetra_Vector> vectorbasedon2x2)
{
  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvector(3);

  split_vector(orgvector, splitvector);
  // build up the reduced residual
  Core::LinAlg::Export(*(splitvector[0]), *vectorbasedon2x2);
  Core::LinAlg::Export(*(splitvector[1]), *vectorbasedon2x2);

  return;
}



/*-------------------------------------------------------*/
/*  Set the flag for multifield problems     wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshtying::IsMultifield(
    const Core::LinAlg::MultiMapExtractor& interface,  ///< interface maps for split of ale matrix
    bool ismultifield                                  ///< flag for multifield problems
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
void ALE::Meshtying::MshtSplit(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat)
{
  if (is_multifield_)
  {
    // generate map for blockmatrix
    std::vector<Teuchos::RCP<const Epetra_Map>> alemaps;
    alemaps.push_back(gndofrowmap_);
    alemaps.push_back(gmdofrowmap_);
    alemaps.push_back(gsdofrowmap_);

    Core::LinAlg::MultiMapExtractor extractor;

    extractor.Setup(*dofrowmap_, alemaps);

    // check, if extractor maps are valid
    extractor.check_for_valid_map_extractor();

    // allocate 3x3 block sparse matrix with the interface split strategy
    // the interface split strategy speeds up the assembling process,
    // since the information, which nodes are part of the interface, is available
    // -------------------
    // | knn | knm | kns |
    // | kmn | kmm | kms |
    // | ksn | ksm | kss |
    // -------------------

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat;
    mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>(
        extractor, extractor, 108, false, true));
    // nodes on the interface
    Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->conditioned_element_map(*discret_);
    mat->SetCondElements(condelements);

    sysmat = mat;
  }
}

/*-------------------------------------------------------*/
/*  Use the split of the multifield problem for the      */
/*  sysmat                                   wirtz 01/16 */
/*-------------------------------------------------------*/
void ALE::Meshtying::MultifieldSplit(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat)
{
  if (is_multifield_)
  {
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);

    Teuchos::RCP<Epetra_Vector> ones =
        Teuchos::rcp(new Epetra_Vector(sysmatnew->Matrix(2, 2).RowMap()));
    ones->PutScalar(1.0);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
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

    Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix = sysmatnew->Merge();

    Teuchos::RCP<Core::LinAlg::MapExtractor> extractor =
        Teuchos::rcp(new Core::LinAlg::MapExtractor(
            *multifield_interface_.FullMap(), multifield_interface_.Map(1)));

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<ALE::UTILS::InterfaceSplitStrategy>> mat =
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
  adaptermeshtying_ = Teuchos::rcp(new Core::Adapter::CouplingMortar(
      Global::Problem::Instance()->NDim(), Global::Problem::Instance()->mortar_coupling_params(),
      Global::Problem::Instance()->contact_dynamic_params(),
      Global::Problem::Instance()->spatial_approximation_type()));

  // Setup of meshtying adapter
  adaptermeshtying_->Setup(discret_, discret_, Teuchos::null, coupleddof, "Mortar",
      discret_->Comm(), Global::Problem::Instance()->FunctionManager(), true);
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
Teuchos::RCP<Core::LinAlg::SparseMatrix> ALE::Meshtying::GetMortarMatrixP()
{
  return adaptermeshtying_->GetMortarMatrixP();
}

/*-------------------------------------------------------*/
/*  Condensation operation block matrix      wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void ALE::Meshtying::condensation_operation_block_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat, Teuchos::RCP<Epetra_Vector>& residual,
    Teuchos::RCP<Epetra_Vector>& dispnp)
{
  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  split_vector(residual, splitres);

  /**********************************************************************/
  /* Condensate blockmatrix                                             */
  /**********************************************************************/

  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Condensation Operation");

  // cast Teuchos::RCP<Core::LinAlg::SparseOperator> to a
  // Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);

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

    split_vector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = GetMortarMatrixP();

  /*--------------------------------------------------------------------*/
  // block nm
  /*--------------------------------------------------------------------*/
  // compute modification for block nm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_mod =
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
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 0), false, false, false, true);

  // Add transformation matrix to mn
  sysmatnew->Matrix(1, 0).UnComplete();
  sysmatnew->Matrix(1, 0).Add(*kmn_mod, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // block mm
  /*--------------------------------------------------------------------*/
  // compute modification for block kmm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kss_mod =
      MLMultiply(*P, true, sysmatnew->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_mod =
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
  Core::LinAlg::Export(*fm_mod, *fm_modexp);
  residual->Update(1.0, *fm_modexp, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
  {
    Teuchos::RCP<Epetra_Vector> fn_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
    Core::LinAlg::Export(*dcnm, *fn_exp);
    residual->Update(-1.0, *fn_exp, 1.0);
  }

  // export r_s = zero to residual
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  Core::LinAlg::Export(*fs_mod, *residual);

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
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // split incremental and displacement vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitinc(3);
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdisp(3);
  split_vector(inc, splitinc);
  split_vector(dispnp, splitdisp);

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process

  // split vector containing Dirichlet boundary conditions, if any
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdcmaster(3);
  if (dconmaster_ == true and firstnonliniter_ == true) split_vector(valuesdc_, splitdcmaster);

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = GetMortarMatrixP();

  // define new incremental vector
  Teuchos::RCP<Epetra_Vector> incnew = Core::LinAlg::CreateVector(*dofrowmap, true);

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
  Core::LinAlg::Export(*(splitinc[0]), *fnexp);
  incnew->Update(1.0, *fnexp, 1.0);

  // export master degrees of freedom
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  Core::LinAlg::Export(*(splitinc[1]), *fmexp);
  incnew->Update(1.0, *fmexp, 1.0);

  // export slave degrees of freedom
  Teuchos::RCP<Epetra_Vector> fs_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  Core::LinAlg::Export(*fs_mod, *fs_modexp);
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
int ALE::Meshtying::SolveMeshtying(Core::LinAlg::Solver& solver,
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector>& disi,
    Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3)   Solve meshtying system");

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatsolve =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmatsolve_);

  Teuchos::RCP<Epetra_Vector> res = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dis = Teuchos::null;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix = Teuchos::null;

  res = Core::LinAlg::CreateVector(*mergedmap_, true);
  dis = Core::LinAlg::CreateVector(*mergedmap_, true);

  mergedmatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mergedmap_, 108, false, true));

  int errorcode = 0;

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");
    split_vector_based_on3x3(residual, res);

    // assign blocks to the solution matrix
    sysmatsolve->Assign(0, 0, Core::LinAlg::View, sysmatnew->Matrix(0, 0));
    sysmatsolve->Assign(0, 1, Core::LinAlg::View, sysmatnew->Matrix(0, 1));
    sysmatsolve->Assign(1, 0, Core::LinAlg::View, sysmatnew->Matrix(1, 0));
    sysmatsolve->Assign(1, 1, Core::LinAlg::View, sysmatnew->Matrix(1, 1));
    sysmatsolve->Complete();

    mergedmatrix = sysmatsolve->Merge();
  }

  {
    TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");

    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    errorcode = solver_.Solve(mergedmatrix->EpetraOperator(), dis, res, solver_params);

    Core::LinAlg::Export(*dis, *disi);
    Core::LinAlg::Export(*res, *residual);
    // compute and update slave dof's
    UpdateSlaveDOF(disi, dispnp);
  }
  return errorcode;
}

FOUR_C_NAMESPACE_CLOSE
