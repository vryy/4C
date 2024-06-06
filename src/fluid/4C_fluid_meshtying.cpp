/*-----------------------------------------------------------*/
/*! \file

\brief Methods to apply meshtying to fluid and scatra systems


\level 2

*/
/*-----------------------------------------------------------*/

#define DIRECTMANIPULATION
#define ZEROSYSMAT

#include "4C_fluid_meshtying.hpp"

#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_mortar_node.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


FLD::Meshtying::Meshtying(Teuchos::RCP<Discret::Discretization> dis, Core::LinAlg::Solver& solver,
    int msht, int nsd, const UTILS::MapExtractor* surfacesplitter)
    : discret_(dis),
      solver_(solver),
      msht_(msht),
      myrank_(dis->Comm().MyPID()),
      surfacesplitter_(surfacesplitter),
      dofrowmap_(nullptr),
      problemrowmap_(Teuchos::null),
      gndofrowmap_(Teuchos::null),
      gsmdofrowmap_(Teuchos::null),
      gsdofrowmap_(Teuchos::null),
      gmdofrowmap_(Teuchos::null),
      mergedmap_(Teuchos::null),
      valuesdc_(Teuchos::null),
      adaptermeshtying_(
          Teuchos::rcp(new Core::Adapter::CouplingMortar(Global::Problem::Instance()->NDim(),
              Global::Problem::Instance()->mortar_coupling_params(),
              Global::Problem::Instance()->contact_dynamic_params(),
              Global::Problem::Instance()->spatial_approximation_type()))),
      pcoupled_(true),
      dconmaster_(false),
      firstnonliniter_(false),
      nsd_(nsd),
      multifield_condelements_(Teuchos::null),
      multifield_condelements_shape_(Teuchos::null),
      multifield_splitmatrix_(false),
      is_multifield_(false)
{
}

/*-------------------------------------------------------*/
/*  Setup mesh-tying problem                ehrl (04/11) */
/*-------------------------------------------------------*/

void FLD::Meshtying::setup_meshtying(const std::vector<int>& coupleddof, const bool pcoupled)
{
  // get pointer to dof row map
  dofrowmap_ = discret_->dof_row_map();

  // set whether pressure dof is coupled
  pcoupled_ = pcoupled;

  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  1)   Setup Meshtying");

  // Setup of meshtying adapter
  adaptermeshtying_->Setup(discret_, discret_, Teuchos::null, coupleddof, "Mortar",
      discret_->Comm(), Global::Problem::Instance()->FunctionManager(), true);

  // 4 different systems to solve
  // a) Condensation with a block matrix (condensed_bmat)
  //    system is solved in a 2x2 (n,m) block matrix with the respective solvers

  // b) Condensation with a block matrix merged to a sparse martix (condensed_bmat_merged)
  //    - condensation operation is done in the 2x2 (n,m) block matrix (no splitting operations)
  //      -> graph can be saved resulting in accelerated element assembly
  //         (ifdef: allocation of new matrix, more memory, slower element assembly,
  //                 no block matrix subtraction)
  //    - one's are assigned to the diagonal entries in the ss-block (as for dirichlet conditions)
  //    properties:
  //    - splitting operation is an additional, time consuming operation
  //    - resulting system matrix is easy to solve

  // c) Condensation in a sparse matrix (condensed_smat)
  //    - condensation operation is done in the original 3x3 (n,m,s) sparse matrix
  //      by splitting operations
  //      -> graph can be saved resulting in accelerated element assembly
  //         (ifdef: allocation of new matrix, more memory, slower element assembly,
  //                 no block matrix subtraction)
  //    - one's are assigned to the diagonal entries in the ss-block (as for dirichlet conditions)
  //    properties:
  //    - splitting operation is an additional, time consuming operation
  //    - resulting system matrix is easy to solve

  // these options were deleted (ehrl 19.06.2013)
  // since the implementation was only temporary and not well tested
  // c) Saddle point system sparse matrix (sps_coupled)
  // d) Saddle point system block matrix (sps_pc)

  // number of nodes master < number of nodes slave
  // -> better results, Krylov does not work ??

  if (myrank_ == 0)
  {
    int numdofmaster = (adaptermeshtying_->MasterDofMap())->NumGlobalElements();
    int numdofslave = (adaptermeshtying_->SlaveDofMap())->NumGlobalElements();

    std::cout << std::endl << "number of master dof's:   " << numdofmaster << std::endl;
    std::cout << "number of slave dof's:   " << numdofslave << std::endl << std::endl;

    if (numdofmaster > numdofslave)
    {
      std::cout << "The master side is discretized by more elements than the slave side"
                << std::endl;
    }
    else
      std::cout << "The slave side is discretized by more elements than the master side"
                << std::endl;
  }

  switch (msht_)
  {
    case Inpar::FLUID::condensed_bmat:
    case Inpar::FLUID::condensed_bmat_merged:
    {
      if (!pcoupled_)
      {
        FOUR_C_THROW(
            "The system cannot be solved in a block matrix!! \n"
            "The null space does not have the right length. Fix it or use option Smat");
      }

      // slave dof rowmap
      gsdofrowmap_ = adaptermeshtying_->SlaveDofMap();

      // master dof rowmap
      gmdofrowmap_ = adaptermeshtying_->MasterDofMap();

      // merge dofrowmap for slave and master discretization
      gsmdofrowmap_ = Core::LinAlg::MergeMap(*gmdofrowmap_, *gsdofrowmap_, false);

      // dofrowmap for discretisation without slave and master dofrowmap
      gndofrowmap_ = Core::LinAlg::SplitMap(*dofrowmap_, *gsmdofrowmap_);

      // map for 2x2 (uncoupled dof's & master dof's)
      mergedmap_ = Core::LinAlg::MergeMap(*gndofrowmap_, *gmdofrowmap_, false);

      // std::cout << "number of n dof   " << gndofrowmap_->NumGlobalElements() << std::endl;
      // std::cout << "number of m dof   " << gmdofrowmap_->NumGlobalElements() << std::endl;
      // std::cout << "number of s dof   " << gsdofrowmap_->NumGlobalElements() << std::endl;

      // Important: right way to do it (Tobias W.)
      // allocate 2x2 solution matrix with the default block matrix strategy in order to solve the
      // reduced system memory is not allocated(1), since the matrix gets a Teuchos::RCP on the
      // respective blocks of the 3x3 block matrix
      // ---------------
      // | knn  | knm' |
      // | kmn' | kmm' |
      // ---------------

      Core::LinAlg::MapExtractor mapext(*mergedmap_, gmdofrowmap_, gndofrowmap_);
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
          matsolve = Teuchos::rcp(
              new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                  mapext, mapext, 1, false, true));
      sysmatsolve_ = matsolve;

      // fixing length of nullspace for block matrix (solver/preconditioner ML)
      if (msht_ == Inpar::FLUID::condensed_bmat_merged)
      {
        std::string inv = "BMatMerged";
        const Epetra_Map& oldmap = *(dofrowmap_);
        const Epetra_Map& newmap = *(mergedmap_);
        Core::LinearSolver::Parameters::FixNullSpace(inv.data(), oldmap, newmap, solver_.Params());
        std::cout << std::endl;
      }
      else if (msht_ == Inpar::FLUID::condensed_bmat)
      {
        // fixing length of Inverse1 nullspace (solver/preconditioner ML)
        {
          std::string inv = "Inverse1";
          const Epetra_Map& oldmap = *(dofrowmap_);
          const Epetra_Map& newmap = matsolve->Matrix(0, 0).EpetraMatrix()->RowMap();
          Core::LinearSolver::Parameters::FixNullSpace(
              inv.data(), oldmap, newmap, solver_.Params().sublist("Inverse1"));
          std::cout << std::endl;
        }
        // fixing length of Inverse2 nullspace (solver/preconditioner ML)
        {
          std::string inv = "Inverse2";
          const Epetra_Map& oldmap = *(dofrowmap_);
          const Epetra_Map& newmap = matsolve->Matrix(1, 1).EpetraMatrix()->RowMap();
          Core::LinearSolver::Parameters::FixNullSpace(
              inv.data(), oldmap, newmap, solver_.Params().sublist("Inverse2"));
          std::cout << std::endl;
        }
      }
    }
    break;
    case Inpar::FLUID::condensed_smat:
    {
      // slave dof rowmap
      gsdofrowmap_ = adaptermeshtying_->SlaveDofMap();

      // master dof rowmap
      gmdofrowmap_ = adaptermeshtying_->MasterDofMap();

      // merge dofrowmap for slave and master discretization
      gsmdofrowmap_ = Core::LinAlg::MergeMap(*gmdofrowmap_, *gsdofrowmap_, false);

      // dofrowmap for discretisation without slave and master dofrowmap
      gndofrowmap_ = Core::LinAlg::SplitMap(*dofrowmap_, *gsmdofrowmap_);

      if (myrank_ == 0)
      {
#ifdef ZEROSYSMAT
        std::cout << "Condensation operation takes place in the original sysmat -> graph is saved"
                  << std::endl;
        std::cout << "The sysmat is set to zero and all parts are added -> exact" << std::endl
                  << std::endl;
#else

#ifdef DIRECTMANIPULATION
        std::cout << "Condensation operation takes place in the original sysmat -> graph is saved"
                  << std::endl
                  << std::endl;
#else

        std::cout << "Condensation operation is carried out in a new allocated sparse matrix -> "
                     "graph is not saved"
                  << std::endl
                  << std::endl;
#endif
#endif
      }
    }
    break;
    default:
      FOUR_C_THROW("Choose a correct mesh-tying option");
      break;
  }
}

Teuchos::RCP<Core::LinAlg::SparseOperator> FLD::Meshtying::init_system_matrix() const
{
  switch (msht_)
  {
    case Inpar::FLUID::condensed_bmat:
    case Inpar::FLUID::condensed_bmat_merged:
    {
      // generate map for blockmatrix
      std::vector<Teuchos::RCP<const Epetra_Map>> fluidmaps;
      fluidmaps.push_back(gndofrowmap_);
      fluidmaps.push_back(gmdofrowmap_);
      fluidmaps.push_back(gsdofrowmap_);

      Core::LinAlg::MultiMapExtractor extractor;

      extractor.Setup(*dofrowmap_, fluidmaps);

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

      Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>> mat;
      mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(
          extractor, extractor, 108, false, true));
      // nodes on the interface
      Teuchos::RCP<std::set<int>> condelements =
          surfacesplitter_->conditioned_element_map(*discret_);
      mat->SetCondElements(condelements);

      return mat;
    }
    case Inpar::FLUID::condensed_smat:
    {
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_, 108, false, true));
    }
    default:
    {
      FOUR_C_THROW("Choose a correct mesh-tying option");
      break;
    }
  }
  return Teuchos::null;
}

const Epetra_Map* FLD::Meshtying::GetMergedMap()
{
  const Epetra_Map& newmap = *(mergedmap_);

  return &newmap;
}

/*-----------------------------------------------------*/
/*  Check if there are overlapping BCs    ehrl (08/13) */
/*-----------------------------------------------------*/
void FLD::Meshtying::CheckOverlappingBC(Teuchos::RCP<Epetra_Map> map)
{
  bool overlap = false;

  // loop over all slave row nodes of the interface
  for (int j = 0; j < adaptermeshtying_->Interface()->SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = adaptermeshtying_->Interface()->SlaveRowNodes()->GID(j);
    Core::Nodes::Node* node = adaptermeshtying_->Interface()->Discret().gNode(gid);
    if (!node) FOUR_C_THROW("ERROR: Cannot find node with gid %", gid);
    auto* mtnode = static_cast<Mortar::Node*>(node);

    // check if this node's dofs are in given map
    for (int k = 0; k < mtnode->NumDof(); ++k)
    {
      int currdof = mtnode->Dofs()[k];
      int lid = map->LID(currdof);

      // found slave node intersecting with given map
      if (lid >= 0)
      {
        overlap = true;
        break;
      }
    }
  }

  // print warning message to screen
  if (overlap && myrank_ == 0)
  {
    if (myrank_ == 0)
    {
      FOUR_C_THROW(
          "Slave boundary and volume flow rate boundary conditions overlap!\n"
          "This leads to an over-constraint problem setup");
    }
  }
}

/*---------------------------------------------------*/
/*  Correct slave Dirichlet value       ehrl (12/12) */
/*---------------------------------------------------*/
void FLD::Meshtying::project_master_to_slave_for_overlapping_bc(
    Teuchos::RCP<Epetra_Vector>& velnp, Teuchos::RCP<const Epetra_Map> bmaps)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(bmaps);
  Teuchos::RCP<const Epetra_Map> gmdofrowmap = gmdofrowmap_;
  intersectionmaps.push_back(gmdofrowmap);
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    if (myrank_ == 0)
    {
      std::cout << "number of intersecting elements:  " << intersectionmap->NumGlobalElements()
                << std::endl;
      std::cout << "Overlapping boundary conditions are projected from the master to the slave side"
                << std::endl
                << std::endl;
    }

    // All dofs of the master are projected to the slave nodes. Therefore, all values of slave nodes
    // are overwritten although it is not necessary for the nodes which do not intersect with the
    // overlapping BC. In the case of Dirichlet or Dirichlet-like BC (Wormersly) new values are
    // assigned to the master side of internal interface but not to the slave side.

    // std::cout << "BEFORE projection from master to slave side" << std::endl;
    // OutputVectorSplit(velnp);

    UpdateSlaveDOF(velnp, velnp);

    // std::cout << "AFTER projection from master to slave side" << std::endl;
    // OutputVectorSplit(velnp);
  }
}

/*------------------------------------------------------------------------------*/
/*  Check if Dirichlet BC are defined on the master                ehrl (08/13) */
/*------------------------------------------------------------------------------*/
void FLD::Meshtying::DirichletOnMaster(Teuchos::RCP<const Epetra_Map> bmaps)
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
  //      if(msht_ != Inpar::FLUID::no_meshtying)
  //        meshtying_->project_master_to_slave_for_overlapping_bc(velnp_, dbcmaps_->CondMap());
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
}

/*------------------------------------------------------------------*/
/*  Include Dirichlet BC in condensation operation     ehrl (08/13) */
/*------------------------------------------------------------------*/
void FLD::Meshtying::include_dirichlet_in_condensation(
    const Teuchos::RCP<Epetra_Vector>& velnp, const Teuchos::RCP<Epetra_Vector>& veln)
{
  if (dconmaster_)
  {
    valuesdc_ = Core::LinAlg::CreateVector(*dofrowmap_, true);
    valuesdc_->Update(1.0, *velnp, 1.0);
    valuesdc_->Update(-1.0, *veln, 1.0);

    firstnonliniter_ = true;
  }
}


/*---------------------------------------------------*/
/*  evaluation of matrix P with potential            */
/*  mesh relocation in ALE case             vg 01/14 */
/*---------------------------------------------------*/
void FLD::Meshtying::evaluate_with_mesh_relocation(Teuchos::RCP<Epetra_Vector>& dispnp)
{
  // get ALE discretization
  Teuchos::RCP<Discret::Discretization> aledis = Global::Problem::Instance()->GetDis("ale");

  // call mortar evaluate routine including mesh correction
  adaptermeshtying_->evaluate_with_mesh_relocation(
      discret_, aledis, dispnp, discret_->Comm(), true);
}

/*---------------------------------------------------*/
/*  Prepare Meshtying                    wirtz 02/16 */
/*---------------------------------------------------*/
void FLD::Meshtying::PrepareMeshtying(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& residual, const Teuchos::RCP<Epetra_Vector>& velnp,
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& shapederivatives)
{
  prepare_meshtying_system(sysmat, residual, velnp);
  MultifieldSplit(sysmat);

  if (shapederivatives != Teuchos::null)
  {
    condensation_operation_block_matrix_shape(shapederivatives);
    multifield_split_shape(shapederivatives);
  }
}


/*---------------------------------------------------*/
/*  Prepare Meshtying system            ehrl (04/11) */
/*---------------------------------------------------*/
void FLD::Meshtying::prepare_meshtying_system(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& residual, const Teuchos::RCP<Epetra_Vector>& velnp)
{
  switch (msht_)
  {
    case Inpar::FLUID::condensed_bmat:
    case Inpar::FLUID::condensed_bmat_merged:
      condensation_block_matrix(sysmat, residual, velnp);
      break;
    case Inpar::FLUID::condensed_smat:
      condensation_sparse_matrix(sysmat, residual, velnp);
      break;
    default:
      FOUR_C_THROW("Meshtying algorithm not recognized!");
      break;
  }
}


/*---------------------------------------------------*/
/*---------------------------------------------------*/
void FLD::Meshtying::ApplyPTToResidual(Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat,
    Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Core::LinAlg::KrylovProjector> projector)
{
  // define residual vector for case of block matrix
  Teuchos::RCP<Epetra_Vector> res = Core::LinAlg::CreateVector(*mergedmap_, true);

  // split original residual vector
  split_vector_based_on3x3(residual, res);

  // apply projector
  projector->ApplyPT(*res);

  // export residual back to original vector
  Core::LinAlg::Export(*res, *residual);
}


/*-------------------------------------------------------*/
/*  Krylov projection                      ehrl (04/11)  */
/*-------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::Meshtying::adapt_krylov_projector(Teuchos::RCP<Epetra_Vector> vec)
{
  if (pcoupled_)
  {
    // Remove slave nodes from vec
    Teuchos::RCP<Epetra_Vector> fm_slave = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
    // add fm subvector to feffnew
    Core::LinAlg::Export(*fm_slave, *vec);

    switch (msht_)
    {
      case Inpar::FLUID::condensed_bmat:
      case Inpar::FLUID::condensed_bmat_merged:
      {
        Teuchos::RCP<Epetra_Vector> vec_mesht = Core::LinAlg::CreateVector(*mergedmap_, true);
        split_vector_based_on3x3(vec, vec_mesht);
        vec = vec_mesht;
      }
      break;
      case Inpar::FLUID::condensed_smat:
        break;
      default:
        FOUR_C_THROW("Krylov projection not supported for this meshtying option.");
        break;
    }
  }
  return vec;
}

/*-------------------------------------------------------*/
/*  solve mesh-tying system                ehrl (04/11)  */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::SolveMeshtying(Core::LinAlg::Solver& solver,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& incvel, const Teuchos::RCP<Epetra_Vector>& residual,
    const Teuchos::RCP<Epetra_Vector>& velnp, const int itnum,
    Core::LinAlg::SolverParams& solver_params)
{
  solver_params.refactor = true;
  solver_params.reset = itnum == 1;

  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3)   Solve meshtying system");

  switch (msht_)
  {
    case Inpar::FLUID::condensed_bmat:
    {
      Teuchos::RCP<Epetra_Vector> res = Core::LinAlg::CreateVector(*mergedmap_, true);
      Teuchos::RCP<Epetra_Vector> inc = Core::LinAlg::CreateVector(*mergedmap_, true);

      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatsolve =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmatsolve_);

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");

        split_vector_based_on3x3(residual, res);
        // assign blocks to the solution matrix
        sysmatsolve->Assign(0, 0, Core::LinAlg::View, sysmatnew->Matrix(0, 0));
        sysmatsolve->Assign(0, 1, Core::LinAlg::View, sysmatnew->Matrix(0, 1));
        sysmatsolve->Assign(1, 0, Core::LinAlg::View, sysmatnew->Matrix(1, 0));
        sysmatsolve->Assign(1, 1, Core::LinAlg::View, sysmatnew->Matrix(1, 1));
        sysmatsolve->Complete();
      }

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");
        solver_.Solve(sysmatsolve->EpetraOperator(), inc, res, solver_params);
      }

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Update");

        // Export the computed increment to the global increment
        Core::LinAlg::Export(*inc, *incvel);
        Core::LinAlg::Export(*res, *residual);

        // compute and update slave dof's
        UpdateSlaveDOF(incvel, velnp);
      }
    }
    break;
    case Inpar::FLUID::condensed_bmat_merged:
    {
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatnew =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat);
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmatsolve =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmatsolve_);

      Teuchos::RCP<Epetra_Vector> res = Teuchos::null;
      Teuchos::RCP<Epetra_Vector> inc = Teuchos::null;

      Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix = Teuchos::null;

      res = Core::LinAlg::CreateVector(*mergedmap_, true);
      inc = Core::LinAlg::CreateVector(*mergedmap_, true);

      mergedmatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mergedmap_, 108, false, true));

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

        solver_.Solve(mergedmatrix->EpetraOperator(), inc, res, solver_params);

        Core::LinAlg::Export(*inc, *incvel);
        Core::LinAlg::Export(*res, *residual);
        // compute and update slave dof's
        UpdateSlaveDOF(incvel, velnp);
      }
    }
    break;
    case Inpar::FLUID::condensed_smat:
    {
      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Solve");
        solver_.Solve(sysmat->EpetraOperator(), incvel, residual, solver_params);
      }

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Update");
        // compute and update slave dof's
        UpdateSlaveDOF(incvel, velnp);
      }
    }
    break;
    default:
      FOUR_C_THROW("");
      break;
  }
}


/*-------------------------------------------------------*/
/*  Condensation Sparse Matrix              ehrl (04/11) */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::condensation_sparse_matrix(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& residual, const Teuchos::RCP<Epetra_Vector>& velnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation sparse matrix");

  /**********************************************************************/
  /* Split sysmat and residual                                          */
  /**********************************************************************/

  // container for split matrix and vector
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> splitmatrix;
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvel(3);

  split_matrix(sysmat, splitmatrix);
  split_vector(residual, splitres);
  split_vector(velnp, splitvel);

  /**********************************************************************/
  /* Condensate sparse matrix                                           */
  /**********************************************************************/

  condensation_operation_sparse_matrix(sysmat, residual, splitmatrix, splitres, splitvel);
}


/*-------------------------------------------------------*/
/*  Condensation Block Matrix               ehrl (04/11) */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::condensation_block_matrix(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& residual, const Teuchos::RCP<Epetra_Vector>& velnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation block matrix");

  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitres(3);
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvel(3);
  split_vector(residual, splitres);
  split_vector(velnp, splitvel);

  /**********************************************************************/
  /* Condensate blockmatrix                                             */
  /**********************************************************************/

  condensation_operation_block_matrix(sysmat, residual, splitres, splitvel);
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | split sparse global system matrix into 3x3 block sparse matrix associated with interior, master,
 and slave dofs   fang 08/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
void FLD::Meshtying::split_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator>
        matrix,  //!< original sparse global system matrix before split
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>&
        splitmatrix  //!< resulting block sparse matrix after split
)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Split Matrix");

  // initialize map extractor for matrix splitting
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidmaps;
  fluidmaps.push_back(gndofrowmap_);
  fluidmaps.push_back(gmdofrowmap_);
  fluidmaps.push_back(gsdofrowmap_);
  Core::LinAlg::MultiMapExtractor extractor(*dofrowmap_, fluidmaps);
  extractor.check_for_valid_map_extractor();

  // perform matrix splitting
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------
  splitmatrix = Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(matrix)
                    ->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(extractor, extractor);

  // finalize resulting block sparse matrix
  splitmatrix->Complete();
}

/*-------------------------------------------------------*/
/*  Split Vector                           ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::split_vector(
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
}


/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
void FLD::Meshtying::split_vector_based_on3x3(
    Teuchos::RCP<Epetra_Vector> orgvector, Teuchos::RCP<Epetra_Vector> vectorbasedon2x2)
{
  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvector(3);

  split_vector(orgvector, splitvector);
  // build up the reduced residual
  Core::LinAlg::Export(*(splitvector[0]), *vectorbasedon2x2);
  Core::LinAlg::Export(*(splitvector[1]), *vectorbasedon2x2);
}


/*-------------------------------------------------------*/
/*  Condensation operation sparse matrix    ehrl (04/11) */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::condensation_operation_sparse_matrix(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>&
        sysmat,                                   ///> sysmat established by the element routine
    const Teuchos::RCP<Epetra_Vector>& residual,  ///> residual established by the element routine
    const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>&
        splitmatrix,  ///> container with split original sysmat
    const std::vector<Teuchos::RCP<Epetra_Vector>>&
        splitres,  ///> container with split original residual
    const std::vector<Teuchos::RCP<Epetra_Vector>>&
        splitvel  ///> container with split velocity vector
)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.3)   - Condensation Operation");

  /**********************************************************************/
  /* Build the final sysmat                                             */
  /**********************************************************************/

  // ---------------------        ------------------
  // | nn | nm | ns | 0  |        | nn  | nm' | 0  |
  // | mn | mm | ms | D  |   =    | mn' | mm' | 0  |
  // | sn | sm | ss | -M |        |  0  |  0  | 1  |
  // |  0 | DT |-MT | 0  |        ------------------
  // ---------------------
  // solved system
  // ------------------
  // | nn  | nm' | 0  |
  // | mn' | mm' | 0  |
  // |  0  |  0  | 1  |
  // ------------------

  // splitmatrix
  // -------------------
  // | knn | knm | kns |
  // | kmn | kmm | kms |
  // | ksn | ksm | kss |
  // -------------------

  // DIRECTMANIPULATION:
  // the sysmat is manipulated directly with out changing the graph
  // -> subtract blocks to get zeros in the slave blocks
  // -> graph is saved -> fast element assembly
  // -> less memory is needed since everything is done with the original system matrix
  // -> is it dangerous to subtract blocks to get zeros
  //
  // not DIRECTMANIPULATION:
  // a new matrix is allocated
  // -> more memory is required and element time is slower since graph cannot be saved
  // -> there are zeros in the slave block by definition
  //
  // both methods work with a 3x3 (n,m,s) system matrix

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process
  Teuchos::RCP<Epetra_Vector> dcnm = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dcmm = Teuchos::null;
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdcmaster(3);

  if (dconmaster_ and firstnonliniter_)
  {
    dcnm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
    dcmm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));

    split_vector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = adaptermeshtying_->GetMortarMatrixP();

  /**********************************************************************/
  /* Condensation operation for the sysmat                              */
  /**********************************************************************/

  // the sysmat is manipulated directly with out changing the graph
  // (subtract blocks to get zeros in the slave blocks)

#ifdef ZEROSYSMAT
  sysmat->UnComplete();

  /*--------------------------------------------------------------------*/
  // Part nn
  /*--------------------------------------------------------------------*/
  sysmat->Add(splitmatrix->Matrix(0, 0), false, 1.0, 0.0);

  /*--------------------------------------------------------------------*/
  // Part nm
  /*--------------------------------------------------------------------*/
  // knm: add kns*P
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gndofrowmap_, 100));
  knm_mod->Add(splitmatrix->Matrix(0, 1), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_add =
      MLMultiply(splitmatrix->Matrix(0, 2), false, *P, false, false, false, true);
  knm_mod->Add(*knm_add, false, 1.0, 1.0);
  knm_mod->Complete(splitmatrix->Matrix(0, 1).DomainMap(), splitmatrix->Matrix(0, 1).RowMap());

  sysmat->Add(*knm_mod, false, 1.0, 1.0);

  if (dconmaster_ and firstnonliniter_) knm_add->Multiply(false, *(splitdcmaster[1]), *dcnm);

  /*--------------------------------------------------------------------*/
  // Part mn
  /*--------------------------------------------------------------------*/
  // kmn: add P^T*ksn
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmn_mod->Add(splitmatrix->Matrix(1, 0), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_add =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 0), false, false, false, true);
  kmn_mod->Add(*kmn_add, false, 1.0, 1.0);
  kmn_mod->Complete(splitmatrix->Matrix(1, 0).DomainMap(), splitmatrix->Matrix(1, 0).RowMap());

  sysmat->Add(*kmn_mod, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // Part mm
  /*--------------------------------------------------------------------*/
  // kms: add P^T*kss, kmm: add kms*P + kmm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmm_mod->Add(splitmatrix->Matrix(1, 1), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kms =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_add =
      MLMultiply(*kms, false, *P, false, false, false, true);
  kmm_mod->Add(*kmm_add, false, 1.0, 1.0);
  kmm_mod->Complete(splitmatrix->Matrix(1, 1).DomainMap(), splitmatrix->Matrix(1, 1).RowMap());

  sysmat->Add(*kmm_mod, false, 1.0, 1.0);

  if (dconmaster_ and firstnonliniter_) kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag;
  ones->PutScalar(1.0);
  onesdiag = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
  onesdiag->Complete();

  sysmat->Add(*onesdiag, false, 1.0, 1.0);

  sysmat->Complete();
#else
#ifdef DIRECTMANIPULATION
  sysmat->UnComplete();

  /*--------------------------------------------------------------------*/
  // Part nm
  /*--------------------------------------------------------------------*/
  // knm: add kns*P
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_add =
      MLMultiply(splitmatrix->Matrix(0, 2), false, *P, false, false, false, true);
  knm_add->Complete(splitmatrix->Matrix(0, 1).DomainMap(), splitmatrix->Matrix(0, 1).RowMap());
  sysmat->Add(*knm_add, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    knm_add->Multiply(false, *(splitdcmaster[1]), *dcnm);

  /*--------------------------------------------------------------------*/
  // Part mn
  /*--------------------------------------------------------------------*/
  // kmn: add P^T*ksn
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_add =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 0), false, false, false, true);
  kmn_add->Complete(splitmatrix->Matrix(1, 0).DomainMap(), splitmatrix->Matrix(1, 0).RowMap());
  sysmat->Add(*kmn_add, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // Part mm
  /*--------------------------------------------------------------------*/
  // kms: add P^T*kss, kmm: add kms*P + kmm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kms =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_add =
      MLMultiply(*kms, false, *P, false, false, false, true);
  kmm_mod->Add(*kmm_add, false, 1.0, 1.0);
  kmm_mod->Complete(splitmatrix->Matrix(1, 1).DomainMap(), splitmatrix->Matrix(1, 1).RowMap());

  sysmat->Add(*kmm_mod, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  // Dangerous??: Get zero in block ... by subtracting
  sysmat->Add(splitmatrix->Matrix(0, 2), false, -1.0, 1.0);
  sysmat->Add(splitmatrix->Matrix(1, 2), false, -1.0, 1.0);
  sysmat->Add(splitmatrix->Matrix(2, 0), false, -1.0, 1.0);
  sysmat->Add(splitmatrix->Matrix(2, 1), false, -1.0, 1.0);
  sysmat->Add(splitmatrix->Matrix(2, 2), false, -1.0, 1.0);

  Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag;
  // build identity matrix for slave dofs
  ones->PutScalar(1.0);
  // Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag = Teuchos::rcp(new
  // Core::LinAlg::SparseMatrix(*ones));
  onesdiag = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
  onesdiag->Complete();

  sysmat->Add(*onesdiag, false, 1.0, 1.0);

  sysmat->Complete();

  // output_sparse_matrix_split(sysmat);

#else
  // the sysmat is manipulated indirectly via a second sparse matrix
  // and therefore, the graph changes
  Teuchos::RCP<Core::LinAlg::SparseOperator> sysmatnew = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*dofrowmapsysmatnew->Matrix(0, 2) _, 81, true, false));

  /*--------------------------------------------------------------------*/
  // Part nn
  /*--------------------------------------------------------------------*/
  sysmatnew->Add(splitmatrix->Matrix(0, 0), false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // Part nm
  /*--------------------------------------------------------------------*/
  // knm: add kns*P
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gndofrowmap_, 100));
  knm_mod->Add(splitmatrix->Matrix(0, 1), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_add =
      MLMultiply(splitmatrix->Matrix(0, 2), false, *P, false, false, false, true);
  knm_mod->Add(*knm_add, false, 1.0, 1.0);
  knm_mod->Complete(splitmatrix->Matrix(0, 1).DomainMap(), splitmatrix->Matrix(0, 1).RowMap());

  sysmatnew->Add(*knm_mod, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    knm_add->Multiply(false, *(splitdcmaster[1]), *dcnm);

  /*--------------------------------------------------------------------*/
  // Part mn
  /*--------------------------------------------------------------------*/
  // kmn: add P^T*ksn
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmn_mod->Add(splitmatrix->Matrix(1, 0), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_add =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 0), false, false, false, true);
  kmn_mod->Add(*kmn_add, false, 1.0, 1.0);
  kmn_mod->Complete(splitmatrix->Matrix(1, 0).DomainMap(), splitmatrix->Matrix(1, 0).RowMap());

  sysmatnew->Add(*kmn_mod, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // Part mm
  /*--------------------------------------------------------------------*/
  // kms: add P^T*kss, kmm: add kms*P + kmm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_mod =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gmdofrowmap_, 100));
  kmm_mod->Add(splitmatrix->Matrix(1, 1), false, 1.0, 1.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kms =
      MLMultiply(*P, true, splitmatrix->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_add =
      MLMultiply(*kms, false, *P, false, false, false, true);
  kmm_mod->Add(*kmm_add, false, 1.0, 1.0);
  kmm_mod->Complete(splitmatrix->Matrix(1, 1).DomainMap(), splitmatrix->Matrix(1, 1).RowMap());

  sysmatnew->Add(*kmm_mod, false, 1.0, 1.0);

  if (dconmaster_ == true and firstnonliniter_ == true)
    kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag;
  ones->PutScalar(1.0);
  onesdiag = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
  onesdiag->Complete();

  sysmatnew->Add(*onesdiag, false, 1.0, 1.0);

  sysmatnew->Complete();

  sysmat = sysmatnew;
#endif
#endif

  //*************************************************
  //  condensation operation for the residual
  //*************************************************
  // splitres[ii]
  // r_n [0]
  // r_m [1]
  // r_s [2]

  Teuchos::RCP<Epetra_Vector> resnew = Core::LinAlg::CreateVector(*dofrowmap_, true);

  // r_m: add P^T*r_s
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  P->Multiply(true, *(splitres[2]), *fm_mod);

  // r_m: add P^T*K_ss*vp_i^s
  Teuchos::RCP<Epetra_Vector> fm_mod_ss = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  kms->Multiply(false, *(splitvel[2]), *fm_mod_ss);
  fm_mod->Update(1.0, *fm_mod_ss, 1.0);

  // r_m: subtract P^T*K_ss*P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fm_mod_mm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  kmm_add->Multiply(false, *(splitvel[1]), *fm_mod_mm);
  fm_mod->Update(-1.0, *fm_mod_mm, 1.0);

  // r_m: insert Dirichlet boundary conditions
  if (dconmaster_ and firstnonliniter_) fm_mod->Update(-1.0, *dcmm, 1.0);

  // export additions to r_m subvector to r_new
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*fm_mod, *fm_modexp);
  resnew->Update(1.0, *fm_modexp, 1.0);

  // export r_m subvector to r_new
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*(splitres[1]), *fmexp);
  resnew->Update(1.0, *fmexp, 1.0);


  // r_n: add K_ns*vp_i^s
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(splitmatrix->Matrix(0, 2)));
  // Core::LinAlg::SparseMatrix& knm = splitmatrix->Matrix(0,2);
  Teuchos::RCP<Epetra_Vector> fn_mod = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
  knm->Multiply(false, *(splitvel[2]), *fn_mod);

  // r_n: subtrac K_ns*P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fn_mod_nm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
  knm_add->Multiply(false, *(splitvel[1]), *fn_mod_nm);
  fn_mod->Update(-1.0, *fn_mod_nm, 1.0);

  // export additions to r_n subvector to r_new
  Teuchos::RCP<Epetra_Vector> fn_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*fn_mod, *fn_modexp);
  resnew->Update(1.0, *fn_modexp, 1.0);

  // export r_n subvector to r_new
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*(splitres[0]), *fnexp);
  resnew->Update(1.0, *fnexp, 1.0);

  if (dconmaster_ and firstnonliniter_)
  {
    Teuchos::RCP<Epetra_Vector> fn_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
    Core::LinAlg::Export(*dcnm, *fn_exp);
    resnew->Update(-1.0, *fn_exp, 1.0);
  }

  residual->Update(1.0, *resnew, 0.0);
}

/*-------------------------------------------------------*/
/*  Condensation operation block matrix     ehrl (04/11) */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::condensation_operation_block_matrix(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    const Teuchos::RCP<Epetra_Vector>& residual,
    const std::vector<Teuchos::RCP<Epetra_Vector>>& splitres,
    const std::vector<Teuchos::RCP<Epetra_Vector>>& splitvel)
{
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

  if (dconmaster_ and firstnonliniter_)
  {
    dcnm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
    dcmm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));

    split_vector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = adaptermeshtying_->GetMortarMatrixP();

  /*--------------------------------------------------------------------*/
  // block nm
  /*--------------------------------------------------------------------*/
  // compute modification for block nm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_mod =
      MLMultiply(sysmatnew->Matrix(0, 2), false, *P, false, false, false, true);

  // Add transformation matrix to nm
  sysmatnew->Matrix(0, 1).UnComplete();
  sysmatnew->Matrix(0, 1).Add(*knm_mod, false, 1.0, 1.0);

  if (dconmaster_ and firstnonliniter_) knm_mod->Multiply(false, *(splitdcmaster[1]), *dcnm);

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

  if (dconmaster_ and firstnonliniter_) kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  // complete matrix
  sysmatnew->Complete();

  //*************************************************
  //  condensation operation for the residual
  //*************************************************
  // r_m: add P^T*r_s
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  P->Multiply(true, *(splitres[2]), *fm_mod);

  // r_m: add P^T*K_ss*vp_i^s
  Teuchos::RCP<Epetra_Vector> fm_mod_ss = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  kss_mod->Multiply(false, *(splitvel[2]), *fm_mod_ss);
  fm_mod->Update(1.0, *fm_mod_ss, 1.0);

  // r_m: subtract P^T*K_ss*P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fm_mod_mm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));
  kmm_mod->Multiply(false, *(splitvel[1]), *fm_mod_mm);
  fm_mod->Update(-1.0, *fm_mod_mm, 1.0);

  // r_m: insert Dirichlet boundary conditions
  if (dconmaster_ and firstnonliniter_) fm_mod->Update(-1.0, *dcmm, 1.0);

  // export and add r_m subvector to residual
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*fm_mod, *fm_modexp);
  residual->Update(1.0, *fm_modexp, 1.0);


  // r_n: add K_ns*vp_i^s
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(sysmatnew->Matrix(0, 2)));
  // Core::LinAlg::SparseMatrix& knm = sysmatnew->Matrix(0,2);
  Teuchos::RCP<Epetra_Vector> fn_mod = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
  knm->Multiply(false, *(splitvel[2]), *fn_mod);

  // r_n: subtract K_ns*P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fn_mod_nm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
  knm_mod->Multiply(false, *(splitvel[1]), *fn_mod_nm);
  fn_mod->Update(-1.0, *fn_mod_nm, 1.0);

  // export and add r_n subvector to residual
  Teuchos::RCP<Epetra_Vector> fn_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  Core::LinAlg::Export(*fn_mod, *fn_modexp);
  residual->Update(1.0, *fn_modexp, 1.0);

  if (dconmaster_ and firstnonliniter_)
  {
    Teuchos::RCP<Epetra_Vector> fn_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
    Core::LinAlg::Export(*dcnm, *fn_exp);
    residual->Update(-1.0, *fn_exp, 1.0);
  }

  // export r_s = zero to residual
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  Core::LinAlg::Export(*fs_mod, *residual);
}

/*-------------------------------------------------------*/
/*  Compute and update Slave DOF's          ehrl (04/11) */
/* (including ALE case   vg 01/14)                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::UpdateSlaveDOF(
    const Teuchos::RCP<Epetra_Vector>& inc, const Teuchos::RCP<Epetra_Vector>& velnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.4)   - Update slave DOF");

  // get dof row map
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // split incremental and velocity-pressure vector
  std::vector<Teuchos::RCP<Epetra_Vector>> splitinc(3);
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvel(3);
  split_vector(inc, splitinc);
  split_vector(velnp, splitvel);

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process

  // split vector containing Dirichlet boundary conditions, if any
  std::vector<Teuchos::RCP<Epetra_Vector>> splitdcmaster(3);
  if (dconmaster_ and firstnonliniter_) split_vector(valuesdc_, splitdcmaster);

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = adaptermeshtying_->GetMortarMatrixP();

  // define new incremental vector
  Teuchos::RCP<Epetra_Vector> incnew = Core::LinAlg::CreateVector(*dofrowmap, true);

  // delta_vp^s: add P*delta_vp^m
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  P->Multiply(false, *(splitinc[1]), *fs_mod);

  // delta_vp^s: subtract vp_i^s
  fs_mod->Update(-1.0, *(splitvel[2]), 1.0);

  // delta_vp^s: add P*vp_i^m
  Teuchos::RCP<Epetra_Vector> fs_mod_m = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  P->Multiply(false, *(splitvel[1]), *fs_mod_m);
  fs_mod->Update(1.0, *fs_mod_m, 1.0);

  // set Dirichlet boundary conditions, if any
  if (dconmaster_ and firstnonliniter_)
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
  if (dconmaster_ and firstnonliniter_) firstnonliniter_ = false;

  // define incremental vector to new incremental vector
  inc->Update(1.0, *incnew, 0.0);
}

/*-------------------------------------------------------*/
/*  Output maps and projection matrix      ehrl (04/11)  */
/*-------------------------------------------------------*/

void FLD::Meshtying::OutputSetUp()
{
  if (myrank_ == 0)
  {
    // Output:

    /*std::cout << std::endl << "dof_row_map:" << std::endl;
    std::cout << *(discret_->dof_row_map())<< std::endl << std::endl;
    std::cout << std::endl << "masterDofRowMap:" << std::endl;
    std::cout << *(adaptermeshtying_->MasterDofRowMap())<< std::endl << std::endl;
    std::cout << "slaveDofRowMap:" << std::endl;
    std::cout << *(adaptermeshtying_->SlaveDofRowMap())<< std::endl << std::endl;
   */
    std::cout << "Projection matrix:" << std::endl;
    std::cout << *(adaptermeshtying_->GetMortarMatrixP()) << std::endl << std::endl;
  }

  /* {
   const std::string fname = "c_after.txt";

   std::ofstream f;
   f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
   f << "\n" << "Begin" << "\n";
   f << *c;
   f << "End" << "\n";
   f.close();
   }*/
}

/*-------------------------------------------------------*/
/*  Output: split sparse matrix            ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::output_sparse_matrix_split(Teuchos::RCP<Core::LinAlg::SparseOperator> conmat)
{
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> splitmatrix;

  split_matrix(conmat, splitmatrix);

  std::cout << "Teil nn " << std::endl << splitmatrix->Matrix(0, 0) << std::endl;
  std::cout << "Teil nm: " << std::endl << splitmatrix->Matrix(0, 1) << std::endl;
  std::cout << "Teil ns: " << std::endl << splitmatrix->Matrix(0, 2) << std::endl;

  std::cout << "Teil mn: " << std::endl << splitmatrix->Matrix(1, 0) << std::endl;
  std::cout << "Teil mm: " << std::endl << splitmatrix->Matrix(1, 1) << std::endl;
  std::cout << "Teil ms: " << std::endl << splitmatrix->Matrix(1, 2) << std::endl;

  std::cout << "Teil sn: " << std::endl << splitmatrix->Matrix(2, 0) << std::endl;
  std::cout << "Teil sm: " << std::endl << splitmatrix->Matrix(2, 1) << std::endl;
  std::cout << "Teil ss: " << std::endl << splitmatrix->Matrix(2, 2) << std::endl;

  FOUR_C_THROW("Matrix output finished");
}

/*-------------------------------------------------------*/
/*  Output: block matrix                   ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::OutputBlockMatrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> blockmatrix, Teuchos::RCP<Epetra_Vector> residual)
{
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockmatrixnew =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(blockmatrix);

  Core::LinAlg::SparseMatrix sysmat0 = blockmatrixnew->Matrix(0, 0);
  Core::LinAlg::SparseMatrix sysmat1 = blockmatrixnew->Matrix(0, 1);
  // Core::LinAlg::SparseMatrix sysmat2 = blockmatrixnew->Matrix(0,2);

  Core::LinAlg::SparseMatrix sysmat3 = blockmatrixnew->Matrix(1, 0);
  Core::LinAlg::SparseMatrix sysmat4 = blockmatrixnew->Matrix(1, 1);
  // Core::LinAlg::SparseMatrix sysmat5 = blockmatrixnew->Matrix(1,2);
  /*
    Core::LinAlg::SparseMatrix sysmat6 = blockmatrixnew->Matrix(2,0);
    Core::LinAlg::SparseMatrix sysmat7 = blockmatrixnew->Matrix(2,1);
    Core::LinAlg::SparseMatrix sysmat8 = blockmatrixnew->Matrix(2,2);*/

  std::cout << "Block nn" << *(sysmat0.EpetraMatrix()) << std::endl;
  std::cout << "Block nm" << *(sysmat1.EpetraMatrix()) << std::endl;
  // std::cout << "Block ns" << *(sysmat2.EpetraMatrix()) << std::endl;

  std::cout << "Block mn" << *(sysmat3.EpetraMatrix()) << std::endl;
  std::cout << "Block mm" << *(sysmat4.EpetraMatrix()) << std::endl;
}

/*-------------------------------------------------------*/
/*  Output: vector                         ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::OutputVectorSplit(Teuchos::RCP<Epetra_Vector> vector)
{
  std::vector<Teuchos::RCP<Epetra_Vector>> splitvector(3);
  split_vector(vector, splitvector);

  // std::cout << "vector " << std::endl << *vector << std::endl << std::endl;

  std::cout.precision(20);
  // std::cout << "Teil fn " << std::endl << *(splitvector[0]) << std::endl << std::endl;
  std::cout << "Teil fm: " << std::endl << *(splitvector[1]) << std::endl << std::endl;
  std::cout << "Teil fs: " << std::endl << *(splitvector[2]) << std::endl;
}

/*-------------------------------------------------------*/
/*  Output: Analyze matrix                 ehrl (11/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::AnalyzeMatrix(Teuchos::RCP<Core::LinAlg::SparseMatrix> sparsematrix)
{
  double localmatrixentries = 0.0;
  double parmatrixentries = 0.0;
  Teuchos::RCP<Epetra_CrsMatrix> matrix = sparsematrix->EpetraMatrix();
  {
    // number of row elements
    const int numdofrows = sparsematrix->RowMap().NumMyElements();

    for (int i = 0; i < numdofrows; ++i)
    {
      // max. number of non-zero values
      int maxnumentries = matrix->MaxNumEntries();
      int numOfNonZeros = 0;
      std::vector<int> indices(maxnumentries, 0);
      std::vector<double> values(maxnumentries, 0.0);

      int error =
          matrix->ExtractMyRowCopy(i, maxnumentries, numOfNonZeros, values.data(), indices.data());
      if (error != 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractMyRowCopy returned err=%d", error);

      for (int ii = 0; ii < numOfNonZeros; ii++)
      {
        localmatrixentries += values[ii];
      }
    }

    discret_->Comm().SumAll(&localmatrixentries, &parmatrixentries, 1);
  }
  double normfrob = matrix->NormFrobenius();
  double norminf = matrix->NormInf();
  double normone = matrix->NormOne();
  double matrixsize = matrix->NumGlobalRows() * matrix->NumGlobalCols();
  double nonzero = matrix->NumGlobalNonzeros();

  if (myrank_ == 0)
  {
    {
      std::cout.precision(20);
      std::cout << std::endl;
      std::cout << "-------------- Analyze Matrix ----------------------" << std::endl;
      std::cout << "| global matrix size:          " << matrixsize << std::endl;
      std::cout << "| number of global non-zeros:  " << nonzero << std::endl;
      std::cout << "| Matrix norm (Frobenius):     " << normfrob << std::endl;
      std::cout << "| Matrix norm (Inf):           " << norminf << std::endl;
      std::cout << "| Matrix norm (One):           " << normone << std::endl;
      std::cout << "| sum of all matrix entries:   " << parmatrixentries << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
    }
  }
}  // end AnalyzeMatrix()

// -------------------------------------------------------------------
// check absolut velinc norm                           ehrl   11/2011
// -------------------------------------------------------------------
/*
void FLD::FluidImplicitTimeInt::PrintAbsoluteL2Norm(Teuchos::RCP<Epetra_Vector>&   vector)
{
  double incvelnorm_L2;
  double incprenorm_L2;

  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_->ExtractOtherVector(vector);
  onlyvel->Norm2(&incvelnorm_L2);

  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_->ExtractCondVector(vector);
  onlypre->Norm2(&incprenorm_L2);

  printf("+------------+-------------------+--------------+\n");
  printf("| %10.14E   | %10.14E   |",
         incvelnorm_L2,incprenorm_L2);
  printf(")\n");
  printf("+------------+-------------------+--------------+\n");

  return;
}
  */

/*-------------------------------------------------------*/
/*  Set the flag for multifield problems     wirtz 01/16 */
/*                                                       */
/*-------------------------------------------------------*/
void FLD::Meshtying::IsMultifield(
    Teuchos::RCP<std::set<int>> condelements,           ///< conditioned elements of fluid
    const Core::LinAlg::MultiMapExtractor& domainmaps,  ///< domain maps for split of fluid matrix
    const Core::LinAlg::MultiMapExtractor& rangemaps,   ///< range maps for split of fluid matrix
    Teuchos::RCP<std::set<int>> condelements_shape,     ///< conditioned elements
    const Core::LinAlg::MultiMapExtractor&
        domainmaps_shape,  ///< domain maps for split of shape deriv. matrix
    const Core::LinAlg::MultiMapExtractor&
        rangemaps_shape,  ///< domain maps for split of shape deriv. matrix
    bool splitmatrix,     ///< flag for split of matrices
    bool ismultifield     ///< flag for multifield problems
)
{
  multifield_condelements_ = condelements;
  multifield_domainmaps_ = domainmaps;
  multifield_rangemaps_ = rangemaps;
  multifield_condelements_shape_ = condelements_shape;
  multifield_domainmaps_shape_ = domainmaps_shape;
  multifield_rangemaps_shape_ = rangemaps_shape;
  multifield_splitmatrix_ = splitmatrix;
  is_multifield_ = ismultifield;
}

/*-------------------------------------------------------*/
/*  Use the split of the fluid mesh tying for the sysmat */
/*                                           wirtz 01/16 */
/*-------------------------------------------------------*/
void FLD::Meshtying::MshtSplit(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat,
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& shapederivatives)
{
  if (is_multifield_)
  {
    // generate map for blockmatrix
    std::vector<Teuchos::RCP<const Epetra_Map>> fluidmaps;
    fluidmaps.push_back(gndofrowmap_);
    fluidmaps.push_back(gmdofrowmap_);
    fluidmaps.push_back(gsdofrowmap_);

    Core::LinAlg::MultiMapExtractor extractor;

    extractor.Setup(*dofrowmap_, fluidmaps);

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

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>> mat;
    mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(
        extractor, extractor, 108, false, true));
    // nodes on the interface
    Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->conditioned_element_map(*discret_);
    mat->SetCondElements(condelements);

    sysmat = mat;

    if (shapederivatives != Teuchos::null) MshtSplitShape(shapederivatives);
  }
}

/*-------------------------------------------------------*/
/*  Use the split of the fluid mesh tying for the shape  */
/*  derivatives                              wirtz 01/16 */
/*-------------------------------------------------------*/
void FLD::Meshtying::MshtSplitShape(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& shapederivatives)
{
  // generate map for blockmatrix
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidmaps;
  fluidmaps.push_back(gndofrowmap_);
  fluidmaps.push_back(gmdofrowmap_);
  fluidmaps.push_back(gsdofrowmap_);

  Core::LinAlg::MultiMapExtractor extractor;

  extractor.Setup(*dofrowmap_, fluidmaps);

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

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>> mat;
  mat = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(
      extractor, extractor, 108, false, true));
  // nodes on the interface
  Teuchos::RCP<std::set<int>> condelements = surfacesplitter_->conditioned_element_map(*discret_);
  mat->SetCondElements(condelements);

  shapederivatives = mat;
}

/*-------------------------------------------------------*/
/*  Use the split of the multifield problem for the      */
/*  sysmat                                   wirtz 01/16 */
/*-------------------------------------------------------*/
void FLD::Meshtying::MultifieldSplit(Teuchos::RCP<Core::LinAlg::SparseOperator>& sysmat)
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

    if (multifield_splitmatrix_)
    {
      Teuchos::RCP<Core::LinAlg::MapExtractor> extractor =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(
              *multifield_domainmaps_.FullMap(), multifield_domainmaps_.Map(1)));
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>> mat =
          mergedmatrix->Split<FLD::UTILS::InterfaceSplitStrategy>(*extractor, *extractor);
      mat->SetCondElements(multifield_condelements_);
      mat->Complete();

      sysmat = mat;
    }
    else
    {
      sysmat = mergedmatrix;
    }
  }
}

/*-------------------------------------------------------*/
/*  Use the split of the multifield problem for the      */
/*  shape derivatives                        wirtz 01/16 */
/*-------------------------------------------------------*/
void FLD::Meshtying::multifield_split_shape(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& shapederivatives)
{
  if (is_multifield_)
  {
    Teuchos::RCP<Epetra_Vector> ones =
        Teuchos::rcp(new Epetra_Vector(shapederivatives->Matrix(2, 2).RowMap()));
    ones->PutScalar(1.0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> onesdiag =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*ones));
    onesdiag->Complete();

    shapederivatives->Matrix(0, 2).UnComplete();
    shapederivatives->Matrix(0, 2).Zero();

    shapederivatives->Matrix(1, 2).UnComplete();
    shapederivatives->Matrix(1, 2).Zero();

    shapederivatives->Matrix(2, 2).UnComplete();
    shapederivatives->Matrix(2, 2).Zero();
    shapederivatives->Matrix(2, 2).Add(*onesdiag, false, 1.0, 1.0);

    shapederivatives->Matrix(2, 0).UnComplete();
    shapederivatives->Matrix(2, 0).Zero();

    shapederivatives->Matrix(2, 1).UnComplete();
    shapederivatives->Matrix(2, 1).Zero();

    shapederivatives->Complete();

    Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedshapederivatives = shapederivatives->Merge();

    Teuchos::RCP<Core::LinAlg::MapExtractor> extractor =
        Teuchos::rcp(new Core::LinAlg::MapExtractor(
            *multifield_domainmaps_shape_.FullMap(), multifield_domainmaps_shape_.Map(1)));
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>>
        matshapederivatives = mergedshapederivatives->Split<FLD::UTILS::InterfaceSplitStrategy>(
            *extractor, *extractor);
    matshapederivatives->SetCondElements(multifield_condelements_shape_);
    matshapederivatives->Complete();
    shapederivatives = matshapederivatives;
  }
}

/*-------------------------------------------------------*/
/*  Prepare condensation of the shape derivatives        */
/*                                           wirtz 01/16 */
/*-------------------------------------------------------*/
void FLD::Meshtying::condensation_operation_block_matrix_shape(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& shapederivatives)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Condensation Operation");

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

  if (dconmaster_ and firstnonliniter_)
  {
    dcnm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_, true));
    dcmm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_, true));

    split_vector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> P = adaptermeshtying_->GetMortarMatrixP();

  /*--------------------------------------------------------------------*/
  // block nm
  /*--------------------------------------------------------------------*/
  // compute modification for block nm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> knm_mod =
      MLMultiply(shapederivatives->Matrix(0, 2), false, *P, false, false, false, true);

  // Add transformation matrix to nm
  shapederivatives->Matrix(0, 1).UnComplete();
  shapederivatives->Matrix(0, 1).Add(*knm_mod, false, 1.0, 1.0);

  if (dconmaster_ and firstnonliniter_) knm_mod->Multiply(false, *(splitdcmaster[1]), *dcnm);  //???

  /*--------------------------------------------------------------------*/
  // block mn
  /*--------------------------------------------------------------------*/
  // compute modification for block kmn
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmn_mod =
      MLMultiply(*P, true, shapederivatives->Matrix(2, 0), false, false, false, true);

  // Add transformation matrix to mn
  shapederivatives->Matrix(1, 0).UnComplete();
  shapederivatives->Matrix(1, 0).Add(*kmn_mod, false, 1.0, 1.0);

  /*--------------------------------------------------------------------*/
  // block mm
  /*--------------------------------------------------------------------*/
  // compute modification for block kmm
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kss_mod =
      MLMultiply(*P, true, shapederivatives->Matrix(2, 2), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm_mod =
      MLMultiply(*kss_mod, false, *P, false, false, false, true);

  // Add transformation matrix to mm
  shapederivatives->Matrix(1, 1).UnComplete();
  shapederivatives->Matrix(1, 1).Add(*kmm_mod, false, 1.0, 1.0);

  if (dconmaster_ and firstnonliniter_) kmm_mod->Multiply(false, *(splitdcmaster[1]), *dcmm);

  // complete matrix
  shapederivatives->Complete();
}

FOUR_C_NAMESPACE_CLOSE
