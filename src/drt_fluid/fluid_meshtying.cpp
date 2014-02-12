/*!----------------------------------------------------------------------
\file fluid_meshtying.cpp

\brief Methods needed to apply rotationally symmetric periodic boundary
       conditions for fluid problems

<pre>
Maintainer: Andreas Ehrl
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>



*----------------------------------------------------------------------*/

#define DIRECTMANIPULATION
#define ZEROSYSMAT

#include "fluid_meshtying.H"
#include "fluid_utils.H"
#include "fluid_utils_mapextractor.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_fluid_ele/fluid_ele_parameter.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_krylov_projector.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>


FLD::Meshtying::Meshtying(Teuchos::RCP<DRT::Discretization>      dis,
                          LINALG::Solver&               solver,
                          int                           msht,
                          int                           nsd,
                          const UTILS::MapExtractor*    surfacesplitter
                          ):
  discret_(dis),
  solver_(solver),
  msht_(msht),
  surfacesplitter_(surfacesplitter),
  dofrowmap_(discret_->DofRowMap()),
  problemrowmap_(Teuchos::null),
  gndofrowmap_(Teuchos::null),
  gsmdofrowmap_(Teuchos::null),
  gsdofrowmap_(Teuchos::null),
  gmdofrowmap_(Teuchos::null),
  mergedmap_(Teuchos::null),
  valuesdc_(Teuchos::null),
  pcoupled_ (true),
  dconmaster_(false),
  firstnonliniter_(false),
  nsd_(nsd)
{
  // get the processor ID from the communicator
  myrank_  = discret_->Comm().MyPID();

  adaptermeshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
}

/*-------------------------------------------------------*/
/*  Setup mesh-tying problem                ehrl (04/11) */
/*-------------------------------------------------------*/

RCP<LINALG::SparseOperator> FLD::Meshtying::Setup(std::vector<int> coupleddof)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  1)   Setup Meshtying");
  if(coupleddof[nsd_]==0)
    pcoupled_=false;
  // Setup of meshtying adapter
 adaptermeshtying_->Setup(discret_,
                          discret_,
                          Teuchos::null,
                          coupleddof,
                          "Mortar",
                          discret_->Comm(),
                          false);

  //OutputSetUp();
  //AnalyzeMatrix(adaptermeshtying_->GetMortarTrafo());

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
  // since the implementation was only temporarily and not well tested
  // c) Saddle point system sparse matrix (sps_coupled)
  // d) Saddle point system block matrix (sps_pc)

  // number of nodes master < number of nodes slave
  // -> better results, Krylov does not work ??

  if(myrank_==0)
  {
    int numdofmaster = (adaptermeshtying_->MasterDofRowMap())->NumGlobalElements();
    int numdofslave = (adaptermeshtying_->SlaveDofRowMap())->NumGlobalElements();

    std::cout << std::endl << "number of master dof's:   " << numdofmaster << std::endl;
    std::cout << "number of slave dof's:   " << numdofslave << std::endl << std::endl;

    if(numdofmaster > numdofslave)
      std::cout << "The master side is discretized by more elements than the slave side" << std::endl;
    else
      std::cout << "The slave side is discretized by more elements than the master side" << std::endl;
  }

  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  case INPAR::FLUID::condensed_bmat_merged:
  {
    if (pcoupled_ == false)
      dserror("The system cannot be solved in a block matrix!! \n"
          "The null space does not have the right length. Fix it or use option Smat");

    // slave dof rowmap
    gsdofrowmap_ = adaptermeshtying_->SlaveDofRowMap();

    // master dof rowmap
    gmdofrowmap_ = adaptermeshtying_->MasterDofRowMap();

    // merge dofrowmap for slave and master discretization
    gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_,*gsdofrowmap_,false);

    // dofrowmap for discretisation without slave and master dofrowmap
    gndofrowmap_ = LINALG::SplitMap(*dofrowmap_, *gsmdofrowmap_);

    // map for 2x2 (uncoupled dof's & master dof's)
    mergedmap_ = LINALG::MergeMap(*gndofrowmap_,*gmdofrowmap_,false);

    //std::cout << "number of n dof   " << gndofrowmap_->NumGlobalElements() << std::endl;
    //std::cout << "number of m dof   " << gmdofrowmap_->NumGlobalElements() << std::endl;
    //std::cout << "number of s dof   " << gsdofrowmap_->NumGlobalElements() << std::endl;

    // generate map for blockmatrix
    std::vector<Teuchos::RCP<const Epetra_Map> > fluidmaps;
    fluidmaps.push_back(gndofrowmap_);
    fluidmaps.push_back(gmdofrowmap_);
    fluidmaps.push_back(gsdofrowmap_);

    LINALG::MultiMapExtractor extractor;

    extractor.Setup(*dofrowmap_,fluidmaps);

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

    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat;
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(extractor,extractor,108,false,true));
    // nodes on the interface
    Teuchos::RCP<std::set<int> > condelements = surfacesplitter_->ConditionedElementMap(*discret_);
    mat->SetCondElements(condelements);

    // Important: right way to do it (Tobias W.)
    // allocate 2x2 solution matrix with the default block matrix strategy in order to solve the reduced system
    // memory is not allocated(1), since the matrix gets a Teuchos::RCP on the respective blocks
    // of the 3x3 block matrix
    // ---------------
    // | knn  | knm' |
    // | kmn' | kmm' |
    // ---------------

    LINALG::MapExtractor rowmapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
    LINALG::MapExtractor dommapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > matsolve
      = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> (dommapext,rowmapext,1,false,true));
    sysmatsolve_ = matsolve;

    // fixing length of nullspace for block matrix (solver/preconditioner ML)
    if(msht_ ==INPAR::FLUID::condensed_bmat_merged)
    {
      string inv="BMatMerged";
      const Epetra_Map& oldmap = *(dofrowmap_);
      const Epetra_Map& newmap = *(mergedmap_);
      solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params());
      std::cout << std::endl;
    }
    else if (msht_ ==INPAR::FLUID::condensed_bmat)
    {
      // fixing length of Inverse1 nullspace (solver/preconditioner ML)
      {
        std::string inv="Inverse1";
        const Epetra_Map& oldmap = *(dofrowmap_);
        const Epetra_Map& newmap = matsolve->Matrix(0,0).EpetraMatrix()->RowMap();
        solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params().sublist("Inverse1"));
        std::cout << std::endl;
      }
      // fixing length of Inverse2 nullspace (solver/preconditioner ML)
      {
        std::string inv="Inverse2";
        const Epetra_Map& oldmap = *(dofrowmap_);
        const Epetra_Map& newmap = matsolve->Matrix(1,1).EpetraMatrix()->RowMap();
        solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params().sublist("Inverse2"));
        std::cout << std::endl;
      }
    }

    return mat;
  }
  break;
  case INPAR::FLUID::condensed_smat:
  {
    // slave dof rowmap
    gsdofrowmap_ = adaptermeshtying_->SlaveDofRowMap();

    // master dof rowmap
    gmdofrowmap_ = adaptermeshtying_->MasterDofRowMap();

    // merge dofrowmap for slave and master discretization
    gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_,*gsdofrowmap_,false);

    // dofrowmap for discretisation without slave and master dofrowmap
    gndofrowmap_ = LINALG::SplitMap(*dofrowmap_, *gsmdofrowmap_);

    if(myrank_==0)
    {
#ifdef ZEROSYSMAT
      std::cout << "Condensation operation takes place in the original sysmat -> graph is saved" << std::endl;
      std::cout << "The sysmat is set to zero and all parts are added -> exact" << std::endl << std::endl;
#else

#ifdef DIRECTMANIPULATION
      std::cout << "Condensation operation takes place in the original sysmat -> graph is saved" << std::endl << std::endl;
#else

      std::cout << "Condensation operation is carried out in a new allocated sparse matrix -> graph is not saved" << std::endl << std::endl;
#endif
#endif
    }

    return Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,108,false,true));
  }
  break;
  default:
     dserror("Choose a correct mesh-tying option");
     break;
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
void FLD::Meshtying::CheckOverlappingBC(
    Teuchos::RCP<Epetra_Map> map)
{
  bool overlap = false;

  // loop over all slave row nodes of the interface
  for (int j=0;j<adaptermeshtying_->Interface()->SlaveRowNodes()->NumMyElements();++j)
  {
    int gid = adaptermeshtying_->Interface()->SlaveRowNodes()->GID(j);
    DRT::Node* node = adaptermeshtying_->Interface()->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // check if this node's dofs are in given map
    for (int k=0;k<mtnode->NumDof();++k)
    {
      int currdof = mtnode->Dofs()[k];
      int lid = map->LID(currdof);

      // found slave node intersecting with given map
      if (lid>=0)
      {
        overlap = true;
        break;
      }
    }
  }

  // print warning message to screen
  if (overlap && myrank_==0)
  {
    if(myrank_==0)
    {
      dserror("Slave boundary and volume flow rate boundary conditions overlap!\n"
              "This leads to an over-constraint problem setup");
    }
  }

  return;
}

/*---------------------------------------------------*/
/*  Correct slave Dirichlet value       ehrl (12/12) */
/*---------------------------------------------------*/
void FLD::Meshtying::ProjectMasterToSlaveForOverlappingBC(
    Teuchos::RCP<Epetra_Vector>&  velnp,
    Teuchos::RCP<const Epetra_Map>  bmaps)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(bmaps);
  Teuchos::RCP<const Epetra_Map> gmdofrowmap = gmdofrowmap_;
  intersectionmaps.push_back(gmdofrowmap);
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    if(myrank_==0)
    {
      std::cout << "number of intersecting elements:  " << intersectionmap->NumGlobalElements() << std::endl;
      std::cout << "Overlapping boundary conditions are projected from the master to the slave side"
                << std::endl << std::endl;
    }

    // All dofs of the master are projected to the slave nodes. Therefore, all values of slave nodes
    // are overwritten although it is not necessary for the nodes which do not intersect with the
    // overlapping BC. In the case of Dirichlet or Dirichlet-like BC (Wormersly) new values are
    // assigned to the master side of internal interface but not to the slave side.

    //std::cout << "BEFORE projection from master to slave side" << std::endl;
    //OutputVectorSplit(velnp);

    UpdateSlaveDOF(velnp);

    //std::cout << "AFTER projection from master to slave side" << std::endl;
    //OutputVectorSplit(velnp);

  }

  return;
}

/*------------------------------------------------------------------------------*/
/*  Check if Dirichlet BC are defined on the master                ehrl (08/13) */
/*------------------------------------------------------------------------------*/
void FLD::Meshtying::DirichletOnMaster(
    Teuchos::RCP<const Epetra_Map>  bmaps)
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
  //      if(msht_ != INPAR::FLUID::no_meshtying)
  //        meshtying_->ProjectMasterToSlaveForOverlappingBC(velnp_, dbcmaps_->CondMap());
  //
  // (c)  DC are included in the condensation process (-> actual strategy)

  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(bmaps);
  Teuchos::RCP<const Epetra_Map> gmdofrowmap = gmdofrowmap_;
  intersectionmaps.push_back(gmdofrowmap);
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    dconmaster_ = true;
    if(myrank_==0)
    {
      std::cout << "Dirichlet or Dirichlet-like boundary condition defined on master side of the internal interface!\n "
                << "These conditions has to be also included at the slave side of the internal interface"
                << std::endl << std::endl;
    }

    if(msht_==INPAR::FLUID::condensed_smat)
    {
      if(myrank_==0)
      {
        dserror("Dirichlet BC are only support for the options bmat and bmat_merged so far!! \n"
                "The implementation for the option smat is straight forward. Just do it!!");
      }
    }
  }

  return;
}

/*------------------------------------------------------------------*/
/*  Include Dirichlet BC in condensation operation     ehrl (08/13) */
/*------------------------------------------------------------------*/
void FLD::Meshtying::IncludeDirichletInCondensation(
    Teuchos::RCP<Epetra_Vector>&  velnp,
    Teuchos::RCP<Epetra_Vector>&  veln)
{
  if(dconmaster_==true)
  {
    valuesdc_ = LINALG::CreateVector(*dofrowmap_,true);
    valuesdc_->Update(1.0,*velnp,1.0);
    valuesdc_->Update(-1.0,*veln,1.0);

    firstnonliniter_=true;
  }

  return;
}


/*---------------------------------------------------*/
/*  Prepare Meshtying system            ehrl (04/11) */
/*---------------------------------------------------*/
void FLD::Meshtying::PrepareMeshtyingSystem(
    Teuchos::RCP<LINALG::SparseOperator>&    sysmat,
    Teuchos::RCP<Epetra_Vector>&             residual)
{
  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  case INPAR::FLUID::condensed_bmat_merged:
    CondensationBlockMatrix(sysmat, residual);
    break;
  case INPAR::FLUID::condensed_smat:
    CondensationSparseMatrix(sysmat, residual);
    break;
  default:
    dserror("");
    break;
  }

  return;
}

void FLD::Meshtying::ApplyPTToResidual(
       Teuchos::RCP<LINALG::SparseOperator>     sysmat,
       Teuchos::RCP<Epetra_Vector>              residual,
       Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  Teuchos::RCP<Epetra_Vector> res=Teuchos::null;

  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  case INPAR::FLUID::condensed_bmat_merged:
     res      = LINALG::CreateVector(*mergedmap_,true);
    CondensationBlockMatrix(sysmat, residual);
    SplitVectorBasedOn3x3(residual, res);
    if(projector!=Teuchos::null)
      projector->ApplyPT(*res);
    LINALG::Export(*res,*residual);
    break;
  case INPAR::FLUID::condensed_smat:
    CondensationSparseMatrix(sysmat, residual);
    break;
  default:
    dserror("");
    break;
  }

  return;
}


/*-------------------------------------------------------*/
/*  Krylov projection                      ehrl (04/11)  */
/*-------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::Meshtying::AdaptKrylovProjector(Teuchos::RCP<Epetra_Vector> vec)
{

  if(pcoupled_)
  {
    // Remove slave nodes from vec
    Teuchos::RCP<Epetra_Vector> fm_slave = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
    // add fm subvector to feffnew
    LINALG::Export(*fm_slave,*vec);

    switch(msht_)
    {
    case INPAR::FLUID::condensed_bmat:
    case INPAR::FLUID::condensed_bmat_merged:
    {
      Teuchos::RCP<Epetra_Vector> vec_mesht  = LINALG::CreateVector(*mergedmap_,true);
      SplitVectorBasedOn3x3(vec, vec_mesht);
      vec=vec_mesht;

    }
    break;
    case INPAR::FLUID::condensed_smat:
      break;
    default:
      dserror("Krylov projection not supported for this meshtying option.");
      break;
    }
  }
  return vec;
}

/*-------------------------------------------------------*/
/*  solve mesh-tying system                ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::SolveMeshtying(
    LINALG::Solver&                 solver,
    Teuchos::RCP<LINALG::SparseOperator>     sysmat,
    Teuchos::RCP<Epetra_Vector>&             incvel,
    Teuchos::RCP<Epetra_Vector>              residual,
    int                             itnum,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3)   Solve meshtying system");

  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  {
    Teuchos::RCP<Epetra_Vector> res      = LINALG::CreateVector(*mergedmap_,true);
    Teuchos::RCP<Epetra_Vector> inc      = LINALG::CreateVector(*mergedmap_,true);

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatsolve = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmatsolve_);

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");

      SplitVectorBasedOn3x3(residual, res);
      // assign blocks to the solution matrix
      sysmatsolve->Assign(0,0, View, sysmatnew->Matrix(0,0));
      sysmatsolve->Assign(0,1, View, sysmatnew->Matrix(0,1));
      sysmatsolve->Assign(1,0, View, sysmatnew->Matrix(1,0));
      sysmatsolve->Assign(1,1, View, sysmatnew->Matrix(1,1));
      sysmatsolve->Complete();
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");
      solver_.Solve(sysmatsolve->EpetraOperator(),inc,res,true,itnum==1, projector);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Update");

      // Export the computed increment to the global increment
      LINALG::Export(*inc,*incvel);
      LINALG::Export(*res,*residual);

      // compute and update slave dof's
      UpdateSlaveDOF(incvel);
    }
  }
  break;
  case INPAR::FLUID::condensed_bmat_merged:
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatsolve = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmatsolve_);

      RCP<Epetra_Vector> res      = Teuchos::null;
      RCP<Epetra_Vector> inc      = Teuchos::null;

      RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::null;

      res      = LINALG::CreateVector(*mergedmap_,true);
      inc      = LINALG::CreateVector(*mergedmap_,true);

      mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap_,108,false,true));

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");
        SplitVectorBasedOn3x3(residual, res);


        // assign blocks to the solution matrix
        sysmatsolve->Assign(0,0, View, sysmatnew->Matrix(0,0));
        sysmatsolve->Assign(0,1, View, sysmatnew->Matrix(0,1));
        sysmatsolve->Assign(1,0, View, sysmatnew->Matrix(1,0));
        sysmatsolve->Assign(1,1, View, sysmatnew->Matrix(1,1));
        sysmatsolve->Complete();

        mergedmatrix = sysmatsolve->Merge();

    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");


      solver_.Solve(mergedmatrix->EpetraOperator(),inc,res,true,itnum==1, projector);

      LINALG::Export(*inc,*incvel);
      LINALG::Export(*res,*residual);
      // compute and update slave dof's
      UpdateSlaveDOF(incvel);

    }
  }
  break;
  case INPAR::FLUID::condensed_smat:
    {
      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Solve");
        solver_.Solve(sysmat->EpetraOperator(),incvel,residual,true,itnum==1, projector);
      }

      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.3)   - Update");
        // compute and update slave dof's
        UpdateSlaveDOF(incvel);
      }
    }
    break;
  default:
    dserror("");
    break;
  }
  return;
}


/*-------------------------------------------------------*/
/*  Condensation Sparse Matrix              ehrl (04/11) */
/*-------------------------------------------------------*/
void FLD::Meshtying::CondensationSparseMatrix(
    Teuchos::RCP<LINALG::SparseOperator>& sysmat,
    Teuchos::RCP<Epetra_Vector>&          residual
    )
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation sparse matrix");

  /**********************************************************************/
  /* Split sysmat and residual                                          */
  /**********************************************************************/

  // container for split matrix and vector
  std::vector<Teuchos::RCP<LINALG::SparseMatrix> > splitmatrix(9);
  std::vector<Teuchos::RCP<Epetra_Vector> > splitvector(3);

  SplitMatrix(sysmat,splitmatrix);
  SplitVector(residual, splitvector);

  /**********************************************************************/
  /* Condensate sparse matrix                                           */
  /**********************************************************************/

  CondensationOperationSparseMatrix(sysmat, residual, splitmatrix, splitvector);

  return;
}


/*-------------------------------------------------------*/
/*  Condensation Block Matrix               ehrl (04/11) */
/*-------------------------------------------------------*/
void FLD::Meshtying::CondensationBlockMatrix(Teuchos::RCP<LINALG::SparseOperator>&  sysmat,
                                             Teuchos::RCP<Epetra_Vector>&           residual)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation block matrix");

  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector> > splitvector(3);
  SplitVector(residual, splitvector);

  /**********************************************************************/
  /* Condensate blockmatrix                                             */
  /**********************************************************************/

  CondensationOperationBlockMatrix(sysmat, residual, splitvector);

  return;
}

/*-------------------------------------------------------*/
/*  Split Sparse Matrix                    ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::SplitMatrix(
    Teuchos::RCP<LINALG::SparseOperator>                 matrix,
    std::vector<Teuchos::RCP<LINALG::SparseMatrix> >&  splitmatrix)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Split Matrix");

  // cast Teuchos::RCP<LINALG::SparseOperator> to a Teuchos::RCP<LINALG::SparseMatrix>
  Teuchos::RCP<LINALG::SparseMatrix> matrixnew = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(matrix);

  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

  // split into interface and domain dof's
  LINALG::SplitMatrix2x2(matrixnew,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);

  // tempmap and tempmtx1 are dummy matrixes
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  // tempmap and tempmtx1 are dummy matrixes
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  // splitmatrix[ii]
  // -------------------------------
  // | knn [0] | knm [1] | kns [2] |
  // | kmn [3] | kmm [4] | kms [5] |
  // | ksn [6] | ksm [7] | kss [8] |
  // -------------------------------

  splitmatrix[0]=knn;      // -------------------------------
  splitmatrix[1]=knm;      // | knn (0) | knm (1) | kns (2) |
  splitmatrix[2]=kns;      // | kmn (3) | kmm (4) | kms (5) |
  splitmatrix[3]=kmn;      // | ksn (6) | ksm (7) | kss (8) |
  splitmatrix[4]=kmm;      // -------------------------------
  splitmatrix[5]=kms;
  splitmatrix[6]=ksn;
  splitmatrix[7]=ksm;
  splitmatrix[8]=kss;

  return;
}

/*-------------------------------------------------------*/
/*  Split Vector                           ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::SplitVector(
    Teuchos::RCP<Epetra_Vector>           vector,
    std::vector<Teuchos::RCP<Epetra_Vector> >&  splitvector)
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
   LINALG::SplitVector(*dofrowmap_,*vector,gsmdofrowmap_,fsm,gndofrowmap_,fn);

   // we want to split fsm into 2 groups s,m
   fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
   fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

   // do the vector splitting sm -> s+m
   LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

   // splitvector[ii]
   // fn [0]
   // fm [1]
   // fs [2]

   splitvector[0]=fn;
   splitvector[1]=fm;
   splitvector[2]=fs;

  return;
}

void FLD::Meshtying::SplitVectorBasedOn3x3(Teuchos::RCP<Epetra_Vector>   orgvector,
                                           Teuchos::RCP<Epetra_Vector>   vectorbasedon2x2)
{
  // container for split residual vector
  std::vector<Teuchos::RCP<Epetra_Vector> > splitvector(3);

  SplitVector(orgvector, splitvector);
  // build up the reduced residual
  LINALG::Export(*(splitvector[0]),*vectorbasedon2x2);
  LINALG::Export(*(splitvector[1]),*vectorbasedon2x2);

  return;
}

/*-------------------------------------------------------*/
/*  Condensation operation sparse matrix    ehrl (04/11) */
/*-------------------------------------------------------*/
void FLD::Meshtying::CondensationOperationSparseMatrix(
    Teuchos::RCP<LINALG::SparseOperator>&                sysmat,
    Teuchos::RCP<Epetra_Vector>&                         residual,
    std::vector<Teuchos::RCP<LINALG::SparseMatrix> >&    splitmatrix,
    std::vector<Teuchos::RCP<Epetra_Vector> >&           splitvector)
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

  // splitmatrix[ii]
  // -------------------------------
  // | knn [0] | knm [1] | kns [2] |
  // | kmn [3] | kmm [4] | kms [5] |
  // | ksn [6] | ksm [7] | kss [8] |
  // -------------------------------

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

  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  /**********************************************************************/
  /* Condensation operation for the sysmat                              */
  /**********************************************************************/

  // the sysmat is manipulated directly with out changing the graph
  // (subtract blocks to get zeros in the slave blocks)

#ifdef ZEROSYSMAT
  sysmat->UnComplete();

  // Part nn
  sysmat->Add(*(splitmatrix[0]),false,1.0,0.0);

  // Part nm
  {
    // knm: add kns*mbar
    Teuchos::RCP<LINALG::SparseMatrix> knm_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gndofrowmap_,100));
    knm_mod->Add(*(splitmatrix[1]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> knm_add = MLMultiply(*(splitmatrix[2]),false,*P,false,false,false,true);
    knm_mod->Add(*knm_add,false,1.0,1.0);
    knm_mod->Complete((splitmatrix[1])->DomainMap(),(splitmatrix[1])->RowMap());

    sysmat->Add(*knm_mod,false,1.0,1.0);
  }

  //Part mn
  {
    // kmn: add T(mbar)*ksn
    Teuchos::RCP<LINALG::SparseMatrix> kmn_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmn_mod->Add(*(splitmatrix[3]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmn_add = MLMultiply(*P,true,*(splitmatrix[6]),false,false,false,true);
    kmn_mod->Add(*kmn_add,false,1.0,1.0);
    kmn_mod->Complete((splitmatrix[3])->DomainMap(),(splitmatrix[3])->RowMap());

    sysmat->Add(*kmn_mod,false,1.0,1.0);
  }

  // Part mm
  {
    // kms: add T(mbar)*kss
    Teuchos::RCP<LINALG::SparseMatrix> kms_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kms_mod->Add(*(splitmatrix[5]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kms_add = MLMultiply(*P,true,*(splitmatrix[8]),false,false,false,true);
    kms_mod->Add(*kms_add,false,1.0,1.0);
    kms_mod->Complete((splitmatrix[5])->DomainMap(),(splitmatrix[5])->RowMap());

   // kmm: add T(mbar)*ksm + kmsmod*mbar
    Teuchos::RCP<LINALG::SparseMatrix> kmm_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmm_mod->Add(*(splitmatrix[4]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add = MLMultiply(*P,true,*(splitmatrix[7]),false,false,false,true);
    kmm_mod->Add(*kmm_add,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add2 = MLMultiply(*kms_mod,false,*P,false,false,false,true);
    kmm_mod->Add(*kmm_add2,false,1.0,1.0);
    kmm_mod->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());

    sysmat->Add(*kmm_mod,false,1.0,1.0);
  }

  {
   Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
   Teuchos::RCP<LINALG::SparseMatrix> onesdiag;
   ones->PutScalar(1.0);
   onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
   onesdiag->Complete();

   sysmat->Add(*onesdiag,false,1.0,1.0);
  }

  sysmat->Complete();

#else
#ifdef DIRECTMANIPULATION
  sysmat->UnComplete();

  // Part nm
  {
    // knm: add kns*mbar
    Teuchos::RCP<LINALG::SparseMatrix> knm_add = MLMultiply(*(splitmatrix[2]),false,*P,false,false,false,true);
    knm_add->Complete((splitmatrix[1])->DomainMap(),(splitmatrix[1])->RowMap());
    sysmat->Add(*knm_add,false,1.0,1.0);
  }

  // Part mn
  {
    // kmn: add T(mbar)*ksn
    Teuchos::RCP<LINALG::SparseMatrix> kmn_add = MLMultiply(*P,true,*(splitmatrix[6]),false,false,false,true);
    kmn_add->Complete((splitmatrix[3])->DomainMap(),(splitmatrix[3])->RowMap());
    sysmat->Add(*kmn_add,false,1.0,1.0);
  }

  // Part mm
  {
    // kms: add T(mbar)*kss
    Teuchos::RCP<LINALG::SparseMatrix> kms_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kms_mod->Add(*(splitmatrix[5]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kms_add = MLMultiply(*P,true,*(splitmatrix[8]),false,false,false,true);
    kms_mod->Add(*kms_add,false,1.0,1.0);
    kms_mod->Complete((splitmatrix[5])->DomainMap(),(splitmatrix[5])->RowMap());

    // kmm: add T(mbar)*ksm + kmsmod*mbar
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add = MLMultiply(*P,true,*(splitmatrix[7]),false,false,false,true);
    kmm_add->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add2 = MLMultiply(*kms_mod,false,*P,false,false,false,true);
    kmm_add2->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());

    sysmat->Add(*kmm_add,false,1.0,1.0);
    sysmat->Add(*kmm_add2,false,1.0,1.0);
  }

  // Dangerous??: Get zero in block ... by subtracting
  sysmat->Add(*splitmatrix[2],false,-1.0, 1.0);
  sysmat->Add(*splitmatrix[5],false,-1.0, 1.0);
  sysmat->Add(*splitmatrix[6],false,-1.0, 1.0);
  sysmat->Add(*splitmatrix[7],false,-1.0, 1.0);
  sysmat->Add(*splitmatrix[8],false,-1.0, 1.0);

  {
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag;
    // build identity matrix for slave dofs
    ones->PutScalar(1.0);
    //Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    sysmat->Add(*onesdiag,false,1.0,1.0);
  }

  // the sysmat is manipulated indirectly via a second sparse matrix

  sysmat->Complete();

  //OutputSparseMatrixSplit(sysmat);

#else
  // the sysmat is manipulated indirectly via a second sparse matrix
  // and therefore, the graph changes
  Teuchos::RCP<LINALG::SparseOperator> sysmatnew = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false));

  // Part nn
  sysmatnew->Add(*(splitmatrix[0]),false,1.0,1.0);

  // Part nm
  {
    // knm: add kns*mbar
    Teuchos::RCP<LINALG::SparseMatrix> knm_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gndofrowmap_,100));
    knm_mod->Add(*(splitmatrix[1]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> knm_add = MLMultiply(*(splitmatrix[2]),false,*P,false,false,false,true);
    knm_mod->Add(*knm_add,false,1.0,1.0);
    knm_mod->Complete((splitmatrix[1])->DomainMap(),(splitmatrix[1])->RowMap());

    sysmatnew->Add(*knm_mod,false,1.0,1.0);
  }

  //Part mn
  {
    // kmn: add T(mbar)*ksn
    Teuchos::RCP<LINALG::SparseMatrix> kmn_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmn_mod->Add(*(splitmatrix[3]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmn_add = MLMultiply(*P,true,*(splitmatrix[6]),false,false,false,true);
    kmn_mod->Add(*kmn_add,false,1.0,1.0);
    kmn_mod->Complete((splitmatrix[3])->DomainMap(),(splitmatrix[3])->RowMap());

    sysmatnew->Add(*kmn_mod,false,1.0,1.0);
  }

  // Part mm
  {
    // kms: add T(mbar)*kss
    Teuchos::RCP<LINALG::SparseMatrix> kms_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kms_mod->Add(*(splitmatrix[5]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kms_add = MLMultiply(*P,true,*(splitmatrix[8]),false,false,false,true);
    kms_mod->Add(*kms_add,false,1.0,1.0);
    kms_mod->Complete((splitmatrix[5])->DomainMap(),(splitmatrix[5])->RowMap());

   // kmm: add T(mbar)*ksm + kmsmod*mbar
    Teuchos::RCP<LINALG::SparseMatrix> kmm_mod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmm_mod->Add(*(splitmatrix[4]),false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add = MLMultiply(*P,true,*(splitmatrix[7]),false,false,false,true);
    kmm_mod->Add(*kmm_add,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmm_add2 = MLMultiply(*kms_mod,false,*P,false,false,false,true);
    kmm_mod->Add(*kmm_add2,false,1.0,1.0);
    kmm_mod->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());

    sysmatnew->Add(*kmm_mod,false,1.0,1.0);
  }

  {
   Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
   Teuchos::RCP<LINALG::SparseMatrix> onesdiag;
   ones->PutScalar(1.0);
   onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
   onesdiag->Complete();

   sysmatnew->Add(*onesdiag,false,1.0,1.0);
  }

  sysmatnew->Complete();

  sysmat = sysmatnew;
#endif
#endif

  //*************************************************
  //  condensation operation for the residual
  //*************************************************

  // splitvector[ii]
   // fn [0]
   // fm [1]
   // fs [2]

  Teuchos::RCP<Epetra_Vector> resnew = LINALG::CreateVector(*dofrowmap_,true);

  // fm: add T(mbar)*fs
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  // add fn subvector to residual
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*(splitvector[0]),*fnexp);
  resnew->Update(1.0,*fnexp,1.0);

  // add fn subvector to residual
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*(splitvector[1]),*fmexp);
  resnew->Update(1.0,*fmexp,1.0);

  // add fm subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  resnew->Update(1.0,*fm_modexp,1.0);

  residual=resnew;


  return;
}

/*-------------------------------------------------------*/
/*  Condensation operation block matrix     ehrl (04/11) */
/*-------------------------------------------------------*/
void FLD::Meshtying::CondensationOperationBlockMatrix(
    Teuchos::RCP<LINALG::SparseOperator>&             sysmat,
    Teuchos::RCP<Epetra_Vector>&                      residual,
    std::vector<Teuchos::RCP<Epetra_Vector> >&        splitvector)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Condensation Operation");

  // cast Teuchos::RCP<LINALG::SparseOperator> to a Teuchos::RCP<LINALG::BlockSparseMatrixBase>
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

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
  std::vector<Teuchos::RCP<Epetra_Vector> > splitdcmaster(3);

  if(dconmaster_==true and firstnonliniter_==true)
  {
    dcnm = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_,true));
    dcmm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_,true));

    SplitVector(valuesdc_, splitdcmaster);
  }

  // get transformation matrix
  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  // block nm
  {
    // compute modification for block nm
    Teuchos::RCP<LINALG::SparseMatrix> knm_mod = MLMultiply(sysmatnew->Matrix(0,2),false,*P,false,false,false,true);

    // Add transformation matrix to nm
    sysmatnew->Matrix(0,1).UnComplete();
    sysmatnew->Matrix(0,1).Add(*knm_mod,false,1.0,1.0);

    if(dconmaster_==true and firstnonliniter_==true)
      knm_mod->Multiply(false,*(splitdcmaster[1]),*dcnm);
  }

  // block mn
  {
    // compute modification for block kmn
    Teuchos::RCP<LINALG::SparseMatrix> kmn_mod = MLMultiply(*P,true,sysmatnew->Matrix(2,0),false,false,false,true);

    // Add transformation matrix to mn
    sysmatnew->Matrix(1,0).UnComplete();
    sysmatnew->Matrix(1,0).Add(*kmn_mod,false,1.0,1.0);
  }

  // block mm
  {
    // compute modification for block kmm
    Teuchos::RCP<LINALG::SparseMatrix> kss_mod = MLMultiply(*P,true,sysmatnew->Matrix(2,2),false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> kmm_mod = MLMultiply(*kss_mod,false,*P,false,false,false,true);

    // Add transformation matrix to mm
    sysmatnew->Matrix(1,1).UnComplete();
    sysmatnew->Matrix(1,1).Add(*kmm_mod,false,1.0,1.0);

    if(dconmaster_==true and firstnonliniter_==true)
      kmm_mod->Multiply(false,*(splitdcmaster[1]),*dcmm);
  }

  sysmatnew->Complete();

  //*************************************************
  //  condensation operation for the residual
  //*************************************************

  // fm: add T(mbar)*fs
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  if(dconmaster_==true and firstnonliniter_==true)
    fm_mod->Update(-1.0, *dcmm, 1.0);

  // add fm subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  residual->Update(1.0,*fm_modexp,1.0);

  if(dconmaster_==true and firstnonliniter_==true)
  {
    Teuchos::RCP<Epetra_Vector> fn_exp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
    LINALG::Export(*dcnm,*fn_exp);
    residual->Update(-1.0, *fn_exp, 1.0);
  }

  // fs = zero
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
  LINALG::Export(*fs_mod,*residual);

  return;
}

/*-------------------------------------------------------*/
/*  Compute and update Slave DOF's          ehrl (04/11) */
/*-------------------------------------------------------*/
void FLD::Meshtying::UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>&   inc)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.4)   - Update slave DOF");

  const Epetra_Map*  dofrowmap = discret_->DofRowMap();

  /**********************************************************************/
  /* Split inc into 3 subvectors                                       */
  /**********************************************************************/

  // container for split matrix and vector
  std::vector<Teuchos::RCP<Epetra_Vector> > splitvector(3);

  SplitVector(inc, splitvector);

  // Dirichlet or Dirichlet-like condition on the master side of the internal interface:
  // First time step:
  // coupling condition: u_s - u_m = delta u_m^D
  // instead of          u_s - u_m = 0
  //
  // this has to be considered in the condensation and in update process

  std::vector<Teuchos::RCP<Epetra_Vector> > splitdcmaster(3);

  if(dconmaster_==true and firstnonliniter_==true)
    SplitVector(valuesdc_, splitdcmaster);

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including meshtying)            */
  /**********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  Teuchos::RCP<Epetra_Vector> incnew = LINALG::CreateVector(*dofrowmap,true);

  // fs: add T(mbar)*fs
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
  P->Multiply(false,*(splitvector[1]),*fs_mod);

  //std::cout.precision(20);
  //std::cout << "fsmod:"<< *fs_mod << std::endl;

  if(dconmaster_==true and firstnonliniter_==true)
  {
    // fs: add T(mbar)*fs
    Teuchos::RCP<Epetra_Vector> fsdc_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
    P->Multiply(false,*(splitdcmaster[1]),*fsdc_mod);
    fs_mod->Update(1.0, *fsdc_mod, 1.0);

    //std::cout.precision(20);
    //std::cout << "dirichlet:"<< *fsdc_mod << std::endl;
  }

  // add fn subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitvector[0]),*fnexp);
  incnew->Update(1.0,*fnexp,1.0);

  // add fn subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitvector[1]),*fmexp);
  incnew->Update(1.0,*fmexp,1.0);

  // add fm subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fs_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*fs_mod,*fs_modexp);
  incnew->Update(1.0,*fs_modexp,1.0);

  // these terms are applied only in the first Newton iteration
  if(dconmaster_==true and firstnonliniter_==true)
    firstnonliniter_ = false;

  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/

  inc=incnew;

  return;
}

/*-------------------------------------------------------*/
/*  Output maps and projection matrix      ehrl (04/11)  */
/*-------------------------------------------------------*/

void FLD::Meshtying::OutputSetUp()
{
  if(myrank_==0)
  {
    // Output:

    /*std::cout << std::endl << "DofRowMap:" << std::endl;
    std::cout << *(discret_->DofRowMap())<< std::endl << std::endl;
    std::cout << std::endl << "masterDofRowMap:" << std::endl;
    std::cout << *(adaptermeshtying_->MasterDofRowMap())<< std::endl << std::endl;
    std::cout << "slaveDofRowMap:" << std::endl;
    std::cout << *(adaptermeshtying_->SlaveDofRowMap())<< std::endl << std::endl;
   */
    std::cout << "Projection matrix:" << std::endl;
    std::cout << *(adaptermeshtying_->GetMortarTrafo())<< std::endl << std::endl;
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
void FLD::Meshtying::OutputSparseMatrixSplit(
    Teuchos::RCP<LINALG::SparseOperator>                conmat)
{
  std::vector<Teuchos::RCP<LINALG::SparseMatrix> > matrixsplit(9);

  SplitMatrix(conmat,matrixsplit);

  std::cout << "Teil nn " << std::endl << *(matrixsplit[0]) << std::endl;
  std::cout << "Teil nm: " << std::endl << *(matrixsplit[1]) << std::endl;
  std::cout << "Teil ns: " << std::endl << *(matrixsplit[2]) << std::endl;

  std::cout << "Teil mn: " << std::endl << *(matrixsplit[3]) << std::endl;
  std::cout << "Teil mm: " << std::endl << *(matrixsplit[4]) << std::endl;
  std::cout << "Teil ms: " << std::endl << *(matrixsplit[5]) << std::endl;

  std::cout << "Teil sn: " << std::endl << *(matrixsplit[6]) << std::endl;
  std::cout << "Teil sm: " << std::endl << *(matrixsplit[7]) << std::endl;
  std::cout << "Teil ss: " << std::endl << *(matrixsplit[8]) << std::endl;

  dserror("Matrix output finished");

  return;
}

/*-------------------------------------------------------*/
/*  Output: block matrix                   ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::OutputBlockMatrix(
    Teuchos::RCP<LINALG::SparseOperator>       blockmatrix,
    Teuchos::RCP<Epetra_Vector>                residual)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmatrixnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(blockmatrix);

  LINALG::SparseMatrix sysmat0 = blockmatrixnew->Matrix(0,0);
  LINALG::SparseMatrix sysmat1 = blockmatrixnew->Matrix(0,1);
  //LINALG::SparseMatrix sysmat2 = blockmatrixnew->Matrix(0,2);

  LINALG::SparseMatrix sysmat3 = blockmatrixnew->Matrix(1,0);
  LINALG::SparseMatrix sysmat4 = blockmatrixnew->Matrix(1,1);
  //LINALG::SparseMatrix sysmat5 = blockmatrixnew->Matrix(1,2);
/*
  LINALG::SparseMatrix sysmat6 = blockmatrixnew->Matrix(2,0);
  LINALG::SparseMatrix sysmat7 = blockmatrixnew->Matrix(2,1);
  LINALG::SparseMatrix sysmat8 = blockmatrixnew->Matrix(2,2);*/

  std::cout << "Block nn" << *(sysmat0.EpetraMatrix()) << std::endl;
  std::cout << "Block nm" << *(sysmat1.EpetraMatrix()) << std::endl;
  //std::cout << "Block ns" << *(sysmat2.EpetraMatrix()) << std::endl;

  std::cout << "Block mn" << *(sysmat3.EpetraMatrix()) << std::endl;
  std::cout << "Block mm" << *(sysmat4.EpetraMatrix()) << std::endl;
  /*std::cout << "Block ms" << *(sysmat5.EpetraMatrix()) << std::endl;

  std::cout << "Block sn" << *(sysmat6.EpetraMatrix()) << std::endl;
  std::cout << "Block sm" << *(sysmat7.EpetraMatrix()) << std::endl;
  std::cout << "Block ss" << *(sysmat8.EpetraMatrix()) << std::endl;*/

  //LINALG::PrintMatrixInMatlabFormat("sysmat_BlockMatrix",*sysmat->EpetraMatrix(),true);

  /*    if (sysmat->RowMap().SameAs(residual_->Map()))
          std::cout << "juhu" << std::endl;
        else
          std::cout << "nein" << std::endl;  */

  return;
}

/*-------------------------------------------------------*/
/*  Output: vector                         ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::OutputVectorSplit(
    Teuchos::RCP<Epetra_Vector>                 vector)
{
  std::vector<Teuchos::RCP<Epetra_Vector> > splitvector(3);
  SplitVector(vector, splitvector);

  //std::cout << "vector " << std::endl << *vector << std::endl << std::endl;

  std::cout.precision(20);
  //std::cout << "Teil fn " << std::endl << *(splitvector[0]) << std::endl << std::endl;
  std::cout << "Teil fm: " << std::endl << *(splitvector[1]) << std::endl << std::endl;
  std::cout << "Teil fs: " << std::endl << *(splitvector[2]) << std::endl;
  return;
}

/*-------------------------------------------------------*/
/*  Output: Analyze matrix                 ehrl (11/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::AnalyzeMatrix(
    Teuchos::RCP<LINALG::SparseMatrix>        sparsematrix)
{
  double localmatrixentries = 0.0;
  double parmatrixentries = 0.0;
  Teuchos::RCP< Epetra_CrsMatrix > matrix = sparsematrix->EpetraMatrix ();
  {
    // number of row elements
    const int numdofrows = sparsematrix->RowMap().NumMyElements();

    for (int i=0; i<numdofrows; ++i)
    {
      // max. number of non-zero values
      int maxnumentries = matrix->MaxNumEntries();
      int numOfNonZeros = 0;
      std::vector<int> indices(maxnumentries, 0);
      std::vector<double> values(maxnumentries, 0.0);

      int error = matrix->ExtractMyRowCopy(i, maxnumentries, numOfNonZeros, &values[0], &indices[0]);
      if (error!=0) dserror("Epetra_CrsMatrix::ExtractMyRowCopy returned err=%d",error);

      for (int ii = 0; ii<numOfNonZeros; ii++)
      {
        localmatrixentries += values[ii];
      }
    }

    discret_->Comm().SumAll(&localmatrixentries,&parmatrixentries,1);
  }
  double normfrob = matrix->NormFrobenius();
  double norminf = matrix->NormInf();
  double normone = matrix->NormOne();
  double matrixsize = matrix->NumGlobalRows()*matrix->NumGlobalCols();
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

  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(vector);
  onlyvel->Norm2(&incvelnorm_L2);

  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(vector);
  onlypre->Norm2(&incprenorm_L2);

  printf("+------------+-------------------+--------------+\n");
  printf("| %10.14E   | %10.14E   |",
         incvelnorm_L2,incprenorm_L2);
  printf(")\n");
  printf("+------------+-------------------+--------------+\n");

  return;
}
  */

