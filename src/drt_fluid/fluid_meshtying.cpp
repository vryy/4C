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
#define BLOCKMATRIX_2x2

#include "fluid_meshtying.H"
#include "fluid_utils.H"
#include "fluid_utils_mapextractor.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_fluid_ele/fluid_ele_parameter.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>


FLD::Meshtying::Meshtying(Teuchos::RCP<DRT::Discretization>      dis,
                          LINALG::Solver&               solver,
                          Teuchos::ParameterList&       params,
                          const UTILS::MapExtractor*    surfacesplitter):
  discret_(dis),
  solver_(solver),
  msht_(params.get<int>("mshtoption")),
  surfacesplitter_(surfacesplitter),
  dofrowmap_(discret_->DofRowMap()),
  problemrowmap_(Teuchos::null),
  gndofrowmap_(Teuchos::null),
  gsmdofrowmap_(Teuchos::null),
  gsdofrowmap_(Teuchos::null),
  gmdofrowmap_(Teuchos::null),
  glmdofrowmap_(Teuchos::null),
  mergedmap_(Teuchos::null),
  lag_(Teuchos::null),
  lagold_(Teuchos::null),
  theta_(params.get<double>("theta")),
  pcoupled_ (true)
{
  // get the processor ID from the communicator
  myrank_  = discret_->Comm().MyPID();

  adaptermeshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());
}

/*-------------------------------------------------------*/
/*  Setup mesh-tying problem                ehrl (04/11) */
/*-------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> FLD::Meshtying::Setup()
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  1)   Setup Meshtying");

  // Setup of meshtying adapter
  pcoupled_ = adaptermeshtying_->Setup(*discret_,
                          discret_->Comm(),
                          msht_,
                          false);

  //OutputSetUp();

  // 4 different systems to solve
  // a) Condensation with a block matrix (condensed_bmat)
  // b) Condensation with a sparse matrix (condensed_smat)

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

    cout << endl << "number of master dof's:   " << numdofmaster << endl;
    cout << "number of slave dof's:   " << numdofslave << endl << endl;

    if(numdofmaster > numdofslave)
      // dserror("The master side is discretized by more elements than the slave side!! Do you really want to do it?");
      cout << "The master side is discretized by more elements than the slave side" << endl;
    else
      cout << "The slave side is discretized by more elements than the master side" << endl;
  }

  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  case INPAR::FLUID::condensed_bmat_merged:
  case INPAR::FLUID::coupling_iontransport_laplace:
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

    //cout << "number of n dof   " << gndofrowmap_->NumGlobalElements() << endl;
    //cout << "number of m dof   " << gmdofrowmap_->NumGlobalElements() << endl;
    //cout << "number of s dof   " << gsdofrowmap_->NumGlobalElements() << endl;

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
#ifdef BLOCKMATRIX_2x2
    LINALG::MapExtractor rowmapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
    LINALG::MapExtractor dommapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > matsolve
      = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> (dommapext,rowmapext,1,false,true));
    sysmatsolve_ = matsolve;
#else
    if(msht_==INPAR::FLUID::condensed_bmat)
    {
      LINALG::MapExtractor rowmapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
      LINALG::MapExtractor dommapext(*mergedmap_,gmdofrowmap_,gndofrowmap_);
      Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > matsolve
        = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> (dommapext,rowmapext,1,false,true));
      sysmatsolve_ = matsolve;
    }
    else
      sysmatsolve_ = mat;
#endif

    //Teuchos::RCP<std::vector<double> > test1 = solver.Params().sublist("Inverse1").sublist("ML Parameters").get<Teuchos::RCP<std::vector<double> > >("nullspace");
    //cout << "Length of null space before  " << test1->size() << endl;
    //cout << "address  " << test1 << endl;

    //Teuchos::RCP<std::vector<double> > test2 = solver.Params().sublist("Inverse2").sublist("ML Parameters").get<Teuchos::RCP<std::vector<double> > >("nullspace");
    //cout << "Length of null space before  " << test2->size() << endl;
    //cout << "address  " << test2 << endl;

    // fixing length of Inverse1 nullspace
    if (msht_ ==INPAR::FLUID::condensed_bmat)
    {
      {
        std::string inv="Inverse1";
        const Epetra_Map& oldmap = *(dofrowmap_);
        const Epetra_Map& newmap = matsolve->Matrix(0,0).EpetraMatrix()->RowMap();
        solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params().sublist("Inverse1"));
      }
      // fixing length of Inverse2 nullspace
      {
        std::string inv="Inverse2";
        const Epetra_Map& oldmap = *(dofrowmap_);
        const Epetra_Map& newmap = matsolve->Matrix(1,1).EpetraMatrix()->RowMap();
        solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params().sublist("Inverse2"));
      }
    }
#ifdef BLOCKMATRIX_2x2
    else if(msht_ ==INPAR::FLUID::condensed_bmat_merged)
    {
      std::string inv="BMatMerged";
      const Epetra_Map& oldmap = *(dofrowmap_);
      const Epetra_Map& newmap = *(mergedmap_);
      solver_.FixMLNullspace(&inv[0],oldmap, newmap, solver_.Params());
    }
#endif

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
#ifdef DIRECTMANIPULATION
      cout << "Condensation operation takes place in the original sysmat -> graph is saved" << endl;
      cout << "Warning: Dirichlet on the interface does not work in combination with smat" << endl << endl;
#else
      cout << "Condensation operation is carried out in a new allocated sparse matrix -> graph is not saved" << endl;
      cout << "Warning: Dirichlet on the interface does not work in combination with smat" << endl << endl;
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
  case INPAR::FLUID::coupling_iontransport_laplace:
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

/*-------------------------------------------------------*/
/*  Krylov projection                      ehrl (04/11)  */
/*-------------------------------------------------------*/
void FLD::Meshtying::AdaptKrylovProjector(Teuchos::RCP<Epetra_Vector> vec)
{
  // Remove slave nodes from vec
  Teuchos::RCP<Epetra_Vector> fm_slave = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
  // add fm subvector to feffnew
  LINALG::Export(*fm_slave,*vec);

  switch(msht_)
  {
  case INPAR::FLUID::condensed_bmat:
#ifdef BLOCKMATRIX_2x2
  case INPAR::FLUID::condensed_bmat_merged:
#endif
  case INPAR::FLUID::coupling_iontransport_laplace:
  {
    Teuchos::RCP<Epetra_Vector> vec_mesht  = LINALG::CreateVector(*mergedmap_,true);
    SplitVectorBasedOn3x3(vec, vec_mesht);
  }
  break;
  case INPAR::FLUID::condensed_smat:
    break;
  default:
    dserror("Krylov projection not supported for this meshtying option.");
    break;
  }

  return;
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

      // compute and update slave dof's
      UpdateSlaveDOF(incvel);
    }
  }
  break;
  case INPAR::FLUID::condensed_bmat_merged:
  case INPAR::FLUID::coupling_iontransport_laplace:
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmatsolve = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmatsolve_);
#ifdef BLOCKMATRIX_2x2
    Teuchos::RCP<Epetra_Vector> res      = LINALG::CreateVector(*mergedmap_,true);
    Teuchos::RCP<Epetra_Vector> inc      = LINALG::CreateVector(*mergedmap_,true);

    Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap_,108,false,true));
#else
    Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,108,false,true));
#endif

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Preparation");

#ifdef BLOCKMATRIX_2x2
      SplitVectorBasedOn3x3(residual, res);
#else
      sysmatsolve->Assign(0,2, View, sysmatnew->Matrix(0,2));
      sysmatsolve->Assign(1,2, View, sysmatnew->Matrix(1,2));
      sysmatsolve->Assign(2,0, View, sysmatnew->Matrix(2,0));
      sysmatsolve->Assign(2,1, View, sysmatnew->Matrix(2,1));
      sysmatsolve->Assign(2,2, View, sysmatnew->Matrix(2,2));
#endif

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
#ifdef BLOCKMATRIX_2x2
      solver_.Solve(mergedmatrix->EpetraOperator(),inc,res,true,itnum==1, projector);
#else
      solver_.Solve(mergedmatrix->EpetraOperator(),incvel,residual,true,itnum==1, projector);
#endif

#ifdef BLOCKMATRIX_2x2
      LINALG::Export(*inc,*incvel);
#endif

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

  // TODO: Apply dirichlet condition -> until now it only works for systems without
  //       dirichlet conditions on the interface

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

  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  /**********************************************************************/
  /* Condensation operation for the sysmat                              */
  /**********************************************************************/

  // the sysmat is manipulated directly with out changing the graph
  // (subtract blocks to get zeros in the slave blocks)
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

  // TODO: Apply dirichlet condition -> until now it only works for systems without
  //       dirichlet conditions on the interface

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
  // ---------------
  // | knn  | knm' |
  // | kmn' | kmm' |
  // ---------------

  // get transformation matrix
  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  // block nm
  {
    // compute modification for block nm
    Teuchos::RCP<LINALG::SparseMatrix> knm_mod = MLMultiply(sysmatnew->Matrix(0,2),false,*P,false,false,false,true);

    // Add transformation matrix to nm
    sysmatnew->Matrix(0,1).UnComplete();
    sysmatnew->Matrix(0,1).Add(*knm_mod,false,1.0,1.0);
  }

  // block mm
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
  }

#ifndef BLOCKMATRIX_2x2
  // block ss
  {
    // build identity matrix for slave dofs
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(sysmatnew->Matrix(2,2).RowMap()));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // Fill ss with unit matrix
    sysmatnew->Matrix(2,2).UnComplete();
    sysmatnew->Matrix(2,2).Zero();
    sysmatnew->Matrix(2,2).Add(*onesdiag,false,1.0,1.0);
  }

  // Fill ns with zero's
  sysmatnew->Matrix(0,2).UnComplete();
  sysmatnew->Matrix(0,2).Zero();

  // Fill ms with zero's
  sysmatnew->Matrix(1,2).UnComplete();
  sysmatnew->Matrix(1,2).Zero();

  // Fill sn with zero's
  sysmatnew->Matrix(2,0).UnComplete();
  sysmatnew->Matrix(2,0).Zero();

  // Fill sm with zero's
  sysmatnew->Matrix(2,1).UnComplete();
  sysmatnew->Matrix(2,1).Zero();
#endif

  sysmatnew->Complete();

  //*************************************************
  //  condensation operation for the residual
  //*************************************************

  // fm: add T(mbar)*fs
  Teuchos::RCP<Epetra_Vector> fm_mod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  // add fm subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fm_modexp = Teuchos::rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  residual->Update(1.0,*fm_modexp,1.0);

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

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including meshtying)            */
  /**********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> P = adaptermeshtying_->GetMortarTrafo();

  Teuchos::RCP<Epetra_Vector> incnew = LINALG::CreateVector(*dofrowmap,true);

  // fs: add T(mbar)*fs
  Teuchos::RCP<Epetra_Vector> fs_mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));
  P->Multiply(false,*(splitvector[1]),*fs_mod);

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
    cout << endl << "DofRowMap:" << endl;
    cout << *(discret_->DofRowMap())<< endl << endl;
    cout << endl << "masterDofRowMap:" << endl;
    cout << *(adaptermeshtying_->MasterDofRowMap())<< endl << endl;
    cout << "slaveDofRowMap:" << endl;
    cout << *(adaptermeshtying_->SlaveDofRowMap())<< endl << endl;
    cout << "lmDofRowMap:" << endl;
    cout << *(adaptermeshtying_->LmDofRowMap())<< endl << endl;
    cout << "Projection matrix:" << endl;
    cout << *(adaptermeshtying_->GetMortarTrafo())<< endl << endl;
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

  cout << "Teil nn " << endl << *(matrixsplit[0]) << endl;
  cout << "Teil nm: " << endl << *(matrixsplit[1]) << endl;
  cout << "Teil ns: " << endl << *(matrixsplit[2]) << endl;

  cout << "Teil mn: " << endl << *(matrixsplit[3]) << endl;
  cout << "Teil mm: " << endl << *(matrixsplit[4]) << endl;
  cout << "Teil ms: " << endl << *(matrixsplit[5]) << endl;

  cout << "Teil sn: " << endl << *(matrixsplit[6]) << endl;
  cout << "Teil sm: " << endl << *(matrixsplit[7]) << endl;
  cout << "Teil ss: " << endl << *(matrixsplit[8]) << endl;

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
  LINALG::SparseMatrix sysmat2 = blockmatrixnew->Matrix(0,2);

  LINALG::SparseMatrix sysmat3 = blockmatrixnew->Matrix(1,0);
  LINALG::SparseMatrix sysmat4 = blockmatrixnew->Matrix(1,1);
  //LINALG::SparseMatrix sysmat5 = blockmatrixnew->Matrix(1,2);
/*
  LINALG::SparseMatrix sysmat6 = blockmatrixnew->Matrix(2,0);
  LINALG::SparseMatrix sysmat7 = blockmatrixnew->Matrix(2,1);
  LINALG::SparseMatrix sysmat8 = blockmatrixnew->Matrix(2,2);*/

  cout << "Block nn" << *(sysmat0.EpetraMatrix()) << endl;
  cout << "Block nm" << *(sysmat1.EpetraMatrix()) << endl;
  //cout << "Block ns" << *(sysmat2.EpetraMatrix()) << endl;

  cout << "Block mn" << *(sysmat3.EpetraMatrix()) << endl;
  cout << "Block mm" << *(sysmat4.EpetraMatrix()) << endl;
  /*cout << "Block ms" << *(sysmat5.EpetraMatrix()) << endl;

  cout << "Block sn" << *(sysmat6.EpetraMatrix()) << endl;
  cout << "Block sm" << *(sysmat7.EpetraMatrix()) << endl;
  cout << "Block ss" << *(sysmat8.EpetraMatrix()) << endl;*/

  //LINALG::PrintMatrixInMatlabFormat("sysmat_BlockMatrix",*sysmat->EpetraMatrix(),true);

  /*    if (sysmat->RowMap().SameAs(residual_->Map()))
          cout << "juhu" << endl;
        else
          cout << "nein" << endl;  */

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

  cout << "vector " << endl << *vector << endl << endl;

  cout << "Teil fn " << endl << *(splitvector[0]) << endl << endl;
  cout << "Teil fm: " << endl << *(splitvector[1]) << endl << endl;
  cout << "Teil fs: " << endl << *(splitvector[2]) << endl;
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
      cout.precision(20);
      cout << endl;
      cout << "-------------- Analyze Matrix ----------------------" << endl;
      cout << "| global matrix size:          " << matrixsize << endl;
      cout << "| number of global non-zeros:  " << nonzero << endl;
      cout << "| Matrix norm (Frobenius):     " << normfrob << endl;
      cout << "| Matrix norm (Inf):           " << norminf << endl;
      cout << "| Matrix norm (One):           " << normone << endl;
      cout << "| sum of all matrix entries:   " << parmatrixentries << endl;
      cout << "----------------------------------------------------" << endl;
    }
  }
}  // end AnalyzeMatrix()

/*-------------------------------------------------------------*/
/*  Output: Replace matrix entries         ehrl (11/11)        */
/*  Replace computed identity matrix by a real identity matrix */
/*-------------------------------------------------------------*/
void FLD::Meshtying::ReplaceMatrixEntries(
    Teuchos::RCP<LINALG::SparseMatrix>        sparsematrix)
{
  Teuchos::RCP< Epetra_CrsMatrix > Pmat = sparsematrix->EpetraMatrix ();
  const int numdofrows = sparsematrix->RowMap().NumMyElements();

  for (int i=0; i<numdofrows; ++i)
  {
    // max. number of non-zero values
    int maxnumentries = Pmat->MaxNumEntries();
    int numOfNonZeros = 0;
    std::vector<int> indices(maxnumentries);
    std::vector<double> values(maxnumentries);

    int error = Pmat->ExtractMyRowCopy(i, maxnumentries, numOfNonZeros, &values[0], &indices[0]);
    if (error!=0) dserror("Epetra_CrsMatrix::ExtractMyRowCopy returned err=%d",error);

    // Saftey check: only one nonZero value in specific row
    if (numOfNonZeros!=1)
      dserror("Replace PMat: more than one nonZero value in specific row");

    double unity = 1.0;
    int err = Pmat->ReplaceMyValues(i, 1, &unity, &indices[0]);
    if (err!=0) dserror("Epetra_CrsMatrix::ExtractMyRowCopy returned err=%d",error);
  }

  double normfrob = sparsematrix->NormFrobenius();
  double norminf = sparsematrix->NormInf();
  double normone = sparsematrix->NormOne();

  if (myrank_ == 0)
  {
    {
      cout.precision(16);
      cout << "------------------------------------------------------" << endl;
      cout << "| projection matrix is replaced by a identity matrix:" << endl;
      cout << "| Matrix norm (Frobenius):  " << normfrob << endl;
      cout << "| Matrix norm (Inf):  " << norminf << endl;
      cout << "| Matrix norm (One):  " << normone << endl;
      cout << "------------------------------------------------------" << endl;
    }
  }
  //dserror("Remove dserror from OutputSetup()!!");
} // end ReplaceMatrixEntries()

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

