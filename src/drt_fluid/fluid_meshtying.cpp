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
#ifdef CCADISCRET

//#define DIRECTMANIPULATION

#include "fluid_meshtying.H"

#include <Teuchos_TimeMonitor.hpp>
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include <Epetra_Vector.h>


FLD::Meshtying::Meshtying(RCP<DRT::Discretization>      dis,
                          int       meshtyingoption):
  discret_(dis),
  msht_(meshtyingoption),
  problemrowmap_(Teuchos::null),
  gndofrowmap_(Teuchos::null),
  gsmdofrowmap_(Teuchos::null),
  gsdofrowmap_(Teuchos::null),
  gmdofrowmap_(Teuchos::null),
  glmdofrowmap_(Teuchos::null),
  lag_(Teuchos::null),
  lagold_(Teuchos::null),
  theta_(0.66)
{
  // get the processor ID from the communicator
  myrank_  = discret_->Comm().MyPID();
  dofrowmap_ = discret_->DofRowMap();
}


/*-------------------------------------------------------*/
/*  Setup                                                */
/*-------------------------------------------------------*/

RCP<LINALG::SparseOperator> FLD::Meshtying::Setup()
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  1)   Setup Meshtying");

  // Setup of meshtying adapter
  adaptermeshtying_.Setup(*discret_,
                          discret_->Comm(),
                          msht_,
                          false);

  //OutputSetUp();

  // 4 different systems to solve
  switch (msht_)
  {
  case INPAR::FLUID::condensed_bmat:
  {
    // slave dof rowmap
    gsdofrowmap_ = adaptermeshtying_.SlaveDofRowMap();

    // master dof rowmap
    gmdofrowmap_ = adaptermeshtying_.MasterDofRowMap();

    // merge dofrowmap for slave and master discretization
    gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_,*gsdofrowmap_,false);

    // dofrowmap for discretisation without slave and master dofrowmap
    RCP<Epetra_Map> gndofrowmap = LINALG::SplitMap(*dofrowmap_, *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap, *gmdofrowmap_);

    // generate map for blockmatrix
    std::vector<Teuchos::RCP<const Epetra_Map> > fluidmaps;
    fluidmaps.push_back(gndofrowmap_);
    fluidmaps.push_back(gmdofrowmap_);
    fluidmaps.push_back(gsdofrowmap_);

    //Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMapsKeepOrder(fluidmaps);

    //const std::vector<Teuchos::RCP<const Epetra_Map> > maps = fluidmaps;
    extractor_.Setup(*dofrowmap_,fluidmaps);

    // check, if extractor maps are valid
    extractor_.CheckForValidMapExtractor();

    // allocate block matrix
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat;
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
             (extractor_,extractor_,81,false,false));

    return mat;
  }
  break;
  case INPAR::FLUID::condensed_smat:
  {
    // slave dof rowmap
    gsdofrowmap_ = adaptermeshtying_.SlaveDofRowMap();

    // master dof rowmap
    gmdofrowmap_ = adaptermeshtying_.MasterDofRowMap();

    // merge dofrowmap for slave and master discretization
    gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_,*gsdofrowmap_,false);

    // dofrowmap for discretisation without slave and master dofrowmap
    RCP<Epetra_Map> gndofrowmap = LINALG::SplitMap(*dofrowmap_, *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap, *gmdofrowmap_);
    return Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,108,false,true));
  }
  break;
  case INPAR::FLUID::sps_coupled:
  case INPAR::FLUID::sps_pc:
  {
    lag_ = LINALG::CreateVector(*(adaptermeshtying_.LmDofRowMap()),true);
    lagold_ = LINALG::CreateVector(*(adaptermeshtying_.LmDofRowMap()),true);
    mergedmap_ = LINALG::MergeMap(*dofrowmap_,*(adaptermeshtying_.LmDofRowMap()),false);
    problemrowmap_ = rcp(new Epetra_Map(*(discret_->DofRowMap())));

    // slave dof rowmap
    gsdofrowmap_ = adaptermeshtying_.SlaveDofRowMap();

    // master dof rowmap
    gmdofrowmap_ = adaptermeshtying_.MasterDofRowMap();
    glmdofrowmap_ = adaptermeshtying_.LmDofRowMap();

    // merge dofrowmap for slave and master discretization
    gsmdofrowmap_ = LINALG::MergeMap(*gmdofrowmap_,*gsdofrowmap_,false);

    // dofrowmap for discretisation without slave and master dofrowmap
    RCP<Epetra_Map> gndofrowmap = LINALG::SplitMap(*dofrowmap_, *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap, *gmdofrowmap_);
    return Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_,108,false,true));
  }
  break;
  default:
     dserror("");
  }
  return Teuchos::null;
}

/*-------------------------------------------------------*/
/*  Condensation                                         */
/*-------------------------------------------------------*/

void FLD::Meshtying::Condensation(
    RCP<LINALG::SparseOperator>&    sysmat,
    RCP<Epetra_Vector>&             residual)
{
  if(msht_ == INPAR::FLUID::condensed_bmat)
    CondensationBlockMatrix(sysmat, residual);
  else if(msht_ == INPAR::FLUID::condensed_smat)
    CondensationSparseMatrix(sysmat, residual);
  else
    dserror("");

  return;
}

/*-------------------------------------------------------*/
/*  PrepareSaddlePointSystem                                                */
/*-------------------------------------------------------*/

void FLD::Meshtying::ResidualSaddlePointSystem(
    RCP<Epetra_Vector>      feff)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Scale residuum");

  // add meshtying force terms
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  adaptermeshtying_.GetDMatrix()->Multiply(true,*lag_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  feff->Update(-theta_,*fsexp,1.0);

  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  adaptermeshtying_.GetMMatrix()->Multiply(true,*lag_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  feff->Update(theta_,*fmexp,1.0);

  // add old contact forces (t_n)
  RCP<Epetra_Vector> fsold = rcp(new Epetra_Vector(*gsdofrowmap_));
  adaptermeshtying_.GetDMatrix()->Multiply(true,*lagold_,*fsold);
  RCP<Epetra_Vector> fsoldexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fsold,*fsoldexp);
  feff->Update(-(1.0-theta_),*fsoldexp,1.0);

  RCP<Epetra_Vector> fmold = rcp(new Epetra_Vector(*gmdofrowmap_));
  adaptermeshtying_.GetMMatrix()->Multiply(true,*lagold_,*fmold);
  RCP<Epetra_Vector> fmoldexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmold,*fmoldexp);
  feff->Update((1.0-theta_),*fmoldexp,1.0);

  return;
}

/*-------------------------------------------------------*/
/*  Krylov projection                                      */
/*-------------------------------------------------------*/

void FLD::Meshtying::KrylovProjection(RCP<Epetra_Vector>   c)
{
  // Remove slave nodes from c_
  RCP<Epetra_Vector> fm_slave = rcp(new Epetra_Vector(*gsdofrowmap_,true));
  // add fm subvector to feffnew
  LINALG::Export(*fm_slave,*c);

  return;
}

void FLD::Meshtying::SolveMeshtying(
    LINALG::Solver&                 solver,
    RCP<LINALG::SparseOperator>     sysmat,
    RCP<Epetra_Vector>&              incvel,
    RCP<Epetra_Vector>              residual,
    int                             itnum,
    RCP<Epetra_MultiVector>         w,
    RCP<Epetra_MultiVector>         c,
    bool                            project)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3)   Solve meshtying system");

  switch (msht_)
  {
  case INPAR::FLUID::sps_coupled:
  {
    RCP<LINALG::SparseMatrix>   mergedsysmat    = rcp(new LINALG::SparseMatrix(*mergedmap_,100,false,true));
    RCP<Epetra_Vector>          mergedresidual  = LINALG::CreateVector(*mergedmap_,true);
    RCP<Epetra_Vector>          mergedincvel    = LINALG::CreateVector(*mergedmap_,true);

    PrepareSaddlePointSystem(sysmat, mergedsysmat, residual, mergedresidual);

    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");
      solver.Solve(mergedsysmat->EpetraOperator(),mergedincvel,mergedresidual,true,itnum==1, w, c, project);
    }

    UpdateSaddlePointSystem(incvel,mergedincvel);
  }
  break;
  case INPAR::FLUID::sps_pc:
  {
    // row map (equals domain map) extractor
    LINALG::MapExtractor rowmapext(*mergedmap_,glmdofrowmap_,problemrowmap_);
    LINALG::MapExtractor dommapext(*mergedmap_,glmdofrowmap_,problemrowmap_);

    // build block matrix for SIMPLER
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blocksysmat =
          rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));

    RCP<Epetra_Vector>          mergedresidual  = LINALG::CreateVector(*mergedmap_,true);
    RCP<Epetra_Vector>          mergedincvel    = LINALG::CreateVector(*mergedmap_,true);

    PrepareSaddlePointSystemPC(sysmat, blocksysmat, residual, mergedresidual);

    // make solver SIMPLER-ready
    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");
      solver.PutSolverParamsToSubParams("SIMPLER", DRT::Problem::Instance()->FluidPressureSolverParams());
      solver.Params().sublist("SIMPLER").set<bool>("MESHTYING",true);
      solver.Solve(blocksysmat->EpetraOperator(),mergedincvel,mergedresidual,true,itnum==1);
    }

    UpdateSaddlePointSystem(incvel,mergedincvel);
  }
  break;
  case INPAR::FLUID::condensed_bmat:
  {
    {
      TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");

      // Use Simpler
      solver.PutSolverParamsToSubParams("SIMPLER", DRT::Problem::Instance()->FluidPressureSolverParams());
      solver.Params().sublist("SIMPLER").set<bool>("MESHTYING",true);
      solver.Solve(sysmat->EpetraOperator(),incvel,residual,true,itnum==1);
    }
    UpdateSlaveDOF(incvel);
  }
  break;
  case INPAR::FLUID::condensed_smat:
    {
      {
        TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.2)   - Solve");
        solver.Solve(sysmat->EpetraOperator(),incvel,residual,true,itnum==1, w, c, project);
      }

      UpdateSlaveDOF(incvel);
    }
    break;
  default:
    dserror("");
  }

  return;
}


void FLD::Meshtying::UpdateLag()
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  5)   Time update LM");
  lagold_->Update(1.0,*lag_,0.0);
  return;
}

/*-------------------------------------------------------*/
/*  PrepareSaddlePointSystem:   coupled                  */
/*-------------------------------------------------------*/

void FLD::Meshtying::PrepareSaddlePointSystem(
    RCP<LINALG::SparseOperator>    sysmat,
    RCP<LINALG::SparseMatrix>      mergedsysmat,
    RCP<Epetra_Vector>             residual,
    RCP<Epetra_Vector>             mergedresidual)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - prepare Saddle point system");

  Teuchos::RCP<LINALG::SparseMatrix> conmat
                    = (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat));

  //RCP<Epetra_Vector>        mergedzeros       = LINALG::CreateVector(*mergedmap,true);

  // build merged matrix
  mergedsysmat ->Add(*conmat,false,1.0,1.0);
  mergedsysmat ->Add(*(adaptermeshtying_.GetConMatrix()),false,theta_,1.0);
  mergedsysmat ->Add(*(adaptermeshtying_.GetConMatrix()),true,1.0,1.0);
  mergedsysmat ->Complete();

 /*
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(discret_->DofRowMap())));

  RCP<Epetra_Vector> laginc = rcp(new Epetra_Vector(*(meshtying_.LmDofRowMap())));
  LINALG::MapExtractor mapext(*mergedmap,problemrowmap,meshtying_.LmDofRowMap());
  mapext.ExtractCondVector(mergedincvel,incvel_);
  mapext.ExtractOtherVector(mergedincvel,laginc);
  laginc->ReplaceMap(*(meshtying_.SlaveDofRowMap()));

  //dserror("");
  lag_->Update(1.0,*laginc,0.0);*/

  // add meshtying force terms
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  adaptermeshtying_.GetDMatrix()->Multiply(true,*lag_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  residual->Update(theta_,*fsexp,1.0);

  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  adaptermeshtying_.GetMMatrix()->Multiply(true,*lag_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  residual->Update(-theta_,*fmexp,1.0);

  // build constraint rhs (=empty)
  RCP<Epetra_Vector> constrrhs = LINALG::CreateVector(*adaptermeshtying_.LmDofRowMap());

  // build merged rhs
  RCP<Epetra_Vector> residualexp = rcp(new Epetra_Vector(*mergedmap_));
  LINALG::Export(*residual,*residualexp);
  mergedresidual->Update(1.0,*residualexp,1.0);
  RCP<Epetra_Vector> constrexp = rcp(new Epetra_Vector(*mergedmap_));
  LINALG::Export(*constrrhs,*constrexp);
  mergedresidual->Update(1.0,*constrexp,1.0);

  return;
}


void FLD::Meshtying::PrepareSaddlePointSystemPC(
        RCP<LINALG::SparseOperator>         sysmat,
        RCP<LINALG::BlockSparseMatrixBase>  blocksysmat,
        RCP<Epetra_Vector>                  residual,
        RCP<Epetra_Vector>                  mergedresidual)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.1)   - Prepare matrix");

  RCP<LINALG::SparseMatrix> conmat = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat);

  RCP<LINALG::SparseMatrix> constrmt = rcp(new LINALG::SparseMatrix(*problemrowmap_,100,false,true));
  constrmt->Add(*(adaptermeshtying_.GetConMatrix()),false,theta_,0.0);
  constrmt->Complete(*glmdofrowmap_, *problemrowmap_);

  // build transposed constraint matrix
  RCP<LINALG::SparseMatrix> trconstrmt = rcp(new LINALG::SparseMatrix(*glmdofrowmap_,100,false,true));
  trconstrmt->Add(*(adaptermeshtying_.GetConMatrix()),true,1.0,0.0);
  trconstrmt->Complete(*problemrowmap_,*glmdofrowmap_);

  /* // apply Dirichlet conditions to (0,1) block
   RCP<Epetra_Vector> zeros   = rcp(new Epetra_Vector(*problemrowmap_,true));
   RCP<Epetra_Vector> rhscopy = rcp(new Epetra_Vector(*residual));
   LINALG::ApplyDirichlettoSystem(conmat,sold,rhscopy,zeros,dirichtoggle);
   constrmt->ApplyDirichlet(dirichtoggle,false);*/

  blocksysmat->Assign(0,0,View,*conmat);
  blocksysmat->Assign(0,1,View,*constrmt);
  blocksysmat->Assign(1,0,View,*trconstrmt);
  blocksysmat->Complete();

  // build constraint rhs (=empty)
  RCP<Epetra_Vector> constrrhs = rcp(new Epetra_Vector(*glmdofrowmap_));

  // add meshtying force terms
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  adaptermeshtying_.GetDMatrix()->Multiply(true,*lag_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  residual->Update(theta_,*fsexp,1.0);

  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  adaptermeshtying_.GetMMatrix()->Multiply(true,*lag_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  residual->Update(-theta_,*fmexp,1.0);

  // we also need merged rhs here
  RCP<Epetra_Vector> resexp = rcp(new Epetra_Vector(*mergedmap_));
  LINALG::Export(*residual,*resexp);
  mergedresidual->Update(1.0,*resexp,1.0);
  RCP<Epetra_Vector> constrexp = rcp(new Epetra_Vector(*mergedmap_));
  LINALG::Export(*constrrhs,*constrexp);
  mergedresidual->Update(1.0,*constrexp,1.0);

 /* // apply Dirichlet B.C. to mergedrhs and mergedsol
  RCP<Epetra_Vector> dirichtoggleexp = rcp(new Epetra_Vector(*mergedmap));
  LINALG::Export(*dirichtoggle,*dirichtoggleexp);*/

  return;
}

/*-------------------------------------------------------*/
/*  PrepareSaddlePointSystem                                                */
/*-------------------------------------------------------*/

void FLD::Meshtying::UpdateSaddlePointSystem(
    RCP<Epetra_Vector>      inc,
    RCP<Epetra_Vector>      mergedinc)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  4)   Newton iteration update LM");

  RCP<Epetra_Vector> laginc = rcp(new Epetra_Vector(*(adaptermeshtying_.LmDofRowMap())));
  LINALG::MapExtractor mapext(*mergedmap_,problemrowmap_,adaptermeshtying_.LmDofRowMap());
  mapext.ExtractCondVector(mergedinc,inc);
  mapext.ExtractOtherVector(mergedinc,laginc);
  laginc->ReplaceMap(*(adaptermeshtying_.SlaveDofRowMap()));

  //dserror("");
  lag_->Update(1.0,*laginc,0.0);

  return;
}

/*-------------------------------------------------------*/
/*  Condensation Sparse Matrix                           */
/*-------------------------------------------------------*/

void FLD::Meshtying::CondensationSparseMatrix(
    RCP<LINALG::SparseOperator>& sysmat,
    RCP<Epetra_Vector>&          residual
    )
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation sparse matrix");

  /**********************************************************************/
  /* Split sysmat and residual                                          */
  /**********************************************************************/

  // container for split matrix and vector
  std::vector<RCP<LINALG::SparseMatrix> > splitmatrix(9);
  std::vector<RCP<Epetra_Vector> > splitvector(3);

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

void FLD::Meshtying::CondensationBlockMatrix(RCP<LINALG::SparseOperator>&  sysmat,
                                             RCP<Epetra_Vector>&           residual)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2)   Condensation block matrix");

  /**********************************************************************/
  /* Split residual into 3 subvectors                                   */
  /**********************************************************************/

  // container for split residual vector
  std::vector<RCP<Epetra_Vector> > splitvector(3);
  SplitVector(residual, splitvector);

  /**********************************************************************/
  /* Condensate blockmatrix                                             */
  /**********************************************************************/

  CondensationOperationBlockMatrix(sysmat, residual, splitvector);

  return;
}

/*-------------------------------------------------------*/
/*  Split Sparse Matrix                                  */
/*-------------------------------------------------------*/

void FLD::Meshtying::SplitMatrix(
    RCP<LINALG::SparseOperator>                 matrix,
    std::vector<RCP<LINALG::SparseMatrix> >&  splitmatrix)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Split Matrix");

  // cast RCP<LINALG::SparseOperator> to a RCP<LINALG::SparseMatrix>
  RCP<LINALG::SparseMatrix> matrixnew = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(matrix);

  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

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
/*  PrepareSaddlePointSystem                                                */
/*-------------------------------------------------------*/

void FLD::Meshtying::SplitVector(
    RCP<Epetra_Vector>                 vector,
    std::vector<RCP<Epetra_Vector> >&   splitvector)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.2)   - Split Vector");

   // we want to split f into 3 groups s.m,n
   RCP<Epetra_Vector> fs, fm, fn;

   // temporarily we need the group sm
   RCP<Epetra_Vector> fsm;

   /**********************************************************************/
   /* Split feff into 3 subvectors                                       */
   /**********************************************************************/

   // do the vector splitting smn -> sm+n
   LINALG::SplitVector(*dofrowmap_,*vector,gsmdofrowmap_,fsm,gndofrowmap_,fn);

   // we want to split fsm into 2 groups s,m
   fs = rcp(new Epetra_Vector(*gsdofrowmap_));
   fm = rcp(new Epetra_Vector(*gmdofrowmap_));

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

/*-------------------------------------------------------*/
/*  PrepareSaddlePointSystem                                                */
/*-------------------------------------------------------*/

void FLD::Meshtying::CondensationOperationSparseMatrix(
    RCP<LINALG::SparseOperator>&                 sysmat,
    RCP<Epetra_Vector>&                        residual,
    std::vector<RCP<LINALG::SparseMatrix> >&  splitmatrix,
    std::vector<RCP<Epetra_Vector> >&         splitvector)
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

  // splitmatrix[ii]
  // -------------------------------
  // | knn [0] | knm [1] | kns [2] |
  // | kmn [3] | kmm [4] | kms [5] |
  // | ksn [6] | ksm [7] | kss [8] |
  // -------------------------------

  RCP<LINALG::SparseMatrix> P = adaptermeshtying_.GetMortarTrafo();

#ifdef  DIRECTMANIPULATION

  /**********************************************************************/
  /* Condensation operation for the sysmat                              */
  /**********************************************************************/

  sysmat->UnComplete();

  // Part nm
  {
    // knm: add kns*mbar
    RCP<LINALG::SparseMatrix> knm_add = MLMultiply(*(splitmatrix[2]),false,*P,false,false,false,true);
    knm_add->Complete((splitmatrix[1])->DomainMap(),(splitmatrix[1])->RowMap());
    sysmat->Add(*knm_add,false,1.0,1.0);
  }

  // Part mn
  {
    // kmn: add T(mbar)*ksn
    RCP<LINALG::SparseMatrix> kmn_add = MLMultiply(*P,true,*(splitmatrix[6]),false,false,false,true);
    kmn_add->Complete((splitmatrix[3])->DomainMap(),(splitmatrix[3])->RowMap());
    sysmat->Add(*kmn_add,false,1.0,1.0);
  }

  // Part mm
  {
    // kms: add T(mbar)*kss
    RCP<LINALG::SparseMatrix> kms_mod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kms_mod->Add(*(splitmatrix[5]),false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kms_add = MLMultiply(*P,true,*(splitmatrix[8]),false,false,false,true);
    kms_mod->Add(*kms_add,false,1.0,1.0);
    kms_mod->Complete((splitmatrix[5])->DomainMap(),(splitmatrix[5])->RowMap());

    // kmm: add T(mbar)*ksm + kmsmod*mbar
    RCP<LINALG::SparseMatrix> kmm_add = MLMultiply(*P,true,*(splitmatrix[7]),false,false,false,true);
    kmm_add->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());
    RCP<LINALG::SparseMatrix> kmm_add2 = MLMultiply(*kms_mod,false,*P,false,false,false,true);
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
    RCP<Epetra_Vector> ones = rcp (new Epetra_Vector(*gsdofrowmap_));
    RCP<LINALG::SparseMatrix> onesdiag;
    // build identity matrix for slave dofs
    ones->PutScalar(1.0);
    //RCP<LINALG::SparseMatrix> onesdiag = rcp(new LINALG::SparseMatrix(*ones));
    onesdiag = rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    sysmat->Add(*onesdiag,false,1.0,1.0);
  }

  sysmat->Complete();

  OutputSparseMatrixSplit(sysmat);

#else

  RCP<LINALG::SparseOperator> sysmatnew = rcp(new LINALG::SparseMatrix(*dofrowmap_,81,true,false));

  // Part nn
  sysmatnew->Add(*(splitmatrix[0]),false,1.0,1.0);

  // Part nm
  {
    // knm: add kns*mbar
    RCP<LINALG::SparseMatrix> knm_mod = rcp(new LINALG::SparseMatrix(*gndofrowmap_,100));
    knm_mod->Add(*(splitmatrix[1]),false,1.0,1.0);
    RCP<LINALG::SparseMatrix> knm_add = MLMultiply(*(splitmatrix[2]),false,*P,false,false,false,true);
    knm_mod->Add(*knm_add,false,1.0,1.0);
    knm_mod->Complete((splitmatrix[1])->DomainMap(),(splitmatrix[1])->RowMap());

    sysmatnew->Add(*knm_mod,false,1.0,1.0);
  }

  //Part mn
  {
    // kmn: add T(mbar)*ksn
    RCP<LINALG::SparseMatrix> kmn_mod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmn_mod->Add(*(splitmatrix[3]),false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmn_add = MLMultiply(*P,true,*(splitmatrix[6]),false,false,false,true);
    kmn_mod->Add(*kmn_add,false,1.0,1.0);
    kmn_mod->Complete((splitmatrix[3])->DomainMap(),(splitmatrix[3])->RowMap());

    sysmatnew->Add(*kmn_mod,false,1.0,1.0);
  }

  // Part mm
  {
    // kms: add T(mbar)*kss
    RCP<LINALG::SparseMatrix> kms_mod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kms_mod->Add(*(splitmatrix[5]),false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kms_add = MLMultiply(*P,true,*(splitmatrix[8]),false,false,false,true);
    kms_mod->Add(*kms_add,false,1.0,1.0);
    kms_mod->Complete((splitmatrix[5])->DomainMap(),(splitmatrix[5])->RowMap());

   // kmm: add T(mbar)*ksm + kmsmod*mbar
    RCP<LINALG::SparseMatrix> kmm_mod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmm_mod->Add(*(splitmatrix[4]),false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmm_add = MLMultiply(*P,true,*(splitmatrix[7]),false,false,false,true);
    kmm_mod->Add(*kmm_add,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmm_add2 = MLMultiply(*kms_mod,false,*P,false,false,false,true);
    kmm_mod->Add(*kmm_add2,false,1.0,1.0);
    kmm_mod->Complete((splitmatrix[4])->DomainMap(),(splitmatrix[4])->RowMap());

    sysmatnew->Add(*kmm_mod,false,1.0,1.0);
  }

  {
   RCP<Epetra_Vector> ones = rcp (new Epetra_Vector(*gsdofrowmap_));
   RCP<LINALG::SparseMatrix> onesdiag;
   ones->PutScalar(1.0);
   onesdiag = rcp(new LINALG::SparseMatrix(*ones));
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

#ifdef  DIRECTMANIPULATION

  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fm_mod = rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  // add fm subvector to feffnew
  RCP<Epetra_Vector> fm_modexp = rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  residual->Update(1.0,*fm_modexp,1.0);

  // fs = zero
  RCP<Epetra_Vector> fs_mod = rcp(new Epetra_Vector(*gsdofrowmap_,true));
  LINALG::Export(*fs_mod,*residual);

#else
  RCP<Epetra_Vector> resnew = LINALG::CreateVector(*dofrowmap_,true);

  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fm_mod = rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  // add fn subvector to residual
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*(splitvector[0]),*fnexp);
  resnew->Update(1.0,*fnexp,1.0);

  // add fn subvector to residual
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*(splitvector[1]),*fmexp);
  resnew->Update(1.0,*fmexp,1.0);

  // add fm subvector to feffnew
  RCP<Epetra_Vector> fm_modexp = rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  resnew->Update(1.0,*fm_modexp,1.0);

  residual=resnew;
#endif

  return;
}

void FLD::Meshtying::CondensationOperationBlockMatrix(
    RCP<LINALG::SparseOperator>&      sysmat,
    RCP<Epetra_Vector>&                      residual,
    std::vector<RCP<Epetra_Vector> >&        splitvector)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  2.1)   - Condensation Operation");

  // cast RCP<LINALG::SparseOperator> to a RCP<LINALG::BlockSparseMatrixBase>
  RCP<LINALG::BlockSparseMatrixBase> sysmatnew = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat);

  /**********************************************************************/
  /* Build the final sysmat and residual                                */
  /**********************************************************************/

  // ---------------------        ------------------
  // | nn | nm | ns | 0  |        | nn  | nm' | 0  |
  // | mn | mm | ms | D  |   =    | mn' | mm' | 0  |
  // | sn | sm | ss | -M |        |  0  |  0  | 1  |
  // |  0 | DT |-MT | 0  |        ------------------
  // ---------------------

  // get transformation matrix
  RCP<LINALG::SparseMatrix> P = adaptermeshtying_.GetMortarTrafo();

  // block nm
  {
    // compute modification for
    RCP<LINALG::SparseMatrix> knm_mod = MLMultiply(sysmatnew->Matrix(0,2),false,*P,false,false,false,true);

    // Add transformation matrix to nm
    sysmatnew->Matrix(0,1).UnComplete();
    sysmatnew->Matrix(0,1).Add(*knm_mod,false,1.0,1.0);
  }

  // block mm
  {
    // compute modification for block kmn
    RCP<LINALG::SparseMatrix> kmn_mod = MLMultiply(*P,true,sysmatnew->Matrix(2,0),false,false,false,true);

    // Add transformation matrix to mn
    sysmatnew->Matrix(1,0).UnComplete();
    sysmatnew->Matrix(1,0).Add(*kmn_mod,false,1.0,1.0);
  }

  // block mm
  {
    // compute modification for block kmm
    RCP<LINALG::SparseMatrix> kss_mod = MLMultiply(*P,true,sysmatnew->Matrix(2,2),false,false,false,true);
    RCP<LINALG::SparseMatrix> kmm_mod = MLMultiply(*kss_mod,false,*P,false,false,false,true);

    // Add transformation matrix to mm
    sysmatnew->Matrix(1,1).UnComplete();
    sysmatnew->Matrix(1,1).Add(*kmm_mod,false,1.0,1.0);
  }

  // block ss
  {
    // build identity matrix for slave dofs
    RCP<Epetra_Vector> ones = rcp (new Epetra_Vector(sysmatnew->Matrix(2,2).RowMap()));
    ones->PutScalar(1.0);
    RCP<LINALG::SparseMatrix> onesdiag = rcp(new LINALG::SparseMatrix(*ones));
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

  sysmatnew->Complete();

  //*************************************************
  //  condensation operation for the residual
  //*************************************************

  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fm_mod = rcp(new Epetra_Vector(*gmdofrowmap_,true));
  P->Multiply(true,*(splitvector[2]),*fm_mod);

  // add fm subvector to feffnew
  RCP<Epetra_Vector> fm_modexp = rcp(new Epetra_Vector(*dofrowmap_));
  LINALG::Export(*fm_mod,*fm_modexp);
  residual->Update(1.0,*fm_modexp,1.0);

  // fs = zero
  RCP<Epetra_Vector> fs_mod = rcp(new Epetra_Vector(*gsdofrowmap_,true));
  LINALG::Export(*fs_mod,*residual);

  return;
}

/*-------------------------------------------------------*/
/*  Update Slave DOF's                                      */
/*-------------------------------------------------------*/

void FLD::Meshtying::UpdateSlaveDOF(RCP<Epetra_Vector>&   inc)
{
  TEUCHOS_FUNC_TIME_MONITOR("Meshtying:  3.4)   - Update slave DOF");

  const Epetra_Map*  dofrowmap = discret_->DofRowMap();

  /**********************************************************************/
  /* Split inc into 3 subvectors                                       */
  /**********************************************************************/

  // container for split matrix and vector
  std::vector<RCP<Epetra_Vector> > splitvector(3);

  SplitVector(inc, splitvector);

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including meshtying)            */
  /**********************************************************************/

  RCP<LINALG::SparseMatrix> P = adaptermeshtying_.GetMortarTrafo();

  RCP<Epetra_Vector> incnew = LINALG::CreateVector(*dofrowmap,true);

  // fs: add T(mbar)*fs
  RCP<Epetra_Vector> fs_mod = rcp(new Epetra_Vector(*gsdofrowmap_,true));
  P->Multiply(false,*(splitvector[1]),*fs_mod);

  // add fn subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitvector[0]),*fnexp);
  incnew->Update(1.0,*fnexp,1.0);

  // add fn subvector to feffnew
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*(splitvector[1]),*fmexp);
  incnew->Update(1.0,*fmexp,1.0);

  // add fm subvector to feffnew
  RCP<Epetra_Vector> fs_modexp = rcp(new Epetra_Vector(*dofrowmap));
  LINALG::Export(*fs_mod,*fs_modexp);
  incnew->Update(1.0,*fs_modexp,1.0);

  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/

  inc=incnew;

  return;
}

void FLD::Meshtying::OutputSetUp()
{
  if(myrank_==0)
  {
    // Output:
    /*cout << endl << "DofRowMap:" << endl;
    cout << *(discret_->DofRowMap())<< endl << endl;
    cout << endl << "masterDofRowMap:" << endl;
    cout << *(adaptermeshtying_.MasterDofRowMap())<< endl << endl;
    cout << "slaveDofRowMap:" << endl;
    cout << *(adaptermeshtying_.SlaveDofRowMap())<< endl << endl;
    cout << "lmDofRowMap:" << endl;
    cout << *(adaptermeshtying_.LmDofRowMap())<< endl << endl;*/
    cout << "Projection matrix:" << endl;
    cout << *(adaptermeshtying_.GetMortarTrafo())<< endl << endl;
  }
}


void FLD::Meshtying::OutputSparseMatrixSplit(
    RCP<LINALG::SparseOperator>                conmat)
{
  std::vector<RCP<LINALG::SparseMatrix> > matrixsplit(9);

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

void FLD::Meshtying::OutputBlockMatrix(
    RCP<LINALG::BlockSparseMatrixBase>       blockmatrix,
    RCP<Epetra_Vector>                       residual)
{
  LINALG::SparseMatrix sysmat0 = blockmatrix->Matrix(0,0);
  LINALG::SparseMatrix sysmat1 = blockmatrix->Matrix(0,1);
  LINALG::SparseMatrix sysmat2 = blockmatrix->Matrix(0,2);

  LINALG::SparseMatrix sysmat3 = blockmatrix->Matrix(1,0);
  LINALG::SparseMatrix sysmat4 = blockmatrix->Matrix(1,1);
  LINALG::SparseMatrix sysmat5 = blockmatrix->Matrix(1,2);

  LINALG::SparseMatrix sysmat6 = blockmatrix->Matrix(2,0);
  LINALG::SparseMatrix sysmat7 = blockmatrix->Matrix(2,1);
  LINALG::SparseMatrix sysmat8 = blockmatrix->Matrix(2,2);

  cout << "Block nn" << *(sysmat0.EpetraMatrix()) << endl;
  cout << "Block nm" << *(sysmat1.EpetraMatrix()) << endl;
  cout << "Block ns" << *(sysmat2.EpetraMatrix()) << endl;

  cout << "Block mn" << *(sysmat3.EpetraMatrix()) << endl;
  cout << "Block mm" << *(sysmat4.EpetraMatrix()) << endl;
  cout << "Block ms" << *(sysmat5.EpetraMatrix()) << endl;

  cout << "Block sn" << *(sysmat6.EpetraMatrix()) << endl;
  cout << "Block sm" << *(sysmat7.EpetraMatrix()) << endl;
  cout << "Block ss" << *(sysmat8.EpetraMatrix()) << endl;

  //LINALG::PrintMatrixInMatlabFormat("sysmat_BlockMatrix",*sysmat->EpetraMatrix(),true);

  /*    if (sysmat->RowMap().SameAs(residual_->Map()))
          cout << "juhu" << endl;
        else
          cout << "nein" << endl;  */

  return;
}

void FLD::Meshtying::OutputVectorSplit(
    RCP<Epetra_Vector>                conres,
    std::vector<RCP<Epetra_Vector> >&  vectorsplit)
{
  cout << "residual " << endl << *conres << endl << endl;

  cout << "Teil fn " << endl << *(vectorsplit[0]) << endl << endl;
  cout << "Teil fm: " << endl << *(vectorsplit[1]) << endl << endl;
  cout << "Teil fs: " << endl << *(vectorsplit[2]) << endl;
  return;
}

#endif /* CCADISCRET       */
