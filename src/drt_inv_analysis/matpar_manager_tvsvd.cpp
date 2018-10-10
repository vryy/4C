/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_tvsvd.cpp
\brief Creating reduced basis from TV approximation
<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_tvsvd.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

// Some SIGN function in Teuchos_BLAS.hpp clashes with
// our definitions.h macro SIGN. Somehow via pss_table
// it enters the header hierarchy
#ifdef SIGN
#undef SIGN
#endif

/* Anasazi headers*/
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include <fenv.h>

typedef std::map<int, std::vector<int>> PATCHES;

/*----------------------------------------------------------------------*/
INVANA::MatParManagerTVSVD::MatParManagerTVSVD(Teuchos::RCP<DRT::Discretization> discret)
    : MatParManagerPerElement(discret),
      max_num_levels_(1),
      seed_(1),
      eps_(1.0e-02),
      map_restart_file_("none"),
      map_restart_step_(-1),
      nev_(10),
      nblocks_(10),
      bsize_(4)
{
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::Setup()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "MatParManagerTVSVD Setup:" << std::endl;
  }

  // MAP restart stuff
  const Teuchos::ParameterList& invp = Inpar();
  map_restart_step_ = invp.get<int>("MAP_REDUCT_RESTART");
  map_restart_file_ = Teuchos::getNumericStringParameter(invp, "MAP_REDUCT_RESTARTFILE");

  // eigendecomposition stuff
  max_num_levels_ = invp.get<int>("NUM_REDUCT_LEVELS");
  nev_ = invp.get<int>("TVSVD_ANASAZI_NEV");
  nblocks_ = invp.get<int>("TVSVD_ANASAZI_NBLOCKS");
  bsize_ = invp.get<int>("TVSVD_ANASAZI_BSIZE");

  eps_ = invp.get<double>("TVD_EPS");

  if (max_num_levels_ < 1) dserror("Choose at least NUM_LEVELS = 1 for the basis creation!");

  // call setup of the Base class to have all the
  // layout of the elementwise distribution
  MatParManagerPerElement::Setup();
  optparams_elewise_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_, 1, true));
  elewise_map_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // set up random number generator consistently in case of nested parallelity
  util_.SetSeed(seed_ + Comm().MyPID());

  // read map approximation as evaluation point for the TV linearization
  ReadMAPApproximation();


  // create sparse approximation of the MAP solution
  CreateProjection();

  // initialize parameters
  InitParameters();

  // Some user information
  if (Comm().MyPID() == 0) std::cout << std::endl;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // Inject into the elementwise solution space
  Epetra_MultiVector projected(*elewise_map_, 1, true);
  int err = projector_->Multiply(true, *optparams_, projected);
  if (err != 0) dserror("Application of prolongator failed.");

  projected.Update(1.0, *optparams_elewise_, 1.0);

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(*projected(0), k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      params->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid], k, (*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::ApplyParametrization(
    DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals)
{
  // this is ok here since we have a sparse approximation
  Teuchos::RCP<Epetra_CrsMatrix> fullmatrix = matrix.FillMatrix();

  // todo: this is not ok! loop over the single columns of matrix
  // and extract only the diagonal component.
  // matrix * projector_
  Teuchos::RCP<Epetra_CrsMatrix> mr = LINALG::Multiply(fullmatrix, false, projector_, false);
  // projector_'*matrix*projector_
  Teuchos::RCP<Epetra_CrsMatrix> pmr = LINALG::Multiply(projector_, true, mr, false);

  Epetra_Vector diagonal(pmr->RowMap(), true);
  pmr->ExtractDiagonalCopy(diagonal);

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(diagonal, k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      diagonals->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid], k, (*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::InitParameters()
{
  // sanity checks
  if (not projector_->DomainMap().PointSameAs(optparams_elewise_->Map()))
    dserror("Restrictor->DomainMap error.");

  if (not projector_->RangeMap().PointSameAs(optparams_->Map()))
    dserror("Restrictor->RangeMap error");

  // optparams are initialized to zero
  optparams_initial_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
    double val, int elepos, int paraposglobal, int paraposlocal)
{
  if (EleGIDtoparamsLID().find(elepos) == EleGIDtoparamsLID().end())
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  // parameter in the 'elementwise' optimization parameter layout
  int plid = EleGIDtoparamsLID()[elepos].at(paraposlocal);

  // when the algorithm comes in here from the MatParManager::AddEvaluate
  // it comes with elementwise parameters processed locally, so for the
  // 'nachdifferenzieren' one can just pick the local entry plid of each
  // eigenvector
  // check to be sure
  int found = evecs_->Map().GID(plid);
  if (found == -1) dserror("ID not found on this proc. This is fatal!");

  for (int i = 0; i < paramlayoutmap_->NumMyElements(); i++)
  {
    double ival = (*(*evecs_)(i))[plid] * val;
    int success = dfint->SumIntoMyValue(i, 0, ival);
    if (success != 0) dserror("Summation into gradient resulted in %d", success);
  }
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::Finalize(
    Teuchos::RCP<Epetra_MultiVector> source, Teuchos::RCP<Epetra_MultiVector> target)
{
  // sum across processor
  std::vector<double> val(source->MyLength(), 0.0);
  Discret()->Comm().SumAll((*source)(0)->Values(), &val[0], source->MyLength());

  for (int i = 0; i < target->MyLength(); i++) target->SumIntoGlobalValue(i, 0, val[i]);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::ReadMAPApproximation()
{
  // Create the input control file object
  Teuchos::RCP<IO::InputControl> input =
      Teuchos::rcp(new IO::InputControl(map_restart_file_, Discret()->Comm()));

  // and the discretization reader to read from the input file
  IO::DiscretizationReader reader(Discret(), input, map_restart_step_);

  if (not Discret()->Comm().MyPID())
  {
    std::cout << "  Reading MAP approximation: ";
    std::cout << "  step " << map_restart_step_ << " (from: " << input->FileName() << ")"
              << std::endl;
  }

  // MAP solution
  reader.ReadMultiVector(optparams_elewise_, "solution");

  // Read lbfgs matrix storage
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, elewise_map_.get(), true));
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, elewise_map_.get(), true));

  int actsize = reader.ReadInt("storage_size");

  // initialize storage
  sstore->Resize(-actsize + 1, 0, elewise_map_.get(), true);
  ystore->Resize(-actsize + 1, 0, elewise_map_.get(), true);

  Teuchos::RCP<Epetra_MultiVector> storage =
      Teuchos::rcp(new Epetra_MultiVector(*elewise_map_, actsize, false));

  reader.ReadMultiVector(storage, "sstore");
  for (int i = 0; i < actsize; i++) sstore->UpdateSteps(*(*storage)(i));

  storage->Scale(0.0);
  reader.ReadMultiVector(storage, "ystore");
  for (int i = 0; i < actsize; i++) ystore->UpdateSteps(*(*storage)(i));

  // ---- create covariance matrix
  double scalefac = 1.0;
  bool objfuncscal = DRT::INPUT::IntegralValue<bool>(Inpar(), "OBJECTIVEFUNCTSCAL");
  if (objfuncscal) scalefac = Objfunct().GetScaleFac();

  bool initscal = DRT::INPUT::IntegralValue<bool>(Inpar(), "LBFGSINITSCAL");

  fullcovariance_ = Teuchos::rcp(new DcsMatrix(sstore, ystore, initscal, objfuncscal, scalefac));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::CreateProjection()
{
  // Quadratic approximation of the TV functional
  SetupTVOperator();

  // hack
  // LINALG::PrintMatrixInMatlabFormat("lintvop",*lintvop_);

  // Factorization of the linear operator
  Factorize();

  // Set up projection operators
  SetupRandP(max_num_levels_);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::SetupRandP(int numvecs)
{
  int numvecspp = 0;
  if (Comm().MyPID() == 0) numvecspp = numvecs;

  // build paramlayout maps anew
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1, numvecspp, 0, Comm()));
  paramlayoutmap_ = LINALG::AllreduceEMap(*paramlayoutmapunique_);

  Teuchos::RCP<Epetra_Map> colmap = LINALG::AllreduceEMap(graph_->RowMap(), 0);
  int maxbw = colmap->NumGlobalElements();
  projector_ =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *paramlayoutmapunique_, *colmap, maxbw, false));

  // get all the eigenvectors to proc 0
  Epetra_MultiVector evecs(*colmap, nev_, true);
  LINALG::Export(*evecs_, evecs);

  // make unique eigenvectors in all groups
  // BroadcastEigenvectors(evecs);
  MakeConsistentEigenvectors(evecs);

  // insert (unit-) eigenvectors into projection
  for (int i = 0; i < projector_->NumMyRows(); i++)
  {
    int numentries = colmap->NumGlobalElements();
    Epetra_Vector* row = evecs(i);

    int err =
        projector_->InsertGlobalValues(i, numentries, row->Values(), colmap->MyGlobalElements());
    if (err < 0) dserror("Restrictor/Prolongator insertion failed.");
  }

  // Fill Complete
  int err = projector_->FillComplete(graph_->RowMap(), *paramlayoutmapunique_, true);
  if (err != 0) dserror("Restrictor/Prolongator FillComplete failed.");

  // initialize optimization parameters to zero
  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_, 1, true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_, 1, true));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::SetupTVOperator()
{
  if (optparams_elewise_->NumVectors() != 1)
    dserror(
        "Only optimization parameters with NumVectors()==1"
        "are supported in the inverse analysis");

  // -------------------------------------------------------------------
  // Majorization at current map solution:
  Epetra_Vector u(graph_->RowMap(), true);

  // zero out diagonal of the graph
  graph_->ReplaceDiagonalValues(u);

  // initialize the linearized tv operator
  lintvop_ = Teuchos::rcp(new Epetra_CrsMatrix(*graph_));
  lintvop_->FillComplete();
  Teuchos::RCP<Epetra_CrsMatrix> lintvop2 = Teuchos::rcp(new Epetra_CrsMatrix(*graph_));
  lintvop2->FillComplete();

  // communicate theta data from other procs such that every proc can compute sums
  // over adjacent parameters
  Epetra_MultiVector thetacol(graph_->ColMap(), 1, false);
  LINALG::Export(*optparams_elewise_, thetacol);

  for (int i = 0; i < (*optparams_elewise_)(0)->MyLength(); i++)
  {
    // get weights of neighbouring parameters
    int lenindices = graph_->NumMyEntries(i);
    int numindex;
    std::vector<int> indices(lenindices, 0);
    std::vector<double> weights(lenindices, 0);
    graph_->ExtractMyRowCopy(i, lenindices, numindex, &weights[0], &indices[0]);

    // row in local index space of the collayout
    int rowi = thetacol.Map().LID(optparams_elewise_->Map().GID(i));

    double rowsum = 0.0;
    double rowval = (*(*optparams_elewise_)(0))[i];
    for (int j = 0; j < lenindices; j++)
    {
      if (indices[j] != rowi)  // skip substracting from itself
      {
        double colval = (*thetacol(0))[indices[j]];
        rowsum += weights[j] * (colval - rowval) * (colval - rowval);
      }
    }
    // sum over all the rows
    u[i] = 1.0 / sqrt(rowsum + eps_);
  }

  // -------------------------------------------------------------------
  // contributions from nominator i
  // put row sums of the graph on the diagonal
  Epetra_Vector diagsum(optparams_elewise_->Map(), true);
  Epetra_Vector ones(optparams_elewise_->Map(), false);
  ones.PutScalar(1.0);
  int err = graph_->Multiply(false, ones, diagsum);
  if (err != 0) dserror("Matrix-Vector multiplication failed");
  diagsum.Scale(-1.0);
  // subtract ( there should be nothing on the diagonal sofar!)
  lintvop_->ReplaceDiagonalValues(diagsum);
  // scale row_i with approximation u_i
  err = lintvop_->LeftScale(u);
  if (err != 0) dserror("Matrix left scale failed");

  // -------------------------------------------------------------------
  // contributions from nominator j
  err = lintvop2->RightScale(u);
  if (err != 0) dserror("Matrix right scale failed");
  lintvop2->Scale(-1.0);
  lintvop2->InvRowSums(diagsum);  // summing absolute values is exactly correct
  // we needed sums no inverse sums
  for (int i = 0; i < diagsum.MyLength(); i++) diagsum[i] = 1.0 / diagsum[i];
  lintvop2->ReplaceDiagonalValues(diagsum);


  // Add up contributions
  LINALG::Add(lintvop2, false, 1.0, lintvop_, -1.0);

  // regularization helps anasazi
  Epetra_Vector newdiag(optparams_elewise_->Map(), false);
  newdiag.PutScalar(0.00001);
  lintvop_->ExtractDiagonalCopy(diagsum);
  diagsum.Update(1.0, newdiag, 1.0);
  lintvop_->ReplaceDiagonalValues(diagsum);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::Factorize()
{
  //------------------------------------------------
  // Setup the eigenproblem using Anasazi
  evecs_ = Teuchos::rcp(new Epetra_MultiVector(graph_->RowMap(), nev_));

  // prerequ: (no test here since Anasazi complains)
  // numblock*blocksize + maxlocked must be < spacedim
  // maxlocked + blocksize > nev

  Teuchos::ParameterList params;
  std::string which("SM");
  params.set("Which", which);
  params.set("nev", nev_);
  params.set("Block Size", bsize_);
  params.set("Num Blocks", nblocks_);
  params.set("Maximum Restarts", 100);
  params.set("Convergence Tolerance", 1.0e-6);
  params.set("Use Locking", true);
  params.set("Relative Convergence Tolerance", false);
  params.set("Verbosity", Anasazi::Errors);

  // the following seems a little sensible and occasionally fails
  // with floating points exceptions --> Catch them
#ifdef TRAP_FE
  fedisableexcept(FE_ALL_EXCEPT);
#endif

  int iraised = 0;
  int graised = 1;
  while (graised)
  {
    feclearexcept(FE_ALL_EXCEPT);

    try
    {
      AnasaziEigenProblem(lintvop_, evecs_, params);
      // AnasaziEigenProblem(fullcovariance_->FillMatrix(),evecs_,params);

      iraised = fetestexcept(FE_INVALID | FE_DIVBYZERO);

      Comm().SumAll(&iraised, &graised, 1);

      if (graised)
        throw std::runtime_error("AnasaziEigenProblem failed with exception. Try again.");
    }
    catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      graised = 1;
      // Reset seed, since otherwise this group will maybe get evecs_*(-1)
      util_.SetSeed(seed_ + Comm().MyPID());
      Comm().Barrier();
    }
  }

#ifdef TRAP_FE
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::AnasaziEigenProblem(Teuchos::RCP<Epetra_CrsMatrix> A,
    Teuchos::RCP<Epetra_MultiVector> L, Teuchos::ParameterList params)
{
  Anasazi::BasicOutputManager<double> printer;
  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  const int nev = params.get("nev", 1);
  const int blockSize = params.get("Block Size", nev);
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;


  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note: This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp(new Epetra_MultiVector(L->Map(), blockSize));
  Random(*ivec);

  //------------------------------------------------
  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP>> problem =
      Teuchos::rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec));

  problem->setHermitian(true);
  problem->setNEV(nev);

  // done
  bool boolret = problem->setProblem();
  if (boolret != true)
  {
    printer.print(Anasazi::Errors, "Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
    dserror("Anasazi could finalize the problem setup");
  }

  //------------------------------------------------
  // Solve
  Anasazi::BlockDavidsonSolMgr<double, MV, OP> solverman(problem, params);
  Anasazi::ReturnType returnCode = solverman.solve();
  if (returnCode != Anasazi::Converged) dserror("Anasazi didn't converge finding an eigenbasis");

  //------------------------------------------------
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double, MV> sol = problem->getSolution();
  std::vector<Anasazi::Value<double>> evals = sol.Evals;
  L->Update(1.0, *sol.Evecs, 0.0);

  //------------------------------------------------
  // Compute residuals.
  std::vector<double> normR(sol.numVecs);
  if (sol.numVecs > 0)
  {
    Teuchos::SerialDenseMatrix<int, double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector tempAevec(L->Map(), sol.numVecs);
    T.putScalar(0.0);
    for (int i = 0; i < sol.numVecs; i++)
    {
      T(i, i) = evals[i].realpart;
    }
    A->Apply(*L, tempAevec);
    MVT::MvTimesMatAddMv(-1.0, *L, T, 1.0, tempAevec);
    MVT::MvNorm(tempAevec, normR);
  }

  // Print the results
  std::ostringstream os;
  os.setf(std::ios_base::right, std::ios_base::adjustfield);
  os << "Solver manager returned "
     << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << std::endl;
  os << std::endl;
  os << "------------------------------------------------------" << std::endl;
  os << std::setw(16) << "Eigenvalue" << std::setw(18) << "Direct Residual" << std::endl;
  os << "------------------------------------------------------" << std::endl;
  for (int i = 0; i < sol.numVecs; i++)
  {
    os << std::setw(16) << evals[i].realpart << std::setw(18) << normR[i] / evals[i].realpart
       << std::endl;
  }
  os << "------------------------------------------------------" << std::endl;
  printer.print(Anasazi::Errors, os.str());
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::Random(Epetra_MultiVector& randvec)
{
  const int mylength = randvec.MyLength();
  const int numvecs = randvec.NumVectors();
  double** pointers = randvec.Pointers();

  // use particularly seeded Epetra_Util object to fill at random
  for (int i = 0; i < numvecs; i++)
  {
    double* const to = pointers[i];
    for (int j = 0; j < mylength; j++) to[j] = util_.RandomDouble();
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::BroadcastEigenvectors(Epetra_MultiVector& evecs)
{
  // Prerequisite: Only data from global proc 0 to all other local proc 0s will
  // be performed

  // get the necessary communicators
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetNPGroup()->GlobalComm();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();

  // pack the Multivector on global proc0
  std::map<int, std::vector<double>> data;
  if (gcomm->MyPID() == 0)
  {
    const int numvecs = evecs.NumVectors();
    const int mylength = evecs.MyLength();
    double** pointers = evecs.Pointers();

    for (int i = 0; i < numvecs; i++)
    {
      double* const val = pointers[i];
      std::vector<double> dummy(mylength, 0.0);
      data[i] = dummy;
      for (int j = 0; j < mylength; j++) data[i][j] = val[j];
    }
  }

  // create source map
  int sourcegids = 0;
  if (gcomm->MyPID() == 0) sourcegids = evecs.NumVectors();
  Epetra_Map sourcemap(-1, sourcegids, 0, *gcomm);

  // create target map
  int numtargetgids = 0;
  std::vector<int> targetgids;
  if (lcomm->MyPID() == 0)
  {
    numtargetgids = evecs.NumVectors();
    targetgids.resize(numtargetgids);
    for (int i = 0; i < numtargetgids; i++) targetgids[i] = i;
  }
  Epetra_Map targetmap(-1, numtargetgids, &targetgids[0], 0, *gcomm);

  // create exporter for the data and export
  DRT::Exporter ex(sourcemap, targetmap, *gcomm);
  ex.Export(data);

  // zero out target to be sure
  evecs.PutScalar(0.0);

  // Unpack on every local proc 0
  if (lcomm->MyPID() == 0)
  {
    const int numvecs = evecs.NumVectors();
    const int mylength = evecs.MyLength();
    double** pointers = evecs.Pointers();

    for (int i = 0; i < numvecs; i++)
    {
      double* const val = pointers[i];
      for (int j = 0; j < mylength; j++) val[j] = data[i][j];
    }
  }
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerTVSVD::MakeConsistentEigenvectors(Epetra_MultiVector& evecs)
{
  // prerequesite: Evecs must live on local proc 0 already

  Epetra_MultiVector abs_evecs(evecs);
  abs_evecs.Abs(evecs);

  const int numvecs = abs_evecs.NumVectors();
  const int length = abs_evecs.MyLength();

  double** pointers = abs_evecs.Pointers();

  for (int i = 0; i < numvecs; i++)
  {
    // find max abs value
    int maxi = -1;
    double maxv = 0.0;
    double* const val = pointers[i];
    for (int j = 0; j < length; j++)
    {
      if (val[j] > maxv)
      {
        maxi = j;
        maxv = val[j];
      }
    }

    // correct signs in case
    if (maxi != -1)
    {
      double mval = evecs[i][maxi];
      if (mval < 0.0) evecs(i)->Scale(-1.0);
    }
  }
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::MatParManagerTVSVD::InitialCovariance()
{
  // projector * covariance
  Teuchos::RCP<Epetra_CrsMatrix> interm = Teuchos::rcp(new Epetra_CrsMatrix(
      Copy, projector_->RowMap(), projector_->ColMap(), projector_->MaxNumEntries(), false));

  Teuchos::RCP<Epetra_Vector> column =
      Teuchos::rcp(new Epetra_Vector(fullcovariance_->RowMap(), false));
  Teuchos::RCP<Epetra_Vector> col_interm =
      Teuchos::rcp(new Epetra_Vector(projector_->RowMap(), false));
  int maxnumentries = fullcovariance_->MaxNumEntries();
  for (int i = 0; i < maxnumentries; i++)
  {
    // get column
    fullcovariance_->ExtractGlobalColumnCopy(i, *column);

    // project this column
    int err = projector_->Multiply(false, *column, *col_interm);
    if (err != 0) dserror("Projection failed.");

    // put column into intermediate matrix
    double* vals;
    int colind = i;
    col_interm->ExtractView(&vals);
    for (int j = 0; j < col_interm->MyLength(); j++)
    {
      int gid = projector_->RowMap().GID(j);
      interm->InsertGlobalValues(gid, 1, &vals[j], &colind);
    }
  }
  interm->FillComplete(projector_->ColMap(), projector_->RangeMap());

  // interm * projector'
  Teuchos::RCP<Epetra_CrsMatrix> cov = LINALG::Multiply(*interm, false, *projector_, true);
  return cov;
}
