/*----------------------------------------------------------------------------*/
/*! \file

\level 3


\brief Hybrid Additive/Multiplicative Schwarz Block Preconditioner for FSI
*/
/*----------------------------------------------------------------------------*/

// Ifpack
#include "4C_fsi_overlapprec_hybrid.hpp"

#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"  // for debugging: print matrices ToDo (mayr) remove?
#include "4C_linear_solver_preconditioner_ifpack.hpp"

#include <Ifpack_LocalFilter.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixHybridSchwarz::OverlappingBlockMatrixHybridSchwarz(
    const Core::LinAlg::MultiMapExtractor& maps, Adapter::FSIStructureWrapper& structure,
    Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale, bool structuresplit, int symmetric,
    std::vector<std::string>& blocksmoother, std::vector<double>& schuromega,
    std::vector<double>& omega, std::vector<int>& iterations, std::vector<double>& somega,
    std::vector<int>& siterations, std::vector<double>& fomega, std::vector<int>& fiterations,
    std::vector<double>& aomega, std::vector<int>& aiterations, int analyze,
    Inpar::FSI::LinearBlockSolver strategy, std::list<int> interfaceprocs,
    Inpar::FSI::Verbosity verbosity)
    : OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit, symmetric,
          omega[0], iterations[0], somega[0],
          siterations[0] - 1,  // base class counts iterations starting from 0
          fomega[0], fiterations[0] - 1, aomega[0], aiterations[0] - 1),
      strategy_(strategy),
      ifpackprec_(Teuchos::null),
      amgprec_(Teuchos::null),
      interfaceprocs_(interfaceprocs),
      additiveschwarzeverywhere_(true)
{
  if (strategy_ != Inpar::FSI::HybridSchwarz)
    FOUR_C_THROW("Type of LINEARBLOCKSOLVER parameter not recognized by this class");

  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  additiveschwarzeverywhere_ = Core::UTILS::IntegralValue<bool>(fsimono, "HYBRIDFULL");

  Inpar::FSI::LinearBlockSolver innerstrategy =
      Core::UTILS::IntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "INNERPREC");

  if (innerstrategy != Inpar::FSI::PreconditionedKrylov)
    FOUR_C_THROW("Type of INNERPREC parameter not recognized by this class");

  // create 'mulitplicative' part of hybrid preconditioner
  amgprec_ = Teuchos::rcp(new FSI::OverlappingBlockMatrixFSIAMG(maps, structure, fluid, ale,
      structuresplit, symmetric, blocksmoother, schuromega, omega, iterations, somega, siterations,
      fomega, fiterations, aomega, aiterations, analyze, innerstrategy, verbosity, this));

  return;
}

/*----------------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixHybridSchwarz::SetupPreconditioner()
{
  Teuchos::Time timer("FSI SetupPreconditioner", true);

  Teuchos::ParameterList ifpacklist;
  Teuchos::ParameterList azlist;

  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // ---------------------------------------------------------------------------
  // create 'additive' part of hybrid preconditioner
  // ---------------------------------------------------------------------------
  ifpacklist.set<int>("fact: level-of-fill", fsimono.get<int>("HYBRID_FILL_LEVEL"));
  ifpacklist.set<std::string>("amesos: solver type", "Amesos_Umfpack");

  // check whether ILU or LU (Amesos) should be used to calculate processor-local inverse
  Inpar::FSI::HybridASType type =
      Core::UTILS::IntegralValue<Inpar::FSI::HybridASType>(fsimono, "HYBRID_AS_TYPE");

  switch (type)
  {
    case Inpar::FSI::hybrid_as_type_Amesos_LU:
    {
      azlist.set<std::string>("Preconditioner Type", "Amesos");
      break;
    }
    case Inpar::FSI::hybrid_as_type_ILU:
    {
      azlist.set<std::string>("Preconditioner Type", "ILU");
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Unknown type %s for Hybrid Additive Schwarz local solver type "
          "'HYBRID_AS_TYPE'.");
      break;
    }
  }

  ifpackprec_ = Teuchos::rcp(new Core::LinearSolver::IFPACKPreconditioner(ifpacklist, azlist));

  // get blocks of system matrix and save them in 2-dim array
  std::vector<std::vector<Teuchos::RCP<Epetra_CrsMatrix>>> rows;

  const Teuchos::RCP<Epetra_CrsMatrix> Matrix00 = Matrix(0, 0).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix01 = Matrix(0, 1).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix02 = Matrix(0, 2).EpetraMatrix();
  std::vector<Teuchos::RCP<Epetra_CrsMatrix>> cols1;
  cols1.push_back(Matrix00);
  cols1.push_back(Matrix01);
  cols1.push_back(Matrix02);
  rows.push_back(cols1);

  const Teuchos::RCP<Epetra_CrsMatrix> Matrix10 = Matrix(1, 0).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix11 = Matrix(1, 1).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix12 = Matrix(1, 2).EpetraMatrix();
  std::vector<Teuchos::RCP<Epetra_CrsMatrix>> cols2;
  cols2.push_back(Matrix10);
  cols2.push_back(Matrix11);
  cols2.push_back(Matrix12);
  rows.push_back(cols2);

  const Teuchos::RCP<Epetra_CrsMatrix> Matrix20 = Matrix(2, 0).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix21 = Matrix(2, 1).EpetraMatrix();
  const Teuchos::RCP<Epetra_CrsMatrix> Matrix22 = Matrix(2, 2).EpetraMatrix();
  std::vector<Teuchos::RCP<Epetra_CrsMatrix>> cols3;
  cols3.push_back(Matrix20);
  cols3.push_back(Matrix21);
  cols3.push_back(Matrix22);
  rows.push_back(cols3);

  Epetra_Map rowmap = FullRowMap();
  const Epetra_Comm& comm = rowmap.Comm();

  /* How much memory to allocate here?
   * 3*81 works for pressure wave
   * aortic arch requires roughly 1200
   */
  /* ToDo (mayr) Estimation of storage requirement needs to be automated or
   * circumvented. Currently, just take the maximum number of entries per row
   * and multiply it with 3. Works in most cases.
   */
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, rowmap, 1200));

  bool interfaceproc = (std::find(interfaceprocs_.begin(), interfaceprocs_.end(), comm.MyPID()) !=
                        interfaceprocs_.end());

  int numRows = 3;
  int numCols = 3;

  // copy matrix block to matrix for additive Schwarz preconditioner
  for (int j = 0; j < numRows; ++j)
  {
    for (int c = 0; c < numCols; ++c)
    {
      Epetra_Map fieldmap = rows[j][c]->RowMap();
      int numMyElements = fieldmap.NumMyElements();
      int maxNumEntries = rows[j][c]->MaxNumEntries();

      std::vector<double> values(maxNumEntries);
      std::vector<int> indices(maxNumEntries);

      double one = 1.0;

      for (int i = 0; i < numMyElements; ++i)
      {
        int numEntries = 0;
        if (additiveschwarzeverywhere_ or interfaceproc)
        {
          int err = rows[j][c]->ExtractGlobalRowCopy(
              fieldmap.GID(i), maxNumEntries, numEntries, values.data(), indices.data());
          if (err != 0) FOUR_C_THROW("ExtractGlobalRowCopy failed, error = %d!", err);
        }
        else
        {
          if (j == c)
          {
            numEntries = 1;
            values[0] = one;
            indices[0] = fieldmap.GID(i);
          }
        }
        if (numEntries > 0)
        {
          int err =
              A->InsertGlobalValues(fieldmap.GID(i), numEntries, values.data(), indices.data());
          if (err != 0)
            FOUR_C_THROW(
                "InsertGlobalValues failed, error = %d!.\n"
                "Available number of row entries is %d.",
                err, numEntries);
        }
      }
    }
  }

  A->FillComplete();
  A->OptimizeStorage();

  comm.Barrier();
  if (comm.MyPID() == 0)
    std::cout << "Copied matrix in " << timer.totalElapsedTime(true) << " seconds." << std::endl;
  timer.reset();

  /****************************************************************************/


  /****************************************************************************/

  // Core::LinAlg::PrintMatrixInMatlabFormat("precondmat.dat",*A);

  // Teuchos::RCP<Ifpack_LocalFilter> localmat = Teuchos::rcp(new Ifpack_LocalFilter(A));

  Teuchos::RCP<Epetra_MultiVector> x = Teuchos::rcp(new Epetra_MultiVector(rowmap, 1));
  ifpackprec_->setup(true, A.getRawPtr(), x.getRawPtr(), x.getRawPtr());

  comm.Barrier();
  if (comm.MyPID() == 0)
    std::cout << "Built ILU in " << timer.totalElapsedTime(true) << " seconds" << std::endl;
  timer.reset();

  // setup 'multiplicative' part of hybrid preconditioner
  amgprec_->SetupPreconditioner();

  comm.Barrier();
  if (comm.MyPID() == 0) printf("AMG prec in %f seconds.\n", timer.totalElapsedTime(true));
  timer.reset();

  //  /*
  //   * Test if we build the right prec
  //   */
  //  int numproc = comm.NumProc();
  //  int myrank = comm.MyPID();
  //
  //  int numGlobalElements = numproc * 2;
  //  int numMyElements = 2;
  //
  //  int myrows[2];
  //  myrows[0] = myrank * 2;
  //  myrows[1] = myrank * 2 + 1;
  //
  //  const int* rowptr = myrows;
  //
  //  Teuchos::RCP<Epetra_Map> testmap = Teuchos::rcp(new
  //  Epetra_Map(numGlobalElements,numMyElements,rowptr,0,comm));
  //
  //  Teuchos::RCP<Epetra_CrsMatrix> testmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*testmap,0));
  //
  //  int indices[2*numproc];
  //
  //  double values[2*numproc];
  //
  //  for (int i=0; i<2*numproc; ++i){
  //    if (i==myrank*2) // diagonal
  //      values[i] = 2;
  //    else
  //      values[i] = 1;
  //
  //    indices[i] = i;
  //
  //  }
  //
  //  int err = testmat->InsertGlobalValues(myrank * 2, numproc * 2, values, indices);
  //  if (err != 0)
  //    FOUR_C_THROW("InsertGlobalValues failed, error = %d!",err);
  //
  //  for (int i=0; i<2*numproc; ++i){
  //    if (i==myrank*2+1) // diagonal
  //      values[i] = 2;
  //    else
  //      values[i] = 1;
  //
  //    indices[i] = i;
  //
  //  }
  //
  //  err = testmat->InsertGlobalValues(myrank * 2 + 1, numproc * 2, values, indices);
  //  if (err != 0)
  //    FOUR_C_THROW("InsertGlobalValues failed, error = %d!",err);
  //
  //  testmat->fill_complete();
  //  testmat->OptimizeStorage();
  //
  //  Teuchos::RCP<Core::LinAlg::SOLVER::IFPACKPreconditioner> testprec;
  //  testprec = Teuchos::rcp(new
  //  Core::LinAlg::SOLVER::IFPACKPreconditioner(&outfile,ifpacklist,azlist));
  //  Teuchos::RCP<Epetra_MultiVector> xtest = Teuchos::rcp(new Epetra_MultiVector(*testmap,1));
  //  testprec->setup(true, testmat.getRawPtr(), x.getRawPtr(), x.getRawPtr());
  //
  //
  //  Core::LinAlg::PrintMatrixInMatlabFormat("testmat.dat",*testmat);
  //  std::cout<<"\nEpetraMatrix: "<<std::endl;
  //  testmat->print(std::cout);
  //
  //  comm.Barrier();
  //

  //  Teuchos::RCP<Epetra_RowMatrix> Pmatrix = Teuchos::rcp(new Epetra_CrsMatrix(*A));
  //
  //  Ifpack Factory;
  //  directifpackprec_ = Teuchos::rcp( Factory.Create("Amesos",localmat.get(),0) );
  //  int err = directifpackprec_->SetParameters(ifpacklist);
  //  std::cout<<"\nparamserr: "<<err;
  //  if (err != 0)
  //    FOUR_C_THROW("\nSet parameters failed, error code %d",err);
  //  err = directifpackprec_->Initialize();
  //  std::cout<<"\niniterr: "<<err;
  //  if (err != 0)
  //    FOUR_C_THROW("\nInitialization of preconditioner failed, error code %d",err);
  //  err = directifpackprec_->Compute();
  //  std::cout<<"\ncomperr: "<<err;
  //  if (err != 0)
  //    FOUR_C_THROW("\nComputation of preconditioner failed, error code %d",err);

  return;
}

/*----------------------------------------------------------------------------*/
int FSI::OverlappingBlockMatrixHybridSchwarz::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(Copy, X, 0));
  Teuchos::RCP<Epetra_Vector> z = Teuchos::rcp(new Epetra_Vector(Copy, X, 0));

  z->Update(0.0, *b, 0.0);

  // apply 'additive' part
  int err = ifpackprec_->PrecOperator()->ApplyInverse(*b, *z);
  if (err != 0) FOUR_C_THROW("Preconditioning 1 failed.");

  Teuchos::RCP<Epetra_Vector> Az = Teuchos::rcp(new Epetra_Vector(Copy, X, 0));
  Az->Update(0.0, *b, 0.0);

  err = this->Apply(*z, *Az);
  if (err != 0) FOUR_C_THROW("Matrix 1 failed.");

  Teuchos::RCP<Epetra_Vector> tmpb = Teuchos::rcp(new Epetra_Vector(Copy, *b, 0));
  tmpb->Update(-1.0, *Az, 1.0);

  Teuchos::RCP<Epetra_Vector> tmpz = Teuchos::rcp(new Epetra_Vector(Copy, X, 0));
  tmpz->Update(0.0, *b, 0.0);

  // apply 'multiplicative' part
  err = amgprec_->ApplyInverse(*tmpb, *tmpz);
  if (err != 0) FOUR_C_THROW("Preconditioning 2 failed.");

  z->Update(1.0, *tmpz, 1.0);

  Az->Update(0.0, *b, 0.0);
  err = this->Apply(*z, *Az);
  if (err != 0) FOUR_C_THROW("Matrix 2 failed.");

  tmpb->Update(1.0, *b, 0.0);
  tmpb->Update(-1.0, *Az, 1.0);
  tmpz->Update(0.0, *b, 0.0);

  // apply 'additive' part
  err = ifpackprec_->PrecOperator()->ApplyInverse(*tmpb, *tmpz);
  if (err != 0) FOUR_C_THROW("Preconditioning 3 failed.");

  z->Update(1.0, *tmpz, 1.0);

  Y.Update(1.0, *z, 0.0);

  return 0;
}

/*----------------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrixHybridSchwarz::Label() const
{
  return "Unknown strategy in FSI::OverlappingBlockMatrixHybridSchwarz";
}

FOUR_C_NAMESPACE_CLOSE
