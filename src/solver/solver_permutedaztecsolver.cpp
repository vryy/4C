/*
 * solver_permutedaztecsolver.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: tobias
 */

#ifdef HAVE_MueLu
#ifdef HAVE_EXPERIMENTAL_MueLu

#include <MueLu_ConfigDefs.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <MueLu.hpp>
#include <MueLu_FactoryBase.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver.hpp>    // remove me
#include <MueLu_HierarchyHelpers.hpp>
#include <MueLu_VerboseObject.hpp>

// Aztec headers
#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"

#include "Epetra_MultiVector.h"

// BACI headers
#include "../drt_lib/drt_dserror.H"
#include "solver_permutedaztecsolver.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::PermutedAztecSolver::PermutedAztecSolver( const Epetra_Comm & comm,
                                            Teuchos::ParameterList & params,
                                            FILE * outfile )
  : KrylovSolver(comm,params,outfile)
{
  bPermuteLinearSystem_ = true; // always permute linear system // TODO introduce decision
  ncall_ = 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::PermutedAztecSolver::~PermutedAztecSolver()
{
  data_ = Teuchos::null;
  PermFact_ = Teuchos::null;
  preconditioner_ = Teuchos::null;
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::PermutedAztecSolver::Setup( Teuchos::RCP<Epetra_Operator> matrix,
                                          Teuchos::RCP<Epetra_MultiVector> x,
                                          Teuchos::RCP<Epetra_MultiVector> b,
                                          bool refactor,
                                          bool reset,
                                          Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean,
                                          Teuchos::RCP<Epetra_MultiVector> kernel_c,
                                          bool project)
{
#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup_ = 0.;
  Epetra_Time tttcreate(Comm()); // time measurement for creation of preconditioner
#endif

  if (!Params().isSublist("Aztec Parameters"))
    dserror("Do not have aztec parameter list");
  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");
  //int azoutput = azlist.get<int>("AZ_output",0);

  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>( matrix );

  // build
  // permP, permQT and A = permQ^T A permP
  // all variables and information is stored in data_
  BuildPermutationOperator(A);

  // TODO decide whether to permute linear system or not using the information of
  //      the permuted system matrix A

  // find problematic rows/columns which are not optimal after permutation
  Teuchos::RCP<Map> nonDiagMap = FindNonDiagonalDominantRows(0.5); // introduce solver parameter for this??

  Teuchos::ParameterList & precondParams = getPrecondParamterList();
  Teuchos::ParameterList & linSystemProps = precondParams.sublist("Linear System properties");
  linSystemProps.set<Teuchos::RCP<Map> >("non diagonal-dominant row map",nonDiagMap);

  if(bPermuteLinearSystem_) {
    // set
    // b_ = permP * b;
    // A_ = permQ^T * A * permP
    PermuteLinearSystem(A,b);

    // calculate (permQT)^T * b_f where b_f is the fine level null space (multi)vector
    //PermuteNullSpace(A);  // TODO think about this
    // do not permute null space to preserve pattern of null space for transfer operators
    // important e.g. for one pt aggregates?
  } else {
    b_ = b;
    A_ = A;
  }
  x_ = x;
  ////

#ifdef WRITEOUTSTATISTICS
  tttcreate.ResetStartTime();
#endif

  // decide whether we recreate preconditioners
  int  reuse  = azlist.get("reuse",0);
  bool create = reset or not Ncall() or not reuse or ( Ncall() % reuse )==0;
  if ( create )
  {
    ncall_ = 0;
    CreatePreconditioner( azlist, A!=Teuchos::null, weighted_basis_mean, kernel_c, project );
  }

  preconditioner_->Setup( create, &*A_, &*x_, &*b_ );

#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup = tttcreate.ElapsedTime();
#endif
}

void LINALG::SOLVER::PermutedAztecSolver::BuildPermutationOperator(const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xCrsA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(Params().sublist("NodalBlockInformation").get<int>("nv")); // set nBlockSize

  if(data_==Teuchos::null) data_ = Teuchos::rcp(new Level());
  data_->setDefaultVerbLevel(Teuchos::VERB_NONE);
  data_->Set("A",xOp);

  // check, if "SlaveDofMap" information is available in parameter lists
  bool bSlaveDofMapAvailable = false;
  if(Params().isSublist("Aztec Parameters") &&
      Params().sublist("Aztec Parameters").isParameter("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap")) {
    bSlaveDofMapAvailable = true;
  }

  // build permutation factory
  if(bSlaveDofMapAvailable) {
    Teuchos::RCP<Epetra_Map> epSlaveDofMap  = Params().sublist("Aztec Parameters").get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
    Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
    data_->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xSlaveDofMap));  // set map with active dofs

    // define permutation factory for permuting the full matrix A
    PermFact_ = Teuchos::rcp(new PermutationFactory("SlaveDofMap", MueLu::NoFactory::getRCP()));
  }
  else
    // permute full matrix
    PermFact_ = Teuchos::rcp(new PermutationFactory("",Teuchos::null));

  // setup main factory manager
  Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
  M->SetFactory("permQT",          PermFact_);
  M->SetFactory("A",               MueLu::NoFactory::getRCP()); // this is the input matrix
  MueLu::SetFactoryManager SFMFinest(data_, M); // set factory manager for data container

  // prepare building process for permutation operators
  data_->Request("A", PermFact_.get());
  data_->Request("permA", PermFact_.get());
  data_->Request("permP", PermFact_.get());
  data_->Request("permQT", PermFact_.get());
  data_->Request("permScaling", PermFact_.get());
  data_->Request("#RowPermutations", PermFact_.get());
  data_->Request("#ColPermutations", PermFact_.get());
  data_->Request("#WideRangeRowPermutations", PermFact_.get());
  data_->Request("#WideRangeColPermutations", PermFact_.get());

  // build permutation operators
  PermFact_->Build(*data_);

}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::PermutedAztecSolver::PermuteLinearSystem(const Teuchos::RCP<Epetra_CrsMatrix>& A,const Teuchos::RCP<Epetra_MultiVector>& b) {

  if(!data_->IsAvailable("A", PermFact_.get()) ||
      !data_->IsAvailable("permP", PermFact_.get()) ||
      !data_->IsAvailable("permScaling", PermFact_.get()))
    dserror("PermutedAztecSolver: call BuildPermutationOperator before PermuteLinearSystem");

  // extract permutation operators from data_
  Teuchos::RCP<Epetra_CrsMatrix> xEpPermCrsMat = GetOperatorNonConst("A",PermFact_);
  Teuchos::RCP<const Epetra_CrsMatrix> epPermPMatrix  = GetOperator("permP",  PermFact_);             // row permutation matrix
  Teuchos::RCP<const Epetra_CrsMatrix> epPermScalingMatrix = GetOperator("permScaling",PermFact_); // leftScaling matrix

  // P_trafo*b
  Teuchos::RCP<Epetra_MultiVector> btemp1 = Teuchos::rcp(new Epetra_MultiVector(*b));
  epPermPMatrix->Multiply(false, *b, *btemp1);
  // P_scaling * P_trafo * b
  epPermScalingMatrix->Multiply(false, *btemp1, *b);

  // set
  // b_ = permP * b;
  // A_ = permQ^T * A * permP
  b_ = b;   // note b is permuted
  A_ = xEpPermCrsMat;//xEpPermCrsMat->getEpetra_CrsMatrixNonConst();

}

void LINALG::SOLVER::PermutedAztecSolver::PermuteNullSpace(const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xCrsA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(Params().sublist("NodalBlockInformation").get<int>("nv")); // set nBlockSize

  // detect MueLu/ML Paramter list
  std::string MultiGridParameterListName = "";
  if(Params().isSublist("ML Parameters"))
    MultiGridParameterListName = "ML Parameters";
  else if(Params().isSublist("MueLu Parameters"))
    MultiGridParameterListName = "MueLu Parameters";
  else if(Params().isSublist("MueLu (Contact) Parameters"))
    MultiGridParameterListName = "MueLu (Contact) Parameters";
  else if(Params().isSublist("MueLu (Contact2) Parameters"))
    MultiGridParameterListName = "MueLu (Contact2) Parameters";
  else if(Params().isSublist("MueLu (Contact3) Parameters"))
    MultiGridParameterListName = "MueLu (Contact3) Parameters";
  else if(Params().isSublist("MueLu (PenaltyContact) Parameters"))
    MultiGridParameterListName = "MueLu (PenaltyContact) Parameters";

  // retransform nullspace vectors
  if(MultiGridParameterListName == "")
    return; // no nullspace to permute

  Teuchos::RCP<Matrix> xPermQtMatrix = data_->Get<Teuchos::RCP<Matrix> >("permQT", PermFact_.get());
  int numdf = Params().sublist(MultiGridParameterListName).get<int>("PDE equations",-1);
  int dimns = Params().sublist(MultiGridParameterListName).get<int>("null space: dimension",-1);
  if(dimns == -1 || numdf == -1) dserror("PermutedAztecSolver: Error in MueLu/ML parameters: PDE equations or null space dimension wrong.");
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = xOp->getRowMap();

  Teuchos::RCP<MultiVector> nspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
  Teuchos::RCP<std::vector<double> > nsdata = Params().sublist(MultiGridParameterListName).get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

  for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
    Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
    const size_t myLength = nspVector->getLocalLength();
    for(size_t j=0; j<myLength; j++) {
      nspVectori[j] = (*nsdata)[i*myLength+j];
    }
  }

  // calculate transformed nullspace multivector
  Teuchos::RCP<MultiVector> permutedNspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xOp->getDomainMap(),dimns,true);
  // calculate Q * b_f
  xPermQtMatrix->apply(*nspVector, *permutedNspVector, Teuchos::TRANS);

  // write data back
  for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
    Teuchos::ArrayRCP<Scalar> permutedNspVectori = permutedNspVector->getDataNonConst(i);
    const size_t myLength = permutedNspVector->getLocalLength();
      for(size_t j=0; j<myLength; j++) {
        (*nsdata)[i*myLength+j] = permutedNspVectori[j];
      }
  }
  #if 0
    // experiment
    Teuchos::RCP<MultiVector> nspVector2 = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata2 = Params().sublist(MultiGridParameterListName).get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
        Teuchos::ArrayRCP<Scalar> nspVector2i = nspVector2->getDataNonConst(i);
        const size_t myLength = nspVector2->getLocalLength();
        for(size_t j=0; j<myLength; j++) {
                nspVector2i[j] = (*nsdata2)[i*myLength+j];
        }
    }

    Teuchos::RCP<MultiVector> test = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xOp->getRowMap(),dimns,true);
    xPermQtMatrix->apply(*nspVector2, *test, Teuchos::NO_TRANS);
    test->update(-1.0, *nspVector, 1.0);
    std::cout << *test << std::endl;
  #endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::PermutedAztecSolver::FindNonDiagonalDominantRows(double diagDominanceRatio)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  // extract permuted but unscaled matrix
  Teuchos::RCP<Matrix> xPermA = data_->Get<Teuchos::RCP<Matrix> >("permA", PermFact_.get());

  const Teuchos::RCP< const Teuchos::Comm< int > > comm = xPermA->getRowMap()->getComm();

  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xPermA->getRowMap(),true);
  xPermA->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP< const Scalar > diagAVecData = diagAVec->getData(0);
  std::vector<GlobalOrdinal> NonDiagonalDominantGIDs;  // vector with nondiagonaldominant row GIDs on cur proc

  // loop over all local rows in matrix A and keep diagonal entries if corresponding
  // matrix rows are not contained in permRowMap
  for (size_t row = 0; row < xPermA->getRowMap()->getNodeNumElements(); row++) {
    GlobalOrdinal grow = xPermA->getRowMap()->getGlobalElement(row);

    // extract local row information from matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    xPermA->getLocalRowView(row, indices, vals);

    // find column entry with max absolute value
    Scalar maxVal = 0.0;
    for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
      if(std::abs(vals[j]) > maxVal) {
        maxVal = std::abs(vals[j]);
      }
    }

    // check the ratio nof the diagonal entry and the entry with
    // maximal absolute value in the current row
    // if the row is diagonal dominant the ratio is 1.0
    // if the row is not diagonal dominant, the ratio is < 1.0
    // if the row has a zero on the diagonal the ratio is zero
    // if the row is a zero row (-> singular matrix) we divide zero by zero
    if(std::abs(diagAVecData[row])/maxVal < diagDominanceRatio) { // optimal would be 1.0 -> row is diagonal dominant
      NonDiagonalDominantGIDs.push_back(grow);
    }
  }

  const Teuchos::ArrayView<const LocalOrdinal> NonDiagonalDominantGIDs_view(&NonDiagonalDominantGIDs[0],NonDiagonalDominantGIDs.size());

  Teuchos::RCP<Map> NonDiagonalDominantGIDsMap = MapFactory::Build(
      xPermA->getRowMap()->lib(),
      Teuchos::OrdinalTraits<int>::invalid(),
      NonDiagonalDominantGIDs_view,
      0, comm);

  int verbosity = Params().sublist("Aztec Parameters").get<int>("verbosity");
  if(verbosity>0)
    *fos << "PermutedAztecSolver: found " << NonDiagonalDominantGIDsMap->getGlobalNumElements() << " non-digaonal dominant entries." << std::endl;

  return NonDiagonalDominantGIDsMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::PermutedAztecSolver::FindZeroDiagonalEntries()
{
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  // extract permuted but unscaled matrix
  Teuchos::RCP<Matrix> xPermA = data_->Get<Teuchos::RCP<Matrix> >("permA", PermFact_.get());

  const Teuchos::RCP< const Teuchos::Comm< int > > comm = xPermA->getRowMap()->getComm();

  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xPermA->getRowMap(),true);
  xPermA->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP< const Scalar > diagAVecData = diagAVec->getData(0);
  LocalOrdinal lNumZeros = 0;
  std::vector<GlobalOrdinal> zeroGids;
  for(size_t i = 0; i<diagAVec->getMap()->getNodeNumElements(); ++i) {

    if(std::abs(diagAVecData[i]) < 1e-5) { // pick out all rows with very small diagonal entries
      lNumZeros++;
      zeroGids.push_back(diagAVec->getMap()->getGlobalElement(i));
    }
  }

  const Teuchos::ArrayView<const LocalOrdinal> zeroGids_view(&zeroGids[0],zeroGids.size());

  Teuchos::RCP<Map> zeroDiagonalMap = MapFactory::Build(
      xPermA->getRowMap()->lib(),
      Teuchos::OrdinalTraits<int>::invalid(),
      zeroGids_view,
      0, comm);

  int verbosity = Params().sublist("Aztec Parameters").get<int>("verbosity");
  if(verbosity>0)
    *fos << "PermutedAztecSolver: found " << zeroDiagonalMap->getGlobalNumElements() << " (near) zero diagonal entries." << std::endl;

  return zeroDiagonalMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::PermutedAztecSolver::CountZerosOnDiagonalEpetra(const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xCrsA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(Params().sublist("NodalBlockInformation").get<int>("nv")); // set nBlockSize

  return CountZerosOnDiagonal(xOp);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::PermutedAztecSolver::CountZerosOnDiagonal(const Teuchos::RCP<const Matrix>& xOp)
{
  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xOp->getRowMap(),true);
  xOp->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP< const Scalar > diagAVecData = diagAVec->getData(0);
  LocalOrdinal lNumZeros = 0;
  GlobalOrdinal gNumZeros = 0;
  for(size_t i = 0; i<diagAVec->getMap()->getNodeNumElements(); ++i) {
    if(diagAVecData[i] == 0.0) {
      lNumZeros++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  sumAll(diagAVec->getMap()->getComm(), (LocalOrdinal)lNumZeros, gNumZeros);

  return Teuchos::as<int>(gNumZeros);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::PermutedAztecSolver::Solve()
{
#ifdef WRITEOUTSTATISTICS
  Epetra_Time ttt(Comm());       // time measurement for whole routine
  ttt.ResetStartTime();
#endif

  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");

  // Allocate an aztec solver with default parameters
  // We do this every time because reusing the solver object
  // does lead to crashes that are not understood

  // create an aztec solver
  AztecOO aztec;
  aztec.SetAztecDefaults();

  // tell aztec to which stream to write
  aztec.SetOutputStream(std::cout);
  aztec.SetErrorStream(std::cerr);

  // Don't want linear problem to alter our aztec parameters (idiot feature!)
  // this is why we set our list here AFTER the linear problem has been set
  aztec.SetProblem( preconditioner_->LinearProblem() );

  {
    // We don't want to use Aztec's scaling capabilities as we prefer to do
    // the scaling ourselves (so we precisely know what happens)
    // Therefore set scaling parameter to none and reset it after aztec has made
    // its internal copy of the parameter list
    string scaling = azlist.get("scaling","none");
    azlist.set("scaling","none");
    aztec.SetParameters(azlist,false);
    azlist.set("scaling",scaling);
  }

  aztec.SetPrecOperator( preconditioner_->PrecOperator() );

  // iterate on the solution
  int iter = azlist.get("AZ_max_iter",500);
  double tol = azlist.get("AZ_tol",1.0e-6);

  // This hurts! It supresses error messages. This needs to be fixed.
#if 0
  // create an aztec convergence test as combination of
  // L2-norm and Inf-Norm to be both satisfied where we demand
  // L2 < tol and Linf < 10*tol
  {
    Epetra_Operator* op  = aztec.GetProblem()->GetOperator();
    Epetra_Vector*   rhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetRHS());
    Epetra_Vector*   lhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetLHS());
    // max iterations
    aztest_maxiter_ = Teuchos::rcp(new AztecOO_StatusTestMaxIters(iter));
    // L2 norm
    aztest_norm2_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,tol));
    aztest_norm2_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                 AztecOO_StatusTestResNorm::TwoNorm);
    aztest_norm2_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                   AztecOO_StatusTestResNorm::TwoNorm);
    // Linf norm (demanded to be 1.0 times L2-norm now, to become an input parameter?)
    aztest_norminf_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,1.0*tol));
    aztest_norminf_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                   AztecOO_StatusTestResNorm::InfNorm);
    aztest_norminf_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                     AztecOO_StatusTestResNorm::InfNorm);
    // L2 AND Linf
    aztest_combo1_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::SEQ));
    // maxiters OR (L2 AND Linf)
    aztest_combo2_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR));
    aztest_combo1_->AddStatusTest(*aztest_norm2_);
    aztest_combo1_->AddStatusTest(*aztest_norminf_);
    aztest_combo2_->AddStatusTest(*aztest_maxiter_);
    aztest_combo2_->AddStatusTest(*aztest_combo1_);
    // set status test
    aztec.SetStatusTest(aztest_combo2_.get());
  }
#endif

  // if you want to get some information on eigenvalues of the Hessenberg matrix/the
  // estimated condition number of the preconditioned system, uncomment the following
  // line and set AZOUTPUT>0 in your .dat-file
  // aztec_->SetAztecOption(AZ_solver,AZ_gmres_condnum);

  //------------------------------- just do it----------------------------------------
  aztec.Iterate(iter,tol);
  //----------------------------------------------------------------------------------

  preconditioner_->Finish( &*A_, &*x_, &*b_ );

  // check status of solution process
  const double* status = aztec.GetAztecStatus();
#if 0
  AztecOO_StatusType stat = aztest_combo2_->GetStatus();
  if (stat!=Converged)
  {
    bool resolve = false;
    if (stat==Unconverged)
    {
      if (comm_.MyPID()==0) printf("Max iterations reached in AztecOO\n");
    }
    else if (stat==Failed || stat==NaN || stat==PartialFailed)
    {
      if (comm_.MyPID()==0) printf("Numerical breakdown in AztecOO\n");
    }
    else dserror("Aztec returned unknown nonzero status %d",(int)stat);
  }
#else
  if (status[AZ_why] != AZ_normal)
  {
    if (status[AZ_why] == AZ_breakdown)
    {
      if (comm_.MyPID()==0) printf("Numerical breakdown in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (comm_.MyPID()==0) printf("Problem is near singular in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (comm_.MyPID()==0) printf("Max iterations reached in AztecOO\n");
    }
  } // if (status[AZ_why] != AZ_normal)
#endif

#ifdef WRITEOUTSTATISTICS
    if(outfile_)
    {
      fprintf(outfile_,"LinIter %i\tNumGlobalElements %i\tAZ_solve_time %f\tAztecSolveTime %f\tAztecPrecondSetup %f\t\n",
              (int)status[AZ_its],
              A_->OperatorRangeMap().NumGlobalElements(),
              status[AZ_solve_time],
              dtimeprecondsetup_ + ttt.ElapsedTime(),
              dtimeprecondsetup_);
      fflush(outfile_);
    }
#endif

  GlobalOrdinal rowperm  = 0;
  GlobalOrdinal colperm  = 0;
  GlobalOrdinal lrowperm = 0;
  GlobalOrdinal lcolperm = 0;

  if(bPermuteLinearSystem_) {
    // repermutate solution vector
    this->ReTransformSolution();
    rowperm = data_->Get<GlobalOrdinal>("#RowPermutations", PermFact_.get());
    colperm = data_->Get<GlobalOrdinal>("#ColPermutations", PermFact_.get());
    lrowperm = data_->Get<GlobalOrdinal>("#WideRangeRowPermutations", PermFact_.get());
    lcolperm = data_->Get<GlobalOrdinal>("#WideRangeColPermutations", PermFact_.get());
  }

  // print some output if desired
  if (comm_.MyPID()==0 && outfile_)
  {
    fprintf(outfile_,"AztecOO: unknowns/iterations/time/rowpermutations/colpermutations/lrowperm/lcolperm %d  %d  %f %d %d %d %d\n",
            A_->OperatorRangeMap().NumGlobalElements(),(int)status[AZ_its],status[AZ_solve_time],rowperm,colperm,lrowperm,lcolperm);
    fflush(outfile_);
  }

  ncall_ += 1;
}

void LINALG::SOLVER::PermutedAztecSolver::ReTransformSolution() {
  Teuchos::RCP<const Epetra_CrsMatrix> epPermQTMatrix = GetOperator("permQT",PermFact_);
  Teuchos::RCP<Epetra_MultiVector> xtemp = Teuchos::rcp(new Epetra_MultiVector(*x_));
  xtemp->Update(1.0,*x_,0.0);
  epPermQTMatrix->Multiply(false, *xtemp, *x_);
}

Teuchos::ParameterList & LINALG::SOLVER::PermutedAztecSolver::getPrecondParamterList()
{
  if ( Params().isSublist("IFPACK Parameters") )  return Params().sublist("IFPACK Parameters");
  else if ( Params().isSublist("ML Parameters") ) return Params().sublist("ML Parameters");
  else if ( Params().isSublist("MueLu Parameters") ) return Params().sublist("MueLu Parameters");
  else if ( Params().isSublist("MueLu (Contact) Parameters") ) return Params().sublist("MueLu (Contact) Parameters");
  else if ( Params().isSublist("MueLu (Contact2) Parameters") )return Params().sublist("MueLu (Contact2) Parameters");
  else if ( Params().isSublist("MueLu (Contact3) Parameters") )return Params().sublist("MueLu (Contact3) Parameters");
  else if ( Params().isSublist("MueLu (PenaltyContact) Parameters") )return Params().sublist("MueLu (PenaltyContact) Parameters");
  else if ( Params().isSublist("AMGBS Parameters") )return Params().sublist("AMGBS Parameters");
  else return Params();
}

Teuchos::RCP<const Epetra_CrsMatrix> LINALG::SOLVER::PermutedAztecSolver::GetOperator(const std::string name, const Teuchos::RCP<FactoryBase> & fact)
{
  Teuchos::RCP<Matrix> xPermScalOp = data_->Get<Teuchos::RCP<Matrix> >(name, fact.get());
  Teuchos::RCP<CrsMatrixWrap> xPermScalCrsOp = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xPermScalOp);
  Teuchos::RCP<CrsMatrix> xPermScalCrsMat = xPermScalCrsOp->getCrsMatrix();
  Teuchos::RCP<EpetraCrsMatrix> xEpPermScalCrsMat = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xPermScalCrsMat);
  return xEpPermScalCrsMat->getEpetra_CrsMatrix();
}

Teuchos::RCP<Epetra_CrsMatrix> LINALG::SOLVER::PermutedAztecSolver::GetOperatorNonConst(const std::string name, const Teuchos::RCP<FactoryBase> & fact)
{
  Teuchos::RCP<Matrix> xPermScalOp = data_->Get<Teuchos::RCP<Matrix> >(name, fact.get());
  Teuchos::RCP<CrsMatrixWrap> xPermScalCrsOp = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xPermScalOp);
  Teuchos::RCP<CrsMatrix> xPermScalCrsMat = xPermScalCrsOp->getCrsMatrix();
  Teuchos::RCP<EpetraCrsMatrix> xEpPermScalCrsMat = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xPermScalCrsMat);
  return xEpPermScalCrsMat->getEpetra_CrsMatrixNonConst();
}

#endif // HAVE_EXPERIMENTAL_MueLu
#endif // HAVE_MueLu





