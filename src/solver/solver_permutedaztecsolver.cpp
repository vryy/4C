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

  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xCrsA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);

  if(data_==Teuchos::null) data_ = Teuchos::rcp(new Level());
  data_->setDefaultVerbLevel(Teuchos::VERB_NONE);
  data_->Set("A",xOp);

  //Teuchos::RCP<Epetra_Map> epMasterDofMap = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap");
  Teuchos::RCP<Epetra_Map> epSlaveDofMap  = Params().sublist("Aztec Parameters").get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
  //Teuchos::RCP<Epetra_Map> epActiveDofMap = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap");

  // build map extractor from different maps
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));

  data_->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xSlaveDofMap));  // set map with active dofs

  // define permutation factory for permuting the full matrix A
  PermFact_ = Teuchos::rcp(new PermutationFactory("SlaveDofMap", MueLu::NoFactory::getRCP()));

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

  // build permutation operators
  PermFact_->Build(*data_);

  // extract permuted matrix from Level
  Teuchos::RCP<Epetra_CrsMatrix> xEpPermCrsMat = GetOperatorNonConst("A",PermFact_);

  Teuchos::RCP<const Epetra_CrsMatrix> epPermPMatrix  = GetOperator("permP",  PermFact_);             // row permutation matrix
  Teuchos::RCP<const Epetra_CrsMatrix> epPermScalingMatrix = GetOperator("permScaling",PermFact_); // leftScaling matrix
  ////

  Teuchos::RCP<Epetra_MultiVector> btemp1 = Teuchos::rcp(new Epetra_MultiVector(*b));

  // P_trafo*b
  epPermPMatrix->Multiply(false, *b, *btemp1);
  // P_scaling * P_trafo * b
  epPermScalingMatrix->Multiply(false, *btemp1, *b);

  ////

  // retransform nullspace vectors
  // TODO: make this working for ML, too
  if(Params().isSublist("MueLu Parameters")) {
    Teuchos::RCP<Matrix> xPermQtMatrix = data_->Get<Teuchos::RCP<Matrix> >("permQT", PermFact_.get());
    int numdf = Params().sublist("MueLu Parameters").get<int>("PDE equations",-1);
    int dimns = Params().sublist("MueLu Parameters").get<int>("null space: dimension",-1);
    if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = xOp->getRowMap();

    Teuchos::RCP<MultiVector> nspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata = Params().sublist("MueLu Parameters").get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

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
    Teuchos::RCP<std::vector<double> > nsdata2 = Params().sublist("MueLu Parameters").get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

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

  ////

  x_ = x;
  b_ = b;
  A_ = this->GetOperatorNonConst("A",PermFact_);//xEpPermCrsMat->getEpetra_CrsMatrixNonConst();

  //data_->print(*fos,Teuchos::VERB_EXTREME);

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

  // repermutate solution vector

    /*Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    data_->print(*fos,Teuchos::VERB_EXTREME);*/

  // extract permuted matrix from Level
  Teuchos::RCP<Matrix> xPermQT = data_->Get<Teuchos::RCP<Matrix> >("permQT", PermFact_.get());
  Teuchos::RCP<CrsMatrixWrap> xPermCrsQT = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xPermQT);
  Teuchos::RCP<CrsMatrix> xPermCrsQTMat = xPermCrsQT->getCrsMatrix();
  Teuchos::RCP<EpetraCrsMatrix> xEpPermCrsQTMat = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xPermCrsQTMat);

  Teuchos::RCP<const Epetra_CrsMatrix> epPermQTMatrix = xEpPermCrsQTMat->getEpetra_CrsMatrix();

  Teuchos::RCP<Epetra_MultiVector> xtemp = Teuchos::rcp(new Epetra_MultiVector(*x_));
  //Teuchos::RCP<Epetra_MultiVector> xtemp2 = Teuchos::rcp(new Epetra_MultiVector(*x_));
  xtemp->Update(1.0,*x_,0.0);
  epPermQTMatrix->Multiply(false, *xtemp, *x_);

  //x_->Update(1.0,*xtemp2,0.0);

  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  //data_->print(*fos,Teuchos::VERB_EXTREME);

  // print some output if desired
  if (comm_.MyPID()==0 && outfile_)
  {
    fprintf(outfile_,"AztecOO: unknowns/iterations/time %d  %d  %f\n",
            A_->OperatorRangeMap().NumGlobalElements(),(int)status[AZ_its],status[AZ_solve_time]);
    fflush(outfile_);
  }

  ncall_ += 1;
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





