/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_cheapsimple.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_downwindmatrix.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linear_solver_preconditioner_linalg_ana.hpp"

#include <EpetraExt_OperatorOut.h>
#include <Ifpack.h>
#include <ml_MultiLevelPreconditioner.h>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

FOUR_C_NAMESPACE_OPEN

#define SIMPLEC_DIAGONAL 1       // 1: row sums     0: just diagonal
#define CHEAPSIMPLE_ALGORITHM 1  // 1: AMG          0: true solve
#define SIMPLER_ALGORITHM 0      // 1: triple solve 0: double solve
#define SIMPLER_ALPHA 0.8        // simple pressure damping parameter
#define SIMPLER_TIMING 0         // printout timing of setup

using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinearSolver::CheapSimpleBlockPreconditioner::CheapSimpleBlockPreconditioner(
    Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& predict_list,
    const Teuchos::ParameterList& correct_list)
    : predict_solver_list_(predict_list),
      schur_solver_list_(correct_list),
      alpha_(SIMPLER_ALPHA),
      vdw_(false),
      pdw_(false),
      label_(SetupLabel())
{
  // remove the SIMPLER sublist from the predictSolver_list_,
  // otherwise it will try to recursively create a SIMPLE
  // preconditioner when we do the subblock solvers
  // if (predictSolver_list_.isSublist("SIMPLER")) predictSolver_list_.remove("SIMPLER");

  // check for contact, meshtying or constraints
  // (no special functionality yet, only checking)
  //  const int myrank = A->Comm().MyPID();
  //  bool mt = schurSolver_list_.get<bool>("MESHTYING",false);
  //  if (!myrank && mt) cout << "\n**********\nMESHTYING SIMPLER\n**********\n\n";
  //  bool co = schurSolver_list_.get<bool>("CONTACT",false);
  //  if (!myrank && co) cout << "\n**********\nCONTACT SIMPLER\n**********\n\n";
  //  bool cstr = schurSolver_list_.get<bool>("CONSTRAINT",false);
  //  if (!myrank && cstr) cout << "\n**********\nCONSTRAINT SIMPLER\n**********\n\n";

  setup(A, predict_list, correct_list);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinearSolver::CheapSimpleBlockPreconditioner::setup(Teuchos::RCP<Epetra_Operator> A,
    const Teuchos::ParameterList& origvlist, const Teuchos::ParameterList& origplist)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  const int myrank = A->Comm().MyPID();
  Teuchos::Time time("", true);
  Teuchos::Time totaltime("", true);
  const bool visml = predict_solver_list_.isSublist("ML Parameters");
  const bool pisml = schur_solver_list_.isSublist("ML Parameters");
  const bool vismuelu = predict_solver_list_.isSublist("MueLu Parameters");
  const bool pismuelu = schur_solver_list_.isSublist("MueLu Parameters");
  const bool visifpack = predict_solver_list_.isSublist("IFPACK Parameters");
  const bool pisifpack = schur_solver_list_.isSublist("IFPACK Parameters");

  if (!visml && !visifpack && !vismuelu)
    FOUR_C_THROW("Have to use either ML or Ifpack for velocities");
  if (!pisml && !pisifpack && !pismuelu)
    FOUR_C_THROW("Have to use either ML or Ifpack for pressure");

  //-------------------------------------------------------------------------
  // either do manual split or use provided BlockSparseMatrixBase
  //-------------------------------------------------------------------------
  a_ = Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(A);
  if (a_ != Teuchos::null)
  {
    // Make a shallow copy of the block matrix as the preconditioners on the
    // blocks will be reused and the next assembly will replace the block
    // matrices.
    a_ = a_->Clone(Core::LinAlg::View);
    mmex_ = a_->RangeExtractor();
  }

  // if the MESHTYING, CONTACT or CONSTRAINT flag is set,
  // adapt multigrid nullspace for constraint block

  //-------------------------------------------------------------------------
  // Modify lists to reuse subblock preconditioner at least maxiter times
  //-------------------------------------------------------------------------
  {
    int maxiter = predict_solver_list_.sublist("Belos Parameters").get("Maximum Iterations", 1);
    predict_solver_list_.sublist("Belos Parameters").set("reuse", maxiter + 1);
  }

#if SIMPLEC_DIAGONAL
  //-------------------------------------------------------------------------
  // Allocate and compute abs(rowsum(A(0,0))^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0), false);
    Teuchos::RCP<Epetra_CrsMatrix> A00 = (*a_)(0, 0).EpetraMatrix();
    A00->InvRowSums(diag);
    diag_ainv_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(diag));
    diag_ainv_->Complete(*mmex_.Map(0), *mmex_.Map(0));
  }
#else
  //-------------------------------------------------------------------------
  // Allocate and compute diag(A(0,0)^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0), false);
    (*A_)(0, 0).ExtractDiagonalCopy(diag);
    int err = diag.Reciprocal(diag);
    if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal returned %d", err);
    diagAinv_ = Teuchos::rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0), *mmex_.Map(0));
  }
#endif
  if (!myrank && SIMPLER_TIMING)
    printf("--- Time to do diagF^{-1}   %10.3E\n", time.totalElapsedTime(true));
  time.reset();

  //-------------------------------------------------------------------------
  // Allocate and compute approximate Schur complement operator S
  // S = A(1,1) - A(1,0) * diagAinv * A(0,1)
  //-------------------------------------------------------------------------
  {
    Teuchos::Time ltime("", true);
    // with Trilinos Q1/2013 there are some improvements in EpetraExt MM.
    // However, they lead to a crash here -> use MLMultiply instead.
    // S_ = Core::LinAlg::Multiply(*diagAinv_,false,(*A_)(0,1),false,true);
    s_ = Core::LinAlg::MLMultiply(*diag_ainv_, (*a_)(0, 1), true);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = diagAinv * A(0,1) %10.3E\n", ltime.totalElapsedTime(true));
    ltime.reset();
    s_ = Core::LinAlg::MLMultiply((*a_)(1, 0), *s_, false);
    // The Core::LinAlg::Multiply method would consume a HUGE amount of memory!!!
    // So always use Core::LinAlg::MLMultiply in the line above! Otherwise you won't be able
    // to solve any large linear problem since you'll definitely run out of memory.
    // S_ = Core::LinAlg::Multiply((*A_)(1,0),false,*S_,false,false);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = A(1,0) * S (ML)   %10.3E\n", ltime.totalElapsedTime(true));
    ltime.reset();
    s_->Add((*a_)(1, 1), false, 1.0, -1.0);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = A(1,1) - S        %10.3E\n", ltime.totalElapsedTime(true));
    ltime.reset();
    s_->Complete((*a_)(1, 1).DomainMap(), (*a_)(1, 1).RangeMap());
    if (!myrank && SIMPLER_TIMING)
      printf("*** S complete            %10.3E\n", ltime.totalElapsedTime(true));
    ltime.reset();
  }
  if (!myrank && SIMPLER_TIMING)
    printf("--- Time to do S            %10.3E\n", time.totalElapsedTime(true));
  time.reset();

#if CHEAPSIMPLE_ALGORITHM
  {
    Epetra_CrsMatrix* A00 = nullptr;
    Epetra_CrsMatrix* A11 = nullptr;
    A00 = (*a_)(0, 0).EpetraMatrix().get();
    A11 = s_->EpetraMatrix().get();

    //-------------------------------------------------------------------------
    // Allocate preconditioner for pressure and velocity
    //-------------------------------------------------------------------------
    if (visml)
    {
      predict_solver_list_.sublist("ML Parameters").remove("init smoother", false);
      ppredict_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(
          *A00, predict_solver_list_.sublist("ML Parameters"), true));
    }
    else if (vismuelu)
    {
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix = Teuchos::rcp(new Epetra_CrsMatrix(*A00));

      // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
          Teuchos::rcp(new EpetraCrsMatrix(Pmatrix));
      Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

      // Create MueLu preconditioner:  using ML like input or xml files
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = Teuchos::null;
      Teuchos::ParameterList& MueLuList = predict_solver_list_.sublist("MueLu Parameters");
      std::string xmlFileName = MueLuList.get<std::string>("xml file", "none");
      if (xmlFileName == "none")
      {
        MueLu::MLParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(MueLuList);
        H = mueLuFactory.CreateHierarchy();
        H->GetLevel(0)->Set("A", mueluOp);
        mueLuFactory.SetupHierarchy(*H);
      }
      else
      {
        // Exctract info from list
        int numdf = MueLuList.get<int>("PDE equations", -1);
        int dimns = MueLuList.get<int>("null space: dimension", -1);

        Teuchos::RCP<Epetra_MultiVector> nsdata =
            MueLuList.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
        if (numdf < 1 or dimns < 1)
          FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");
        if (nsdata == Teuchos::null) FOUR_C_THROW("Error: null space data is empty");

        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector =
            Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(nsdata));

        Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();
        nspVector->replaceMap(rowMap);

        // Create Hierarchy
        mueluOp->SetFixedBlockSize(numdf);
        MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
            xmlFileName, *(mueluOp->getRowMap()->getComm()));
        H = mueLuFactory.CreateHierarchy();
        H->SetDefaultVerbLevel(MueLu::Extreme);
        H->GetLevel(0)->Set("A", mueluOp);
        H->GetLevel(0)->Set("Nullspace", nspVector);
        H->GetLevel(0)->setlib(Xpetra::UseEpetra);
        H->setlib(Xpetra::UseEpetra);
        mueLuFactory.SetupHierarchy(*H);
      }

      // set preconditioner
      ppredict_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
    else
    {
      std::string type =
          predict_solver_list_.sublist("Belos Parameters").get("Preconditioner Type", "ILU");
      Ifpack factory;
      Ifpack_Preconditioner* prec = factory.Create(type, A00, 0);
      prec->SetParameters(predict_solver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      ppredict_ = Teuchos::rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING)
      printf("--- Time to do P(v)         %10.3E\n", time.totalElapsedTime(true));
    time.reset();

    if (pisml)
    {
      schur_solver_list_.sublist("ML Parameters").remove("init smoother", false);
      pschur_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(
          *A11, schur_solver_list_.sublist("ML Parameters"), true));
    }
    else if (pismuelu)
    {
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix = Teuchos::rcp(new Epetra_CrsMatrix(*A11));

      // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
          Teuchos::rcp(new EpetraCrsMatrix(Pmatrix));
      Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

      // Create MueLu preconditioner:  using ML like input or xml files
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = Teuchos::null;
      Teuchos::ParameterList& MueLuList = schur_solver_list_.sublist("MueLu Parameters");
      std::string xmlFileName = MueLuList.get<std::string>("xml file", "none");
      if (xmlFileName == "none")
      {
        MueLu::MLParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(MueLuList);
        H = mueLuFactory.CreateHierarchy();
        H->GetLevel(0)->Set("A", mueluOp);
        mueLuFactory.SetupHierarchy(*H);
      }
      else
      {
        // Exctract info from list
        int numdf = MueLuList.get<int>("PDE equations", -1);
        int dimns = MueLuList.get<int>("null space: dimension", -1);

        Teuchos::RCP<Epetra_MultiVector> nsdata =
            MueLuList.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
        if (numdf < 1 or dimns < 1)
          FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");
        if (nsdata == Teuchos::null) FOUR_C_THROW("Error: null space data is empty");

        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector =
            Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(nsdata));

        Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();
        nspVector->replaceMap(rowMap);

        // Create Hierarchy
        mueluOp->SetFixedBlockSize(numdf);
        MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
            xmlFileName, *(mueluOp->getRowMap()->getComm()));
        H = mueLuFactory.CreateHierarchy();
        H->SetDefaultVerbLevel(MueLu::Extreme);
        H->GetLevel(0)->Set("A", mueluOp);
        H->GetLevel(0)->Set("Nullspace", nspVector);
        H->GetLevel(0)->setlib(Xpetra::UseEpetra);
        H->setlib(Xpetra::UseEpetra);
        mueLuFactory.SetupHierarchy(*H);
      }

      // set preconditioner
      pschur_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
    else
    {
      Ifpack factory;
      std::string type =
          schur_solver_list_.sublist("Belos Parameters").get("Preconditioner Type", "ILU");
      Ifpack_Preconditioner* prec = factory.Create(type, A11, 0);
      prec->SetParameters(schur_solver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      pschur_ = Teuchos::rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING)
      printf("--- Time to do P(p)         %10.3E\n", time.totalElapsedTime(true));
    time.reset();
  }
#else
  //-------------------------------------------------------------------------
  // Allocate solver for pressure and velocity
  //-------------------------------------------------------------------------
  {
    Teuchos::RCP<Teuchos::ParameterList> vrcplist = Teuchos::rcp(&predictSolver_list_, false);
    vsolver_ = Teuchos::rcp(new Core::LinAlg::Solver(vrcplist, A_->Comm(), outfile_));
    Teuchos::RCP<Teuchos::ParameterList> prcplist = Teuchos::rcp(&schurSolver_list_, false);
    psolver_ = Teuchos::rcp(new Core::LinAlg::Solver(prcplist, A_->Comm(), outfile_));
  }
#endif

  //-------------------------------------------------------------------------
  // Allocate velocity and pressure solution and rhs vectors
  //-------------------------------------------------------------------------
  vx_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(0), false));
  vb_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(0), false));
  px_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(1), false));
  pb_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(1), false));

  //-------------------------------------------------------------------------
  // Allocate working vectors for velocity and pressure
  //-------------------------------------------------------------------------
  vwork1_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(0), false));
  vwork2_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(0), false));
  pwork1_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(1), false));
  pwork2_ = Teuchos::rcp(new Core::LinAlg::Ana::Vector(*mmex_.Map(1), false));

  if (!myrank && SIMPLER_TIMING)
    printf("--- Time to do allocate mem %10.3E\n", time.totalElapsedTime(true));
  if (!myrank && SIMPLER_TIMING)
    printf("=== Total simpler setup === %10.3E\n", totaltime.totalElapsedTime(true));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinearSolver::CheapSimpleBlockPreconditioner::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // note: might pass X and Y as physically identical objects,
  // so we better deep copy here

  // extract initial guess and rhs for velocity and pressure
  mmex_.ExtractVector(X, 0, *vb_);
  mmex_.ExtractVector(X, 1, *pb_);

#if CHEAPSIMPLE_ALGORITHM  // SIMPLE and SIMPLEC but without solve, just AMG
  cheap_simple(*vx_, *px_, *vb_, *pb_);
#else

#if SIMPLER_ALGORITHM
  Simpler(*vx_, *px_, *vb_, *pb_);
#else
  Simple(*vx_, *px_, *vb_, *pb_);
#endif

#endif

  // insert solution for velocity and pressure
  mmex_.InsertVector(*vx_, 0, Y);
  mmex_.InsertVector(*px_, 1, Y);

  return 0;
}


/*----------------------------------------------------------------------*
 | Pernice, M., Tocci, M.D.:                                            |
 | A Multigrid Preconditioned Newton-Krylov method for the incomp.      |
 | Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)  |
 *----------------------------------------------------------------------*/
void Core::LinearSolver::CheapSimpleBlockPreconditioner::simpler(Core::LinAlg::Ana::Vector& vx,
    Core::LinAlg::Ana::Vector& px, Core::LinAlg::Ana::Vector& vb,
    Core::LinAlg::Ana::Vector& pb) const
{
  using namespace Core::LinAlg::Ana;
  Core::LinAlg::SparseMatrix& A00 = (*a_)(0, 0);
  Core::LinAlg::SparseMatrix& A10 = (*a_)(1, 0);
  Core::LinAlg::SparseMatrix& A01 = (*a_)(0, 1);
  Core::LinAlg::SparseMatrix& diagAinv = *diag_ainv_;
  Core::LinAlg::SparseMatrix& S = *s_;

  //-------------------------------------------------- L-solve / U-solve

  px = inverse(S, *psolver_, false) * (pb - A10 * (diagAinv * vb));

  vx = inverse(A00, *vsolver_, false) * (vb - A01 * px);

  //------------------------------------------------ Implicit projection

  if (alpha_ != 1.0) px *= alpha_;

  *vwork2_ = diagAinv * A01 * (inverse(S, *psolver_, false) * (A10 * vx));

  vx -= vwork2_;
}


/*----------------------------------------------------------------------*
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 *----------------------------------------------------------------------*/
void Core::LinearSolver::CheapSimpleBlockPreconditioner::simple(Core::LinAlg::Ana::Vector& vx,
    Core::LinAlg::Ana::Vector& px, Core::LinAlg::Ana::Vector& vb,
    Core::LinAlg::Ana::Vector& pb) const
{
  using namespace Core::LinAlg::Ana;
  Core::LinAlg::SparseMatrix& A00 = (*a_)(0, 0);
  Core::LinAlg::SparseMatrix& A10 = (*a_)(1, 0);
  Core::LinAlg::SparseMatrix& A01 = (*a_)(0, 1);
  Core::LinAlg::SparseMatrix& diagAinv = *diag_ainv_;
  Core::LinAlg::SparseMatrix& S = *s_;


  //------------------------------------------------------------ L-solve

  *vwork1_ = inverse(A00, *vsolver_, false) * vb;

  px = inverse(S, *psolver_, false) * (pb - A10 * vwork1_);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);
}


/*----------------------------------------------------------------------*
 | is a cheaper variation from:                                         |
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 |                                                                      |
 | all solves replaced by single AMG sweeps                             |
 *----------------------------------------------------------------------*/
void Core::LinearSolver::CheapSimpleBlockPreconditioner::cheap_simple(Core::LinAlg::Ana::Vector& vx,
    Core::LinAlg::Ana::Vector& px, Core::LinAlg::Ana::Vector& vb,
    Core::LinAlg::Ana::Vector& pb) const
{
  Core::LinAlg::SparseMatrix& A10 = (*a_)(1, 0);
  Core::LinAlg::SparseMatrix& A01 = (*a_)(0, 1);
  Core::LinAlg::SparseMatrix& diagAinv = *diag_ainv_;

  //------------------------------------------------------------ L-solve
  if (vdw_)
  {
    vdwind_->Permute(&vb, &*vdwin_);
    ppredict_->ApplyInverse(*vdwin_, *vdwout_);
    vdwind_->InvPermute(&*vdwout_, &*vwork1_);
  }
  else
    ppredict_->ApplyInverse(vb, *vwork1_);

  *pwork1_ = pb - A10 * vwork1_;

  if (pdw_)
  {
    pdwind_->Permute(&*pwork1_, &*pdwin_);
    pschur_->ApplyInverse(*pdwin_, *pdwout_);
    pdwind_->InvPermute(&*pdwout_, &px);
  }
  else
    pschur_->ApplyInverse(*pwork1_, px);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);
}

FOUR_C_NAMESPACE_CLOSE
