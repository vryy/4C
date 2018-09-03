/*!----------------------------------------------------------------------
\file solver_cheapsimplepreconditioner.cpp

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

// Trilinos headers
#include <EpetraExt_OperatorOut.h>
#include <Ifpack.h>
#include <ml_MultiLevelPreconditioner.h>
#ifdef HAVE_MueLu
#include <MueLu.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_EpetraOperator.hpp>   // Aztec interface

// some typedefs
typedef Scalar SC;
typedef LocalOrdinal LO;
typedef GlobalOrdinal GO;
typedef Node NO;
#endif

// BACI headers
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"  // helper functions (linear Algebra related)
#include "../linalg/linalg_downwindmatrix.H"
#include "../linalg/linalg_multiply.H"

#include "solver_cheapsimplepreconditioner.H"

#define SIMPLEC_DIAGONAL 1       // 1: row sums     0: just diagonal
#define CHEAPSIMPLE_ALGORITHM 1  // 1: AMG          0: true solve
#define SIMPLER_ALGORITHM 0      // 1: triple solve 0: double solve
#define SIMPLER_ALPHA 0.8        // simple pressure damping parameter
#define SIMPLER_TIMING 0         // printout timing of setup
/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/08|
 *----------------------------------------------------------------------*/
LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::CheapSIMPLE_BlockPreconditioner(
    Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& predict_list,
    const Teuchos::ParameterList& correct_list, FILE* outfile)
    : outfile_(outfile),
      predictSolver_list_(predict_list),
      schurSolver_list_(correct_list),
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

  Setup(A, predict_list, correct_list);

  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Setup(Teuchos::RCP<Epetra_Operator> A,
    const Teuchos::ParameterList& origvlist, const Teuchos::ParameterList& origplist)
{
  const int myrank = A->Comm().MyPID();
  Epetra_Time time(A->Comm());
  Epetra_Time totaltime(A->Comm());
  const bool visml = predictSolver_list_.isSublist("ML Parameters");
  const bool pisml = schurSolver_list_.isSublist("ML Parameters");
  const bool vismuelu = predictSolver_list_.isSublist("MueLu Parameters");
  const bool pismuelu = schurSolver_list_.isSublist("MueLu Parameters");
  const bool visifpack = predictSolver_list_.isSublist("IFPACK Parameters");
  const bool pisifpack = schurSolver_list_.isSublist("IFPACK Parameters");

  if (!visml && !visifpack && !vismuelu) dserror("Have to use either ML or Ifpack for velocities");
  if (!pisml && !pisifpack && !pismuelu) dserror("Have to use either ML or Ifpack for pressure");

  //-------------------------------------------------------------------------
  // either do manual split or use provided BlockSparseMatrixBase
  //-------------------------------------------------------------------------
  A_ = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if (A_ != Teuchos::null)
  {
    // Make a shallow copy of the block matrix as the preconditioners on the
    // blocks will be reused and the next assembly will replace the block
    // matrices.
    A_ = A_->Clone(View);
    mmex_ = A_->RangeExtractor();
  }

  // if the MESHTYING, CONTACT or CONSTRAINT flag is set,
  // adapt multigrid nullspace for constraint block

  //-------------------------------------------------------------------------
  // Modify lists to reuse subblock preconditioner at least maxiter times
  //-------------------------------------------------------------------------
  {
    int maxiter = predictSolver_list_.sublist("Aztec Parameters").get("AZ_max_iter", 1);
    predictSolver_list_.sublist("Aztec Parameters").set("reuse", maxiter + 1);
  }

#if SIMPLEC_DIAGONAL
  //-------------------------------------------------------------------------
  // Allocate and compute abs(rowsum(A(0,0))^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0), false);
    Teuchos::RCP<Epetra_CrsMatrix> A00 = (*A_)(0, 0).EpetraMatrix();
    A00->InvRowSums(diag);
    diagAinv_ = Teuchos::rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0), *mmex_.Map(0));
  }
#else
  //-------------------------------------------------------------------------
  // Allocate and compute diag(A(0,0)^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0), false);
    (*A_)(0, 0).ExtractDiagonalCopy(diag);
    int err = diag.Reciprocal(diag);
    if (err) dserror("Epetra_MultiVector::Reciprocal returned %d", err);
    diagAinv_ = Teuchos::rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0), *mmex_.Map(0));
  }
#endif
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do diagF^{-1}   %10.3E\n", time.ElapsedTime());
  time.ResetStartTime();

  //-------------------------------------------------------------------------
  // Allocate and compute approximate Schur complement operator S
  // S = A(1,1) - A(1,0) * diagAinv * A(0,1)
  //-------------------------------------------------------------------------
  {
    Epetra_Time ltime(A_->Comm());
    // with Trilinos Q1/2013 there are some improvements in EpetraExt MM.
    // However, they lead to a crash here -> use MLMultiply instead.
    // S_ = LINALG::Multiply(*diagAinv_,false,(*A_)(0,1),false,true);
    S_ = LINALG::MLMultiply(*diagAinv_, (*A_)(0, 1), true);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = diagAinv * A(0,1) %10.3E\n", ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_ = LINALG::MLMultiply((*A_)(1, 0), *S_, false);
    // The LINALG::Multiply method would consume a HUGE amount of memory!!!
    // So always use LINALG::MLMultiply in the line above! Otherwise you won't be able
    // to solve any large linear problem since you'll definitely run out of memory.
    // S_ = LINALG::Multiply((*A_)(1,0),false,*S_,false,false);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = A(1,0) * S (ML)   %10.3E\n", ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_->Add((*A_)(1, 1), false, 1.0, -1.0);
    if (!myrank && SIMPLER_TIMING)
      printf("*** S = A(1,1) - S        %10.3E\n", ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_->Complete((*A_)(1, 1).DomainMap(), (*A_)(1, 1).RangeMap());
    if (!myrank && SIMPLER_TIMING)
      printf("*** S complete            %10.3E\n", ltime.ElapsedTime());
    ltime.ResetStartTime();
  }
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do S            %10.3E\n", time.ElapsedTime());
  time.ResetStartTime();

#if CHEAPSIMPLE_ALGORITHM
  {
    Epetra_CrsMatrix* A00 = NULL;
    Epetra_CrsMatrix* A11 = NULL;
    A00 = (*A_)(0, 0).EpetraMatrix().get();
    A11 = S_->EpetraMatrix().get();

    //-------------------------------------------------------------------------
    // Allocate preconditioner for pressure and velocity
    //-------------------------------------------------------------------------
    if (visml)
    {
      predictSolver_list_.sublist("ML Parameters").remove("init smoother", false);
      Ppredict_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(
          *A00, predictSolver_list_.sublist("ML Parameters"), true));
    }
    else if (vismuelu)
    {
#ifdef HAVE_MueLu
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix = Teuchos::rcp(new Epetra_CrsMatrix(*A00));

      // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix));
      Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

      // Create MueLu preconditioner:  using ML like input or xml files
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = Teuchos::null;
      Teuchos::ParameterList& MueLuList = predictSolver_list_.sublist("MueLu Parameters");
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
        Teuchos::RCP<std::vector<double>> nsdata =
            MueLuList.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
        if (numdf < 1 or dimns < 1) dserror("Error: PDE equations or null space dimension wrong.");
        if (nsdata == Teuchos::null) dserror("Error: null space data is empty");
        // Convert null space to muelu format
        Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();
        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector =
            Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
                rowMap, dimns, true);
        for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
        {
          Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
          const size_t myLength = nspVector->getLocalLength();
          for (size_t j = 0; j < myLength; j++)
          {
            nspVectori[j] = (*nsdata)[i * myLength + j];
          }
        }
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
      Ppredict_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
#else
      dserror("BACI has been compiled without MueLu support.");
#endif
    }
    else
    {
      std::string type =
          predictSolver_list_.sublist("Aztec Parameters").get("Preconditioner Type", "ILU");
      Ifpack factory;
      Ifpack_Preconditioner* prec = factory.Create(type, A00, 0);
      prec->SetParameters(predictSolver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      Ppredict_ = Teuchos::rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING)
      printf("--- Time to do P(v)         %10.3E\n", time.ElapsedTime());
    time.ResetStartTime();

    if (pisml)
    {
      schurSolver_list_.sublist("ML Parameters").remove("init smoother", false);
      Pschur_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(
          *A11, schurSolver_list_.sublist("ML Parameters"), true));
    }
    else if (pismuelu)
    {
#ifdef HAVE_MueLu
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix = Teuchos::rcp(new Epetra_CrsMatrix(*A11));

      // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix));
      Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

      // Create MueLu preconditioner:  using ML like input or xml files
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = Teuchos::null;
      Teuchos::ParameterList& MueLuList = schurSolver_list_.sublist("MueLu Parameters");
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
        Teuchos::RCP<std::vector<double>> nsdata =
            MueLuList.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
        if (numdf < 1 or dimns < 1) dserror("Error: PDE equations or null space dimension wrong.");
        if (nsdata == Teuchos::null) dserror("Error: null space data is empty");
        // Convert null space to muelu format
        Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();
        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector =
            Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
                rowMap, dimns, true);
        for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
        {
          Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
          const size_t myLength = nspVector->getLocalLength();
          for (size_t j = 0; j < myLength; j++)
          {
            nspVectori[j] = (*nsdata)[i * myLength + j];
          }
        }
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
      Pschur_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
#else
      dserror("BACI has been compiled without MueLu support.");
#endif
    }
    else
    {
      Ifpack factory;
      std::string type =
          schurSolver_list_.sublist("Aztec Parameters").get("Preconditioner Type", "ILU");
      Ifpack_Preconditioner* prec = factory.Create(type, A11, 0);
      prec->SetParameters(schurSolver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      Pschur_ = Teuchos::rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING)
      printf("--- Time to do P(p)         %10.3E\n", time.ElapsedTime());
    time.ResetStartTime();
  }
#else
  //-------------------------------------------------------------------------
  // Allocate solver for pressure and velocity
  //-------------------------------------------------------------------------
  {
    Teuchos::RCP<Teuchos::ParameterList> vrcplist = Teuchos::rcp(&predictSolver_list_, false);
    vsolver_ = Teuchos::rcp(new LINALG::Solver(vrcplist, A_->Comm(), outfile_));
    Teuchos::RCP<Teuchos::ParameterList> prcplist = Teuchos::rcp(&schurSolver_list_, false);
    psolver_ = Teuchos::rcp(new LINALG::Solver(prcplist, A_->Comm(), outfile_));
  }
#endif

  //-------------------------------------------------------------------------
  // Allocate velocity and pressure solution and rhs vectors
  //-------------------------------------------------------------------------
  vx_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0), false));
  vb_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0), false));
  px_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1), false));
  pb_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1), false));

  //-------------------------------------------------------------------------
  // Allocate working vectors for velocity and pressure
  //-------------------------------------------------------------------------
  vwork1_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0), false));
  vwork2_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0), false));
  pwork1_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1), false));
  pwork2_ = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1), false));

  if (!myrank && SIMPLER_TIMING) printf("--- Time to do allocate mem %10.3E\n", time.ElapsedTime());
  if (!myrank && SIMPLER_TIMING)
    printf("=== Total simpler setup === %10.3E\n", totaltime.ElapsedTime());

  return;
}


/*----------------------------------------------------------------------*
 |  apply const operator (public)                            mwgee 02/08|
 *----------------------------------------------------------------------*/
int LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // note: Aztec might pass X and Y as physically identical objects,
  // so we better deep copy here

  // extract initial guess and rhs for velocity and pressure
  mmex_.ExtractVector(X, 0, *vb_);
  mmex_.ExtractVector(X, 1, *pb_);

#if CHEAPSIMPLE_ALGORITHM  // SIMPLE and SIMPLEC but without solve, just AMG
  CheapSimple(*vx_, *px_, *vb_, *pb_);
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
 |  (private)                                                mwgee 02/08|
 | taken from:                                                          |
 | Pernice, M., Tocci, M.D.:                                            |
 | A Multigrid Preconditioned Newton-Krylov method for the incomp.      |
 | Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)  |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Simpler(LINALG::ANA::Vector& vx,
    LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  using namespace LINALG::ANA;
  SparseMatrix& A00 = (*A_)(0, 0);
  SparseMatrix& A10 = (*A_)(1, 0);
  SparseMatrix& A01 = (*A_)(0, 1);
  SparseMatrix& diagAinv = *diagAinv_;
  SparseMatrix& S = *S_;

  //-------------------------------------------------- L-solve / U-solve

  px = inverse(S, *psolver_, false) * (pb - A10 * (diagAinv * vb));

  vx = inverse(A00, *vsolver_, false) * (vb - A01 * px);

  //------------------------------------------------ Implicit projection

  if (alpha_ != 1.0) px *= alpha_;

  *vwork2_ = diagAinv * A01 * (inverse(S, *psolver_, false) * (A10 * vx));

  vx -= vwork2_;

  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 | taken from:                                                          |
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Simple(LINALG::ANA::Vector& vx,
    LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  using namespace LINALG::ANA;
  LINALG::SparseMatrix& A00 = (*A_)(0, 0);
  SparseMatrix& A10 = (*A_)(1, 0);
  SparseMatrix& A01 = (*A_)(0, 1);
  SparseMatrix& diagAinv = *diagAinv_;
  SparseMatrix& S = *S_;


  //------------------------------------------------------------ L-solve

  *vwork1_ = inverse(A00, *vsolver_, false) * vb;

  px = inverse(S, *psolver_, false) * (pb - A10 * vwork1_);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);

  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 | is a cheaper variation from:                                         |
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 |                                                                      |
 | all solves replaced by single AMG sweeps                             |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::CheapSimple(LINALG::ANA::Vector& vx,
    LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  SparseMatrix& A10 = (*A_)(1, 0);
  SparseMatrix& A01 = (*A_)(0, 1);
  SparseMatrix& diagAinv = *diagAinv_;

  //------------------------------------------------------------ L-solve
  if (vdw_)
  {
    vdwind_->Permute(&vb, &*vdwin_);
    Ppredict_->ApplyInverse(*vdwin_, *vdwout_);
    vdwind_->InvPermute(&*vdwout_, &*vwork1_);
  }
  else
    Ppredict_->ApplyInverse(vb, *vwork1_);

  *pwork1_ = pb - A10 * vwork1_;

  if (pdw_)
  {
    pdwind_->Permute(&*pwork1_, &*pdwin_);
    Pschur_->ApplyInverse(*pdwin_, *pdwout_);
    pdwind_->InvPermute(&*pdwout_, &px);
  }
  else
    Pschur_->ApplyInverse(*pwork1_, px);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);

  return;
}
