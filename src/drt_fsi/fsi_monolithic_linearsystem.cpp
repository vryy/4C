


#include "fsi_monolithic_linearsystem.H"
#include "../drt_lib/drt_dserror.H"

// NOX includes
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Scaling.H"
#include "NOX_Utils.H"

// External include files for Epetra
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "AztecOO_StatusTest.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"


#include "AztecOO.h"

#include <typeinfo>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicLinearSystem::MonolithicLinearSystem(Teuchos::ParameterList& printingParams,
    Teuchos::ParameterList& linearSolverParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<Epetra_Operator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<Epetra_Operator>& M, const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::Epetra::LinearSystemAztecOO(
          printingParams, linearSolverParams, iJac, J, iPrec, M, cloneVector, scalingObject),
      lsparams_(linearSolverParams)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicLinearSystem::~MonolithicLinearSystem() {}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicLinearSystem::applyJacobianInverse(
    Teuchos::ParameterList& p, const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result)
{
  // AGS: Rare option, similar to Max Iters=1 but twice as fast.
  if (p.get("Use Preconditioner as Solver", false))
    return applyRightPreconditioning(false, p, input, result);

  // Aztec crashes in reuses because its buggy
  Teuchos::RCP<Epetra_Operator> prec = solvePrecOpPtr;
  aztecSolverPtr = Teuchos::null;
  aztecSolverPtr = Teuchos::rcp(new AztecOO());
  bool isprec = isPrecConstructed;
  softreset(lsparams_);
  isPrecConstructed = isprec;
  solvePrecOpPtr = prec;

  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess) result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(
      jacPtr.get(), &(result.getEpetraVector()), &(nonConstInput.getEpetraVector()));

  // Set objects in aztec solver.
  // RPP: Don't use "aztecSolverPtr->SetProblem(Problem);", it breaks
  // things if you rebuild the prec.  Also, you MUST set Jac operator
  // before Prec Operator in AztecOO object.
  this->setAztecOOJacobian();
  if (solvePrecOpPtr == Teuchos::null) dserror("Preconditioner is zero");
  aztecSolverPtr->SetPrecOperator(solvePrecOpPtr.get());
  // this->setAztecOOPreconditioner();
  aztecSolverPtr->SetLHS(&(result.getEpetraVector()));
  aztecSolverPtr->SetRHS(&(nonConstInput.getEpetraVector()));


  // ************* Begin linear system scaling *******************
  if (!Teuchos::is_null(scaling))
  {
    if (!manualScaling) scaling->computeScaling(Problem);

    scaling->scaleLinearSystem(Problem);

    if (utils.isPrintType(NOX::Utils::Details))
    {
      utils.out() << *scaling << std::endl;
    }
  }
  // ************* End linear system scaling *******************


  // Make sure preconditioner was constructed if requested
  if (!isPrecConstructed && (precAlgorithm != None_))
  {
    throwError("applyJacobianInverse",
        "Preconditioner is not constructed!  Call createPreconditioner() first.");
  }

  // Get linear solver convergence parameters
  int maxit = p.get("Max Iterations", 400);
  double tol = p.get("Tolerance", 1.0e-6);


  int aztecStatus = -1;

  // RPP: This is a hack to get aztec to reuse the preconditioner.
  // There is a bug in AztecOO in how it stores old
  // preconditioners. The storage bin is based on the aztec options
  // list.  When we call ConstructPreconditioner, the option for
  // az_pre_calc is set to AZ_calc, so the preconditioner is stored
  // with that option.  But then the routine toggles the AZ_pre_calc
  // flag to AZ_reuse.  When we call iterate, the first solve works,
  // but subsequent solves fail to find the preconditioner.
  // Will try get Alan to fix for release 7.0.
  if (precAlgorithm == AztecOO_ && precReusePolicy == PRPT_REUSE)
    aztecSolverPtr->SetAztecOption(AZ_pre_calc, AZ_calc);

  // build the status test
  {
    Epetra_Operator* op = jacPtr.get();
    Epetra_Vector* rhs = &(nonConstInput.getEpetraVector());
    Epetra_Vector* lhs = &(result.getEpetraVector());
    if (!op || !rhs || !lhs) dserror("One of the objects in linear system is NULL");

    // max iterations
    aztest_maxiter_ = Teuchos::rcp(new AztecOO_StatusTestMaxIters(maxit));
    // L2 norm
    aztest_norm2_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op, *lhs, *rhs, tol));
    aztest_norm2_->DefineResForm(
        AztecOO_StatusTestResNorm::Implicit, AztecOO_StatusTestResNorm::TwoNorm);
    aztest_norm2_->DefineScaleForm(
        AztecOO_StatusTestResNorm::NormOfInitRes, AztecOO_StatusTestResNorm::TwoNorm);
    // Linf norm (demanded to be 1 times L2-norm now, to become an input parameter)
    aztest_norminf_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op, *lhs, *rhs, 1.0 * tol));
    aztest_norminf_->DefineResForm(
        AztecOO_StatusTestResNorm::Implicit, AztecOO_StatusTestResNorm::InfNorm);
    aztest_norminf_->DefineScaleForm(
        AztecOO_StatusTestResNorm::NormOfInitRes, AztecOO_StatusTestResNorm::InfNorm);
    // L2 AND Linf
    aztest_combo1_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::AND));
    // maxiters OR (L2 AND Linf)
    aztest_combo2_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR));
    aztest_combo1_->AddStatusTest(*aztest_norm2_);
    aztest_combo1_->AddStatusTest(*aztest_norminf_);
    aztest_combo2_->AddStatusTest(*aztest_maxiter_);
    aztest_combo2_->AddStatusTest(*aztest_combo1_);
    // set status test
    aztecSolverPtr->SetStatusTest(aztest_combo2_.get());
  }

  aztecStatus = aztecSolverPtr->Iterate(maxit, tol);

  // Unscale the linear system
  if (!Teuchos::is_null(scaling)) scaling->unscaleLinearSystem(Problem);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails)
  {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    int prevLinIters = outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = 0;
    double achievedTol = -1.0;
    curLinIters = aztecSolverPtr->NumIters();
    achievedTol = aztecSolverPtr->ScaledResidual();

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  if (aztecStatus != 0) return false;

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLinearSystem::softreset(Teuchos::ParameterList& linearSolverParams)
{
  // === soft means do not touch the preconditioner!!
  // destroyPreconditioner();

  // Set the requested preconditioning.
  std::string prec = linearSolverParams.get("Preconditioner", "None");
  if (prec == "AztecOO")
    precAlgorithm = AztecOO_;
  else if (prec == "Ifpack")
    precAlgorithm = Ifpack_;
  else if (prec == "New Ifpack")
    precAlgorithm = NewIfpack_;
#ifdef HAVE_NOX_ML_EPETRA
  else if (prec == "ML")
    precAlgorithm = ML_;
#endif
  else if (prec == "User Defined")
    precAlgorithm = UserDefined_;
  else if (prec == "None")
    precAlgorithm = None_;
  else
  {
    std::string errorMessage = "Option for \"Preconditioner\" is invalid!";
    throwError("reset()", errorMessage);
  }

  // Make sure the correct objects were supplied for the requested
  // preconditioning choice.
  checkPreconditionerValidity();

  zeroInitialGuess = linearSolverParams.get("Zero Initial Guess", false);

  manualScaling = linearSolverParams.get("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails = linearSolverParams.get("Output Solver Details", true);

  throwErrorOnPrecFailure = linearSolverParams.get("Throw Error on Prec Failure", true);

  // The first time a SetProblem is used on the AztecOO solver
  // it sets all aztec options based on the Epetra_LinearProblem
  // options. Subsequent calls do not.  We call this here so we
  // can set our own parameters and not have them overwritten
  // by the first SetProblem call.
  // RPP: Not any more.  We don't set the solver objects with
  // the problem class.  It cause seg faults when new preconditioners
  // were computed during prec reuse.
  // Epetra_LinearProblem& problem = *(new Epetra_LinearProblem);
  // aztecSolverPtr->SetProblem(problem);

  // Set the Jacobian in the solver. It must be set before
  // a preconditioner can be set.
  //   if ((jacType == EpetraRowMatrix) ||
  //       (jacType == EpetraVbrMatrix) ||
  //       (jacType == EpetraCrsMatrix)) {
  //     aztecSolverPtr->SetUserMatrix(dynamic_cast<Epetra_RowMatrix*>(jacPtr.get()));
  //   }
  //   else
  //     aztecSolverPtr->SetUserOperator(jacPtr.get());

  // Set the major aztec options.  Must be called after the first
  // SetProblem() call.
  setAztecOptions(linearSolverParams, *aztecSolverPtr);

  // Setup the preconditioner reuse policy
  std::string preReusePolicyName = linearSolverParams.get("Preconditioner Reuse Policy", "Rebuild");
  if (preReusePolicyName == "Rebuild")
    precReusePolicy = PRPT_REBUILD;
  else if (preReusePolicyName == "Recompute")
    precReusePolicy = PRPT_RECOMPUTE;
  else if (preReusePolicyName == "Reuse")
    precReusePolicy = PRPT_REUSE;
  else
  {
    std::string errorMessage =
        "Option for \"Preconditioner Reuse Policy\" is invalid! \nPossible options are \"Reuse\", "
        "\"Rebuild\", and \"Recompute\".";
    throwError("reset()", errorMessage);
  }
  maxAgeOfPrec = linearSolverParams.get("Max Age Of Prec", 1);
  precQueryCounter = 0;

#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  linearSolveCount = 0;
#endif
#endif
}
