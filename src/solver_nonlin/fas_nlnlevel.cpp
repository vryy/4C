/*----------------------------------------------------------------------------*/
/*!
\file fas_nlnlevel.cpp

\brief Nonlinear level containing the representation and operations for a single nonlinear multigrid level

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <iostream>

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "fas_nlnlevel.H"
#include "nln_operator_factory.H"
#include "nln_operator.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::FAS::NlnLevel::NlnLevel()
: isinit_(false),
  issetup_(false),
  A_(Teuchos::null),
  rop_(Teuchos::null),
  pop_(Teuchos::null),
  xbar_(Teuchos::null),
  fbar_(Teuchos::null),
  fhat_(Teuchos::null),
  nlnproblem_(Teuchos::null),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  presmoother_(Teuchos::null),
  postsmoother_(Teuchos::null),
  coarsesolver_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Initialization */
void NLNSOL::FAS::NlnLevel::Init(const int levelid,
                                 const int numlevels,
                                 Teuchos::RCP<const Epetra_CrsMatrix> A,
                                 Teuchos::RCP<const Epetra_CrsMatrix> R,
                                 Teuchos::RCP<const Epetra_CrsMatrix> P,
                                 const Epetra_Comm& comm,
                                 const Teuchos::ParameterList& params,
                                 Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem
                                 )
{
  // store input arguments into member variables
  SetLevelID(levelid);
  SetNumLevels(numlevels);
  SetMatrix(A);
  SetROp(R);
  SetPOp(P);

  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  nlnproblem_ = nlnproblem;

  // Init() has been called
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup */
void NLNSOL::FAS::NlnLevel::Setup()
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // Factory for creating all level smoothers
  NLNSOL::NlnOperatorFactory operatorfactory;

  // create coarse level solver if this level is the coarsest level
  if (LevelID() == NumLevels() - 1)
  {
    coarsesolver_ = operatorfactory.Create(Params().sublist("Coarse Level Solver"));
    coarsesolver_->Init(Comm(), Params().sublist("Coarse Level Solver"), NlnProblem());
    coarsesolver_->Setup();
  }
  else // otherwise create pre- and post-smoother
  {
    // create presmoother
    presmoother_ = operatorfactory.Create(Params().sublist("Presmoother"));
    presmoother_->Init(Comm(), Params().sublist("Presmoother"), NlnProblem());
    presmoother_->Setup();

    // create postsmoother
    postsmoother_ = operatorfactory.Create(Params().sublist("Postsmoother"));
    postsmoother_->Init(Comm(), Params().sublist("Postsmoother"), NlnProblem());
    postsmoother_->Setup();
  }

  // Setup() has been called
  issetup_ = true;

  return;
}
/*----------------------------------------------------------------------------*/
/* Set matrix */
void NLNSOL::FAS::NlnLevel::SetMatrix(Teuchos::RCP<const Epetra_CrsMatrix> myA)
{
  if (not myA.is_null())
    A_ = Teuchos::rcp(new Epetra_CrsMatrix(*myA));

  return;
}

/*----------------------------------------------------------------------------*/
/* Set restriction operator */
void NLNSOL::FAS::NlnLevel::SetROp(Teuchos::RCP<const Epetra_CrsMatrix> myR)
{
  if (not myR.is_null())
    rop_ = Teuchos::rcp(new Epetra_CrsMatrix(*myR));

  return;
}

/*----------------------------------------------------------------------------*/
/* Set prolongation operator */
void NLNSOL::FAS::NlnLevel::SetPOp(Teuchos::RCP<const Epetra_CrsMatrix> myP)
{
  if (not myP.is_null())
    pop_ = Teuchos::rcp(new Epetra_CrsMatrix(*myP));

  return;
}

/*----------------------------------------------------------------------------*/
/* Restriction vector to next coarser level */
int NLNSOL::FAS::NlnLevel::RestrictToNextCoarserLevel(Teuchos::RCP<Epetra_MultiVector>& vec) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  int err = 0;

  if (HaveROp() and LevelID() > 0)
  {
    Teuchos::RCP<Epetra_MultiVector> tempvec = Teuchos::rcp(new Epetra_MultiVector(rop_->RowMap(), 1, true));

#ifdef DEBUG
    if (not rop_->DomainMap().PointSameAs(vec->Map()))
      dserror("Maps do not match.");
    if (not rop_->RowMap().PointSameAs(tempvec->Map()))
      dserror("Maps do not match.");
#endif

    err = rop_->Apply(*vec, *tempvec);
    if (err != 0) { dserror("Failed."); }
    vec = Teuchos::rcp(new Epetra_MultiVector(*tempvec));
  }
  else // no restriction on finest level
  {
    dserror("You are on the fine level. There is no restriction operator available.");
    err = 1;
  }

  return err;
}

/*----------------------------------------------------------------------------*/
/* Prolongate vector to next finer level */
int NLNSOL::FAS::NlnLevel::ProlongateToNextFinerLevel(Teuchos::RCP<Epetra_MultiVector>& vec) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  int err = 0;

  if (HavePOp() and (LevelID() > 0))
  {
    Teuchos::RCP<Epetra_MultiVector> tempvec = Teuchos::rcp(new Epetra_MultiVector(pop_->RowMap(), 1, true));

#ifdef DEBUG
    if (not pop_->DomainMap().PointSameAs(vec->Map()))
      dserror("Maps do not match.");
    if (not pop_->RowMap().PointSameAs(tempvec->Map()))
      dserror("Maps do not match.");
#endif

    err = pop_->Apply(*vec, *tempvec);
    if (err != 0) { dserror("Failed."); }
    vec = Teuchos::rcp(new Epetra_MultiVector(*tempvec));
  }
  else // no prolongation on finest level
  {
    dserror("You are on the fine level, already. No further prolongation possible.");
    err = 1;
  }

  return err;
}

/*----------------------------------------------------------------------------*/
/* Print method */
void NLNSOL::FAS::NlnLevel::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    // print level ID
    os << "Level " << LevelID() << " out of " << NumLevels() << ":\n";

    // print status of matrix and transfer operators
    os << "\tMatrix A set: ";
    if (HaveMatrix())
      os << "\t\tyes";
    else
      os << "\t\tno";

    os << "\n";

    os << "\tOperator R set: ";
    if (HaveROp())
      os << "\tyes";
    else
      os << "\tno";

    os << "\n";

    os << "\tOperator P set: ";
    if (HavePOp())
      os << "\tyes";
    else
      os << "\tno";

    os << "\n";

    if (HaveMatrix())
      os << "\tSize:\t\t\t" << A_->NumGlobalRows();

    os << "\n";

    os << "\tNlnProblem Type:\t" << NlnProblem()->Label();
    os << "\n";

#if 0
    // -------------------------------------------------------------------------
    // print A, P, and R to Matlab format
    // -------------------------------------------------------------------------
    if (HaveMatrix())
    {
      os << "A\n";
      A_->Print(os);
      std::ostringstream fname;
      fname << "matrix_" << LevelID() << ".mtl";
      LINALG::PrintMatrixInMatlabFormat(fname.str().c_str(),*A_);
    }

    if (HavePOp())
    {
      os << "POp_\n";
      pop_->Print(os);
      std::ostringstream fname;
      fname << "prolongator_" << LevelID() << ".mtl";
      LINALG::PrintMatrixInMatlabFormat(fname.str().c_str(),*pop_);
    }

    if (HaveROp())
    {
      os << "ROp_\n";
      rop_->Print(os);
      std::ostringstream fname;
      fname << "restrictor_" << LevelID() << ".mtl";
      LINALG::PrintMatrixInMatlabFormat(fname.str().c_str(),*rop_);
    }
    // -------------------------------------------------------------------------
#endif

    // finish it
    os << std::flush;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* Return communicator */
const Epetra_Comm& NLNSOL::FAS::NlnLevel::Comm() const
{
  if (comm_.is_null())
    dserror("Communicator has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/* Return parameter list */
const Teuchos::ParameterList& NLNSOL::FAS::NlnLevel::Params() const
{
  if (params_.is_null())
    dserror("Parameter list has not been set, yet.");

  return *params_;
}

/*----------------------------------------------------------------------------*/
/* Return dof row map of this level */
const Epetra_BlockMap& NLNSOL::FAS::NlnLevel::DofRowMap() const
{
  if (A_.is_null())
    dserror("Matrix 'A_' has not been set.");

  return A_->RowMap();
}

/*----------------------------------------------------------------------------*/
/* Access to nonlinear problem */
Teuchos::RCP<NLNSOL::NlnProblem> NLNSOL::FAS::NlnLevel::NlnProblem() const
{
  if (nlnproblem_.is_null())
    dserror("Nonlinear problem 'nlnproblem_' has not been set, yet.");

  return nlnproblem_;

}

/*----------------------------------------------------------------------------*/
/* Do presmoothing sweeps on this level */
int NLNSOL::FAS::NlnLevel::DoPreSmoothing(const Epetra_MultiVector& fhatbar,
                                          Epetra_MultiVector& x
                                          ) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  if (Comm().MyPID() == 0)
  {
    IO::cout << IO::endl
             << "Do presmoothing on level " << LevelID() << ", now."
             << IO::endl;
  }

  return presmoother_->ApplyInverse(fhatbar, x);
}

/*----------------------------------------------------------------------------*/
/* Do postsmoothing sweeps on this level */
int NLNSOL::FAS::NlnLevel::DoPostSmoothing(const Epetra_MultiVector& fhatbar,
                                           Epetra_MultiVector& x
                                           ) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  if (Comm().MyPID() == 0)
  {
    IO::cout << IO::endl
             << "Do postsmoothing on level " << LevelID() << ", now."
             << IO::endl;
  }

  return postsmoother_->ApplyInverse(fhatbar,x);
}

/*----------------------------------------------------------------------------*/
/* Do coarse level solve on this level */
int NLNSOL::FAS::NlnLevel::DoCoarseLevelSolve(const Epetra_MultiVector& fhatbar,
                                              Epetra_MultiVector& x
                                              ) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  if (coarsesolver_.is_null())
    dserror("There is no coarse level solver on level %d. Perhaps this is not the coarsest level.", LevelID());

  if (Comm().MyPID() == 0)
    IO::cout << "*** Do coarse level solve on level " << LevelID() << ", now." << IO::endl;

  return coarsesolver_->ApplyInverse(fhatbar,x);
}

/*----------------------------------------------------------------------------*/
/* Is this level the coarsest level? */
const bool NLNSOL::FAS::NlnLevel::IsCoarsestLevel() const
{
  if (LevelID() == NumLevels() - 1)
  {
    return true;  // coarsest level
  }
  else if (LevelID() < NumLevels() - 1)
  {
    return false; // fine/medium level
  }
  else // catch non-valid level IDs
  {
    dserror("LevelID %i is not a valid level.", LevelID());
    return false;
  }
}

/*----------------------------------------------------------------------------*/
/* Compute full approximation of coarse level residual */
Teuchos::RCP<Epetra_MultiVector> NLNSOL::FAS::NlnLevel::EvaluateCoarseLevelResidual
(
  const Epetra_MultiVector& x,
  const Epetra_MultiVector& fhatbar
) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  Teuchos::RCP<Epetra_MultiVector> f = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  int err = f->Update(1.0, fhatbar, 1.0);
  if (err != 0) { dserror("Failed!"); }

  return f;
}

/*----------------------------------------------------------------------------*/
/* Check for convergence */
bool NLNSOL::FAS::NlnLevel::ConvergenceCheck(const Epetra_MultiVector& f) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  double fnorm2 = 1.0e+12;

  return ConvergenceCheck(f, fnorm2);
}

/*----------------------------------------------------------------------------*/
/* Check for convergence */
bool NLNSOL::FAS::NlnLevel::ConvergenceCheck(const Epetra_MultiVector& f,
    double& fnorm2
    ) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  return NlnProblem()->ConvergenceCheck(f, fnorm2);
}
