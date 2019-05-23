/*----------------------------------------------------------------------------*/
/*!

\brief Coarse level interface of nonlinear solver to BACI

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "fas_hierarchy.H"
#include "fas_nlnlevel.H"
#include "nln_problem.H"
#include "nln_problem_coarselevel.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblemCoarseLevel::NlnProblemCoarseLevel()
    : hierarchy_(Teuchos::null), fhat_(Teuchos::null), fbar_(Teuchos::null), levelid_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // extract some parameters from parameter list
  hierarchy_ = Params()->get<Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy>>("AMG Hierarchy");
  levelid_ = Params()->get<int>("Level ID");

  // setup the underlying fine level nonlinear problem
  nlnproblem_->Setup();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemCoarseLevel::SetModelEvaluator()
{
  if (not Params()->isParameter("Field Problem"))
  {
    dserror("The parameter 'Field Problem' has not been set.");
    return false;
  }

  if (not Params()->isType<Teuchos::RCP<NLNSOL::NlnProblem>>("Field Problem"))
  {
    dserror(
        "Parameter 'Field Problem' isn't of type "
        "'Teuchos::RCP<NLNSOL::NlnProblem>'.");
    return false;
  }

  nlnproblem_ = Params()->get<Teuchos::RCP<NLNSOL::NlnProblem>>("Field Problem");

  return true;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::ComputeF(
    const Epetra_MultiVector& xc, Epetra_MultiVector& fc) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnProblemCoarseLevel::ComputeF");
  Teuchos::TimeMonitor monitor(*time);

  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // evaluate the plain residual
  Teuchos::RCP<const Epetra_MultiVector> fcoarse = ComputePlainF(xc);

  // copy result to provided output variable
  err = fc.Update(1.0, *fcoarse, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // apply coarse grid residual corrections
  if (fhat_.is_null())
  {
    dserror("Residual correction 'fhat_' not set, yet.");
  }
  if (fbar_.is_null())
  {
    dserror("Residual correction 'fbar_' not set, yet.");
  }
  err = fc.Update(-1.0, *fhat_, 1.0, *fbar_, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  return;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> NLNSOL::NlnProblemCoarseLevel::ComputeF(
    const Epetra_MultiVector& x) const
{
  Teuchos::RCP<Epetra_MultiVector> f = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  ComputeF(x, *f);

  return f;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> NLNSOL::NlnProblemCoarseLevel::ComputePlainF(
    const Epetra_MultiVector& xc) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnProblemCoarseLevel::ComputePlainF");
  Teuchos::TimeMonitor monitor(*time);

  // prolongate current solution to fine level
  Teuchos::RCP<const Epetra_MultiVector> xf = Hierarchy().ProlongateToFineLevel(xc, levelid_);

  // call evaluate on the fine level
  Teuchos::RCP<Epetra_MultiVector> ffine = Teuchos::rcp(new Epetra_MultiVector(xf->Map(), true));
  nlnproblem_->ComputeF(*xf, *ffine);

  // restrict fine level residual to current coarse level
  Teuchos::RCP<const Epetra_MultiVector> fcoarse =
      Hierarchy().RestrictToCoarseLevel(*ffine, levelid_);

  if (fcoarse.is_null())
  {
    dserror("fcoarse is a null pointer.");
  }
  if (not fcoarse.is_valid_ptr())
  {
    dserror("fcoarse is a not a valid pointer.");
  }

  return fcoarse;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::ComputeJacobian() const
{
  dserror(
      "This is algebraic multigrid. We cannot update the Jacobian on a "
      "coarse level. \nYou should re-build the RAPs with a updated fine level "
      "Jacobian.");

  return;
}

/*----------------------------------------------------------------------------*/
const NLNSOL::FAS::AMGHierarchy& NLNSOL::NlnProblemCoarseLevel::Hierarchy() const
{
  if (hierarchy_.is_null()) dserror("AMG-FAS hierarchy object 'hierarchy_' has not been set, yet.");

  return *hierarchy_;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NLNSOL::NlnProblemCoarseLevel::GetJacobianOperator() const
{
  // extract coarse level Jacobian operator from multigrid hierarchy
  Teuchos::RCP<const Epetra_Operator> jacop = Hierarchy().NlnLevel(LevelID())->GetMatrix();

  // check if Jacobian operator has already been set
  if (jacop.is_null()) dserror("Jacobian operator 'jac_' has not been initialized, yet.");

  return Teuchos::rcp_const_cast<Epetra_Operator>(jacop);
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::SetupResidualModification(
    Teuchos::RCP<const Epetra_MultiVector> xbar, Teuchos::RCP<const Epetra_MultiVector> fbar)
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnProblemCoarseLevel::SetupResidualModification");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // free the outdated vectors
  fbar_ = Teuchos::null;
  fhat_ = Teuchos::null;

  // copy 'fbar'
  fbar_ = Teuchos::rcp(new Epetra_MultiVector(*fbar));

  // compute 'fhat'
  fhat_ = ComputePlainF(*xbar);

  // safety check
  if (fhat_.is_null())
  {
    dserror("fhat_ is Teuchos::null.");
  }
  if (not fhat_.is_valid_ptr())
  {
    dserror("fhat_ is not a valid pointer.");
  }
  if (fbar_.is_null())
  {
    dserror("fbar_ is Teuchos::null.");
  }
  if (not fbar_.is_valid_ptr())
  {
    dserror("fbar_ is not a valid pointer.");
  }

  // safety check
  dsassert(not fhat_.is_null(), "fhat_ is Teuchos::null.");
  dsassert(fhat_.is_valid_ptr(), "fhat_ is not a valid pointer.");
  dsassert(not fbar_.is_null(), "fbar_ is Teuchos::null.");
  dsassert(fbar_.is_valid_ptr(), "fbar_ is not a valid pointer.");

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec,
    const std::string& description, const IO::VectorType vt) const
{
  //  if (not vec->Map().PointSameAs(DofRowMap()))
  //    dserror("Map of vector does not match map of this level.");

  // Prolongate to fine level
  Teuchos::RCP<const Epetra_MultiVector> vecfine =
      Hierarchy().ProlongateToFineLevel(*vec, LevelID());

  // Write debug output on fine level
  nlnproblem_->WriteVector(vecfine, description, vt);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::WriteVector(
    const Epetra_MultiVector& vec, const std::string& description, const IO::VectorType vt) const
{
  WriteVector(Teuchos::rcp(&vec, false), description, vt);

  return;
}
