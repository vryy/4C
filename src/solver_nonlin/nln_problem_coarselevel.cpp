/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_coarselevel.cpp

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
#include "nln_operator_fas.H"
#include "nln_problem_coarselevel.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblemCoarseLevel::NlnProblemCoarseLevel()
: hierarchy_(Teuchos::null),
  fhat_(Teuchos::null),
  fbar_(Teuchos::null),
  levelid_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::Setup()
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  hierarchy_ = Params().get<Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> >(
      "AMG Hierarchy");
  levelid_ = Params().get<int>("Level ID");

  // call base class
  NLNSOL::NlnProblem::Setup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::ComputeF(const Epetra_MultiVector& xc,
    Epetra_MultiVector& fc) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnProblemCoarseLevel::ComputeF");
  Teuchos::TimeMonitor monitor(*time);

  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // prolongate current solution to fine level
  Teuchos::RCP<Epetra_MultiVector> xf =
      Hierarchy().ProlongateToFineLevel(xc, levelid_);

  // call evaluate on the fine level
  Teuchos::RCP<Epetra_MultiVector> ffine =
      Teuchos::rcp(new Epetra_MultiVector(xf->Map(), true));
  NLNSOL::NlnProblem::ComputeF(*xf,*ffine);

  // restrict fine level residual to current coarse level
  Teuchos::RCP<Epetra_MultiVector> fcoarse =
      Hierarchy().RestrictToCoarseLevel(*ffine, levelid_);

  // residual correction on coarse level
  if (fhat_.is_null()) { dserror("Residual correction 'fhat_' not set, yet."); }
  if (fbar_.is_null()) { dserror("Residual correction 'fbar_' not set, yet."); }
  err = fcoarse->Update(-1.0, *fhat_, 1.0, *fbar_, 1.0);
  if (err != 0) { dserror("Update failed."); }

  // copy result to provided output variable
  err = fc.Update(1.0, *fcoarse, 0.0);
  if (err != 0) { dserror("Update failed."); }

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::ComputeJacobian() const
{
  dserror("This is algebraic multigrid. We cannot update the Jacobian on a "
      "coarse level. \nYou should re-build the RAPs with a updated fine level "
      "Jacobian.");

  return;
}

/*----------------------------------------------------------------------------*/
const NLNSOL::FAS::AMGHierarchy&
NLNSOL::NlnProblemCoarseLevel::Hierarchy() const
{
  if (hierarchy_.is_null())
    dserror("AMG-FAS hierarchy object 'hierarchy_' has not been set, yet.");

  return *hierarchy_;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::SetFHatFBar(
    Teuchos::RCP<const Epetra_MultiVector> fhat,
    Teuchos::RCP<const Epetra_MultiVector> fbar)
{
//  fhat_ = fhat; //ToDo (mayr) switch back to pointer assignment instead of allocating new vectors?
//  fbar_ = fbar;

  fhat_ = Teuchos::rcp(new Epetra_MultiVector(*fhat));
  fbar_ = Teuchos::rcp(new Epetra_MultiVector(*fbar));

  dsassert(not fhat_.is_null(), "fhat_ is Teuchos::null.");
  dsassert(fhat_.is_valid_ptr(), "fhat_ is not a valid pointer.");
  dsassert(not fbar_.is_null(), "fbar_ is Teuchos::null.");
  dsassert(fbar_.is_valid_ptr(), "fbar_ is not a valid pointer.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator>
NLNSOL::NlnProblemCoarseLevel::GetJacobianOperator() const
{
  // extract coarse level Jacobian operator from multigrid hierarchy
  Teuchos::RCP<const Epetra_Operator> jacop =
      Hierarchy().NlnLevel(LevelID())->GetMatrix();

  // check if Jacobian operator has already been set
  if (jacop.is_null())
    dserror("Jacobian operator 'jac_' has not been initialized, yet.");

  return Teuchos::rcp_const_cast<Epetra_Operator>(jacop);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::WriteVector(
    Teuchos::RCP<const Epetra_MultiVector> vec, const std::string& description,
    const IO::DiscretizationWriter::VectorType vt) const
{
  if (not vec->Map().PointSameAs(DofRowMap()))
    dserror("Map of vector does not match map of this level.");

  // Prolongate to fine level
  Teuchos::RCP<Epetra_MultiVector> vecfine =
      Hierarchy().ProlongateToFineLevel(*vec, LevelID());

  // Write debug output on fine level
  NLNSOL::NlnProblem::WriteVector(vecfine, description, vt);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemCoarseLevel::WriteVector(
    const Epetra_MultiVector& vec, const std::string& description,
    const IO::DiscretizationWriter::VectorType vt) const
{
  WriteVector(Teuchos::rcp(&vec, false), description, vt);

  return;
}
