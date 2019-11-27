/*-----------------------------------------------------------*/
/*! \file

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_direction_defaultsteptest.H"
#include "nox_nln_group.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <NOX_Utils.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::checkTest(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp)
{
  computeElementVolumes(dir, grp);

  gnum_bad_eles_ = 0;
  identifyBadElements(grp, gnum_bad_eles_);
  computePrimalDirectionMeasures(dir, grp);

  return isAccepted();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::initAndCheckTest(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp)
{
  my_bad_dofs_.clear();
  return checkTest(dir, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::isAccepted()
{
  utils_->out() << '\n'
                << NOX::Utils::fill(40, '=') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n";

  const bool isaccepted =
      ((isValidDirectionLength() or isValidElementVolumes()) and isPositiveDefinite());

  utils_->out() << NOX::Utils::fill(40, '=') << "\n\n";

  if (isaccepted) dirdir_last_ = dirdir_;

  return isaccepted;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::Test::VolumeChange::computePrimalDirectionMeasures(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp)
{
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);
  const NOX::Epetra::Vector& dir_epetra = dynamic_cast<const NOX::Epetra::Vector&>(dir);

  // range map and domain map are expected to coincide
  const Epetra_Map& rangemap = nln_grp.getJacobianRangeMap(0, 0);
  Teuchos::RCP<Epetra_Vector> primal_dir =
      LINALG::ExtractMyVector(dir_epetra.getEpetraVector(), rangemap);
  NOX::Epetra::Vector nox_primal_dir(primal_dir, NOX::Epetra::Vector::CreateView);

  Teuchos::RCP<NOX::Epetra::Vector> result_ptr;
  NOX::Abstract::Group::ReturnType status =
      nln_grp.applyJacobianBlock(nox_primal_dir, result_ptr, 0, 0);

  if (status != NOX::Abstract::Group::Ok) dserror("applyJacobianBlock failed");

  dirdir_ = nox_primal_dir.innerProduct(nox_primal_dir);
  dirres_ = nox_primal_dir.innerProduct(*result_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Direction::Test::VolumeChange::computeElementVolumes(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp)
{
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);

  const NOX::Abstract::Group::ReturnType eval_status = nln_grp.computeElementVolumes(ref_ele_vols_);

  if (eval_status != NOX::Abstract::Group::Ok)
    dserror("The evaluation of the reference volumes failed!");

  return nln_grp.computeTrialElementVolumes(curr_ele_vols_, dir, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::Test::VolumeChange::identifyBadElements(
    NOX::Abstract::Group& grp, int& gnew_num_bad_eles)
{
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);
  const Epetra_Comm& comm = ref_ele_vols_->Map().Comm();

  gnew_num_bad_eles = fillMyBadDofs(nln_grp);

  int gnum_bad_dofs = 0;
  int lnum_bad_dofs = my_bad_dofs_.size();
  comm.SumAll(&lnum_bad_dofs, &gnum_bad_dofs, 1);

  utils_->out() << "Number of bad dofs = " << gnum_bad_dofs << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::Direction::Test::VolumeChange::getCurrentDiagonal(
    const NOX::Abstract::Group& grp) const
{
  return getCurrentDiagonal(dynamic_cast<const NLN::Group&>(grp));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::Direction::Test::VolumeChange::getCurrentDiagonal(
    const NOX::NLN::Group& grp) const
{
  Teuchos::RCP<Epetra_Vector> diagonal = getEmptyDiagonal(grp);
  fillDiagonalAtBadDofs(*diagonal);

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::Test::VolumeChange::fillDiagonalAtBadDofs(Epetra_Vector& diagonal) const
{
  for (int i : my_bad_dofs_)
  {
    const int dof_lid = diagonal.Map().LID(i);
    if (dof_lid == -1) continue;

    diagonal[dof_lid] = 1.0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::Direction::Test::VolumeChange::getEmptyDiagonal(
    const NOX::NLN::Group& grp) const
{
  const Epetra_Map& jac_rmap = grp.getJacobianRangeMap(0, 0);
  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::rcp(new Epetra_Vector(jac_rmap, true));

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::Direction::Test::VolumeChange::fillMyBadDofs(NOX::NLN::Group& grp)
{
  const unsigned num_my_eles = ref_ele_vols_->Map().NumMyElements();
  const int* my_egids = ref_ele_vols_->Map().MyGlobalElements();
  const Epetra_Comm& comm = ref_ele_vols_->Map().Comm();

  std::vector<int> my_bad_elements;
  my_bad_elements.reserve(ref_ele_vols_->Map().NumMyElements());

  int lnew_num_bad_eles = 0;

  for (unsigned elid = 0; elid < num_my_eles; ++elid)
  {
    /* Detect elements with rapidly changing volume. Note that the first
     * condition also identifies all elements for which the evaluation failed,
     * since the corresponding volume is set to a negative value. */
    if ((*curr_ele_vols_)[elid] / (*ref_ele_vols_)[elid] < 0.5 or
        (*curr_ele_vols_)[elid] / (*ref_ele_vols_)[elid] > 2.0)
    {
      my_bad_elements.push_back(my_egids[elid]);
      ++lnew_num_bad_eles;
    }
  }

  int gnew_num_bad_eles = 0;
  comm.SumAll(&lnew_num_bad_eles, &gnew_num_bad_eles, 1);

  grp.getDofsFromElements(my_bad_elements, my_bad_dofs_);

  return gnew_num_bad_eles;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::isValidElementVolumes() const
{
  const bool isvalidev = (gnum_bad_eles_ == 0);

  utils_->out() << NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "gnum_bad_eles      = " << gnum_bad_eles_ << "\n";
  utils_->out() << NOX::Utils::fill(40, '-') << "\n";

  return isvalidev;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::isValidDirectionLength() const
{
  bool isvalidlength = (dirres_ > dirdir_ and dirdir_ < dirdir_last_);

  utils_->out() << NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "dirdir             = " << dirdir_ << "\n"
                << "dirdir_last        = " << dirdir_last_ << "\n"
                << "dirres             = " << dirres_ << "\n"
                << "IsValidLength      = " << (isvalidlength ? "TRUE" : "FALSE") << "\n";
  utils_->out() << NOX::Utils::fill(40, '-') << "\n";


  return (isvalidlength);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Test::VolumeChange::isPositiveDefinite() const
{
  const double fac = 10 * std::numeric_limits<double>::epsilon();
  const bool isposdef = dirres_ > fac * dirdir_;

  utils_->out() << NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "dirdir             = " << dirdir_ << "\n"
                << "dirres             = " << dirres_ << "\n"
                << "IsPositiveDefinite = " << (isposdef ? "TRUE" : "FALSE") << "\n";
  utils_->out() << NOX::Utils::fill(40, '-') << "\n";

  return (isposdef);
}
