/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_direction_defaultsteptest.hpp"

#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::check_test(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp)
{
  compute_element_volumes(dir, grp);

  gnum_bad_eles_ = 0;
  identify_bad_elements(grp, gnum_bad_eles_);
  compute_primal_direction_measures(dir, grp);

  return is_accepted();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::init_and_check_test(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp)
{
  my_bad_dofs_.clear();
  return check_test(dir, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::is_accepted()
{
  utils_->out() << '\n'
                << ::NOX::Utils::fill(40, '=') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n";

  const bool isaccepted =
      ((is_valid_direction_length() or is_valid_element_volumes()) and is_positive_definite());

  utils_->out() << ::NOX::Utils::fill(40, '=') << "\n\n";

  if (isaccepted) dirdir_last_ = dirdir_;

  return isaccepted;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::Test::VolumeChange::compute_primal_direction_measures(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp)
{
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);
  const ::NOX::Epetra::Vector& dir_epetra = dynamic_cast<const ::NOX::Epetra::Vector&>(dir);

  // range map and domain map are expected to coincide
  const Epetra_Map& rangemap = nln_grp.get_jacobian_range_map(0, 0);
  Teuchos::RCP<Epetra_Vector> primal_dir =
      Core::LinAlg::ExtractMyVector(dir_epetra.getEpetraVector(), rangemap);
  ::NOX::Epetra::Vector nox_primal_dir(primal_dir, ::NOX::Epetra::Vector::CreateView);

  Teuchos::RCP<::NOX::Epetra::Vector> result_ptr;
  ::NOX::Abstract::Group::ReturnType status =
      nln_grp.apply_jacobian_block(nox_primal_dir, result_ptr, 0, 0);

  if (status != ::NOX::Abstract::Group::Ok) FOUR_C_THROW("apply_jacobian_block failed");

  dirdir_ = nox_primal_dir.innerProduct(nox_primal_dir);
  dirres_ = nox_primal_dir.innerProduct(*result_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Direction::Test::VolumeChange::compute_element_volumes(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp)
{
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);

  const ::NOX::Abstract::Group::ReturnType eval_status =
      nln_grp.compute_element_volumes(ref_ele_vols_);

  if (eval_status != ::NOX::Abstract::Group::Ok)
    FOUR_C_THROW("The evaluation of the reference volumes failed!");

  return nln_grp.compute_trial_element_volumes(curr_ele_vols_, dir, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::Test::VolumeChange::identify_bad_elements(
    ::NOX::Abstract::Group& grp, int& gnew_num_bad_eles)
{
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);
  const Epetra_Comm& comm = ref_ele_vols_->Map().Comm();

  gnew_num_bad_eles = fill_my_bad_dofs(nln_grp);

  int gnum_bad_dofs = 0;
  int lnum_bad_dofs = my_bad_dofs_.size();
  comm.SumAll(&lnum_bad_dofs, &gnum_bad_dofs, 1);

  utils_->out() << "Number of bad dofs = " << gnum_bad_dofs << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::Nln::Direction::Test::VolumeChange::get_current_diagonal(
    const ::NOX::Abstract::Group& grp) const
{
  return get_current_diagonal(dynamic_cast<const Nln::Group&>(grp));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::Nln::Direction::Test::VolumeChange::get_current_diagonal(
    const NOX::Nln::Group& grp) const
{
  Teuchos::RCP<Epetra_Vector> diagonal = get_empty_diagonal(grp);
  fill_diagonal_at_bad_dofs(*diagonal);

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::Test::VolumeChange::fill_diagonal_at_bad_dofs(
    Epetra_Vector& diagonal) const
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
Teuchos::RCP<Epetra_Vector> NOX::Nln::Direction::Test::VolumeChange::get_empty_diagonal(
    const NOX::Nln::Group& grp) const
{
  const Epetra_Map& jac_rmap = grp.get_jacobian_range_map(0, 0);
  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::rcp(new Epetra_Vector(jac_rmap, true));

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::Nln::Direction::Test::VolumeChange::fill_my_bad_dofs(NOX::Nln::Group& grp)
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

  grp.get_dofs_from_elements(my_bad_elements, my_bad_dofs_);

  return gnew_num_bad_eles;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::is_valid_element_volumes() const
{
  const bool isvalidev = (gnum_bad_eles_ == 0);

  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << ::NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "gnum_bad_eles      = " << gnum_bad_eles_ << "\n";
  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n";

  return isvalidev;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::is_valid_direction_length() const
{
  bool isvalidlength = (dirres_ > dirdir_ and dirdir_ < dirdir_last_);

  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << ::NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "dirdir             = " << dirdir_ << "\n"
                << "dirdir_last        = " << dirdir_last_ << "\n"
                << "dirres             = " << dirres_ << "\n"
                << "IsValidLength      = " << (isvalidlength ? "TRUE" : "FALSE") << "\n";
  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n";


  return (isvalidlength);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Test::VolumeChange::is_positive_definite() const
{
  const double fac = 10 * std::numeric_limits<double>::epsilon();
  const bool isposdef = dirres_ > fac * dirdir_;

  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n"
                << "Direction::Test::VolumeChange::" << __FUNCTION__ << "\n"
                << ::NOX::Utils::fill(40, '-') << "\n";
  utils_->out() << "dirdir             = " << dirdir_ << "\n"
                << "dirres             = " << dirres_ << "\n"
                << "IsPositiveDefinite = " << (isposdef ? "TRUE" : "FALSE") << "\n";
  utils_->out() << ::NOX::Utils::fill(40, '-') << "\n";

  return (isposdef);
}

FOUR_C_NAMESPACE_CLOSE
