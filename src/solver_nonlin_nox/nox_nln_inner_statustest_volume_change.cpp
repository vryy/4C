/*-----------------------------------------------------------*/
/*! \file

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_volume_change.H"
#include "nox_nln_group.H"
#include "nox_nln_inner_statustest_interface_required.H"

#include "../drt_lib/drt_dserror.H"

#include <NOX_Utils.H>
#include <NOX_Solver_Generic.H>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::VolumeChange::CheckStatus(
    const Interface::Required& interface, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  status_ = status_unevaluated;
  const int iter_ls = interface.GetNumIterations();

  if (iter_ls == 0)
  {
    if (SetElementVolumes(grp, ref_ele_vols_) != NOX::Abstract::Group::Ok)
      dserror("The evaluation of the reference volumes failed!");
    return status_unevaluated;
  }
  else if (checkType == NOX::StatusTest::None)
  {
    curr_ele_vols_ = Teuchos::null;
    return status_unevaluated;
  }

  const NOX::Abstract::Group::ReturnType eval_status = SetElementVolumes(grp, curr_ele_vols_);

  num_bad_eles_ = NumberOfBadElements();
  if (num_bad_eles_ == 0 and eval_status == NOX::Abstract::Group::Ok)
    status_ = status_converged;
  else
    status_ = status_step_too_long;

  return status_;
}



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::VolumeChange::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::VolumeChange::Print(
    std::ostream& stream, int indent) const
{
  const std::string indent_str(indent, ' ');
  const std::string big_indent_str(13 + indent, ' ');

  stream << indent_str << status_ << "Inner Volume Change Test\n"
         << big_indent_str << NOX::Utils::sciformat(params_.lower_bound_, 3) << " < "
         << NOX::Utils::sciformat(min_vol_change_, 3) << " AND "
         << NOX::Utils::sciformat(max_vol_change_, 3) << " < "
         << NOX::Utils::sciformat(params_.upper_bound_, 3) << "\n";
  stream << big_indent_str << "Total number of bad/failing elements = " << num_bad_eles_ << " / "
         << num_failing_eles_ << " of " << ref_ele_vols_->Map().NumGlobalElements()
         << " tested elements." << std::endl;

  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::INNER::StatusTest::VolumeChange::SetElementVolumes(
    const NOX::Abstract::Group& grp, Teuchos::RCP<Epetra_Vector>& ele_vols) const
{
  const NOX::NLN::Group& nln_grp = dynamic_cast<const NOX::NLN::Group&>(grp);

  ele_vols = Teuchos::null;
  return nln_grp.computeElementVolumes(ele_vols);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::INNER::StatusTest::VolumeChange::NumberOfBadElements()
{
  min_vol_change_ = 1.0;
  max_vol_change_ = 1.0;
  num_failing_eles_ = 0;

  std::vector<int> lbadfailcount(2, 0);
  std::vector<int> gbadfailcount(2, 0);
  double lmin = 1.0;
  double lmax = 1.0;

  const unsigned num_my_eles = ref_ele_vols_->Map().NumMyElements();
  const int* my_egids = ref_ele_vols_->Map().MyGlobalElements();
  const Epetra_Comm& comm = ref_ele_vols_->Map().Comm();

  for (unsigned elid = 0; elid < num_my_eles; ++elid)
  {
    if ((*ref_ele_vols_)[elid] <= 0.0)
      dserror("Negative or zero reference volume detected for element #d!", my_egids[elid]);

    const double ratio = (*curr_ele_vols_)[elid] / (*ref_ele_vols_)[elid];

    /* Detect elements with rapidly changing volume. Note that the first
     * condition also identifies all elements for which the evaluation failed,
     * since the corresponding volume is set to a negative value. */
    if (ratio < params_.lower_bound_ or ratio > params_.upper_bound_)
    {
      ++lbadfailcount[0];
      if (ratio < 0.0) ++lbadfailcount[1];
    }

    if (ratio < lmin)
      lmin = ratio;
    else if (ratio > lmax)
      lmax = ratio;
  }

  comm.SumAll(lbadfailcount.data(), gbadfailcount.data(), 2);
  comm.MaxAll(&lmax, &max_vol_change_, 1);
  comm.MinAll(&lmin, &min_vol_change_, 1);

  num_failing_eles_ = gbadfailcount[1];

  return gbadfailcount[0];
}
