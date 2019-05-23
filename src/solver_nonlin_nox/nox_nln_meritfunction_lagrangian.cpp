/*-----------------------------------------------------------*/
/*!

\brief Implementation of the Lagrangian merit function for
       constrained problems.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_meritfunction_lagrangian.H"  // class header
#include "nox_nln_constraint_group.H"

#include "../linalg/linalg_serialdensevector.H"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Lagrangian::Lagrangian(
    const std::string& identifier, const Teuchos::RCP<NOX::Utils>& u)
    : lagrangian_type_(mrtfct_vague), meritFunctionName_()
{
  SetType(identifier);
  meritFunctionName_ = MeritFuncName2String(Type());

  utils_ = u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computef(const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computef()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
  {
    throwError("computef()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  // Get the primary contribution and constraint contributions
  return constr_grp_ptr->GetModelValue(Type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSlope(
    const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computeSlope()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
    throwError("computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");

  // compute the slope
  return constr_grp_ptr->GetLinearizedModelTerms(dir, Type(), linorder_first, lin_wrt_all_dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeMixed2ndOrderTerms(
    const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computeSlope()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
    throwError("computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");

  // compute the slope
  return constr_grp_ptr->GetLinearizedModelTerms(dir, Type(), linorder_second, lin_wrt_mixed_dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(const double& stepPV,
    const double& stepLM, const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computeSaddlePointModel()",
        "The current function value was not "
        "computed yet. Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
  {
    throwError("computeModel()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  const NOX::NLN::CONSTRAINT::Group& constr_grp = *constr_grp_ptr;

  // compute the function value
  double model = 0.0;

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the primary degrees of freedom
  // --------------------------------------------
  model += stepPV *
           constr_grp.GetLinearizedModelTerms(dir, Type(), linorder_first, lin_wrt_primary_dofs);

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier degrees of freedom
  // --------------------------------------------
  model += stepLM * constr_grp.GetLinearizedModelTerms(
                        dir, Type(), linorder_first, lin_wrt_lagrange_multiplier_dofs);

  // --------------------------------------------
  // Get the 2-nd order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier AND primary degrees of freedom
  // --------------------------------------------
  model += stepLM * stepPV *
           constr_grp.GetLinearizedModelTerms(dir, Type(), linorder_second, lin_wrt_mixed_dofs);

  return model;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(
    const double& step, const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp) const
{
  return computeSaddlePointModel(step, step, dir, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::NLN::MeritFunction::Lagrangian::name() const { return meritFunctionName_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::map<std::string, NOX::NLN::MeritFunction::MeritFctName>
NOX::NLN::MeritFunction::Lagrangian::GetSupportedTypeList() const
{
  std::map<std::string, MeritFctName> type_names;

  type_names["Lagrangian"] = mrtfct_lagrangian;
  type_names["Lagrangian Active"] = mrtfct_lagrangian_active;

  return type_names;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Lagrangian::SetType(const std::string& identifier)
{
  static const std::map<std::string, MeritFctName> supported_type_names = GetSupportedTypeList();

  auto cit = supported_type_names.cbegin();
  while (cit != supported_type_names.cend())
  {
    if (boost::iequals(identifier, cit->first))
    {
      lagrangian_type_ = cit->second;
      break;
    }
    ++cit;
  }

  if (cit == supported_type_names.cend())
  {
    std::cout << "\n\n=====================================================\n";
    std::cout << "Supported Lagrangian type names:\n"
                 "EXPECTED INPUT [= deduced merit function type]\n";
    for (const auto& supported_pair : supported_type_names)
      std::cout << supported_pair.first << " [= " << MeritFuncName2String(supported_pair.second)
                << "]\n";

    dserror("Unknown type name: \"%s\"", identifier.c_str());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Lagrangian::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_->isPrintType(NOX::Utils::Error))
  {
    utils_->out() << "ERROR - NOX::NLN::MeritFunction::Lagrangian::" << functionName << " - "
                  << errorMsg << std::endl;
  }
  throw "NOX Error";
}
