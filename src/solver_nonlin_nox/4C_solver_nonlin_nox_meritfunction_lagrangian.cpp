/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of the Lagrangian merit function for
       constrained problems.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_meritfunction_lagrangian.hpp"  // class header

#include "4C_linalg_serialdensevector.hpp"
#include "4C_solver_nonlin_nox_constraint_group.hpp"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Lagrangian::Lagrangian(
    const std::string& identifier, const Teuchos::RCP<::NOX::Utils>& u)
    : lagrangian_type_(mrtfct_vague), merit_function_name_()
{
  set_type(identifier);
  merit_function_name_ = MeritFuncName2String(Type());

  utils_ = u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computef(const ::NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throw_error("computef()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
  {
    throw_error("computef()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  // Get the primary contribution and constraint contributions
  return constr_grp_ptr->GetModelValue(Type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSlope(
    const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throw_error("computeSlope()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
    throw_error("computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");

  // compute the slope
  return constr_grp_ptr->get_linearized_model_terms(dir, Type(), linorder_first, lin_wrt_all_dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::compute_mixed2nd_order_terms(
    const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throw_error("computeSlope()",
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
    throw_error("computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");

  // compute the slope
  return constr_grp_ptr->get_linearized_model_terms(
      dir, Type(), linorder_second, lin_wrt_mixed_dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::compute_saddle_point_model(const double& stepPV,
    const double& stepLM, const ::NOX::Abstract::Vector& dir,
    const ::NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throw_error("compute_saddle_point_model()",
        "The current function value was not "
        "computed yet. Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr)
  {
    throw_error("computeModel()", "Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  const NOX::NLN::CONSTRAINT::Group& constr_grp = *constr_grp_ptr;

  // compute the function value
  double model = 0.0;

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the primary degrees of freedom
  // --------------------------------------------
  model += stepPV *
           constr_grp.get_linearized_model_terms(dir, Type(), linorder_first, lin_wrt_primary_dofs);

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier degrees of freedom
  // --------------------------------------------
  model += stepLM * constr_grp.get_linearized_model_terms(
                        dir, Type(), linorder_first, lin_wrt_lagrange_multiplier_dofs);

  // --------------------------------------------
  // Get the 2-nd order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier AND primary degrees of freedom
  // --------------------------------------------
  model += stepLM * stepPV *
           constr_grp.get_linearized_model_terms(dir, Type(), linorder_second, lin_wrt_mixed_dofs);

  return model;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::compute_saddle_point_model(
    const double& step, const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const
{
  return compute_saddle_point_model(step, step, dir, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::NLN::MeritFunction::Lagrangian::name() const
{
  return merit_function_name_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::map<std::string, NOX::NLN::MeritFunction::MeritFctName>
NOX::NLN::MeritFunction::Lagrangian::get_supported_type_list() const
{
  std::map<std::string, MeritFctName> type_names;

  type_names["Lagrangian"] = mrtfct_lagrangian;
  type_names["Lagrangian Active"] = mrtfct_lagrangian_active;

  return type_names;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Lagrangian::set_type(const std::string& identifier)
{
  static const std::map<std::string, MeritFctName> supported_type_names = get_supported_type_list();

  auto cit = supported_type_names.cbegin();
  while (cit != supported_type_names.cend())
  {
    if (identifier == cit->first)
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

    FOUR_C_THROW("Unknown type name: \"%s\"", identifier.c_str());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Lagrangian::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_->isPrintType(::NOX::Utils::Error))
  {
    utils_->out() << "ERROR - NOX::NLN::MeritFunction::Lagrangian::" << functionName << " - "
                  << errorMsg << std::endl;
  }
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
