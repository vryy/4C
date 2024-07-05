/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of the infeasibility merit function for
       constrained problems. Especially useful for the filter method.



\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_solver_nonlin_nox_meritfunction_infeasibility.hpp"

#include "4C_solver_nonlin_nox_constraint_group.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::MeritFunction::Infeasibility::Infeasibility(
    const Teuchos::ParameterList& params, const ::NOX::Utils& u)
    : /* utils_( u ), */
      infeasibility_type_(mrtfct_vague)
{
  const std::string& type_name = params.get<std::string>("Type");
  set_type(type_name);

  merit_function_name_ = MeritFuncName2String(Type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::map<std::string, NOX::Nln::MeritFunction::MeritFctName>
NOX::Nln::MeritFunction::Infeasibility::get_supported_type_list() const
{
  std::map<std::string, MeritFctName> type_names;

  type_names["Two Norm"] = mrtfct_infeasibility_two_norm;
  type_names["Two Norm Active"] = mrtfct_infeasibility_two_norm_active;

  return type_names;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::MeritFunction::Infeasibility::set_type(const std::string& type_name)
{
  static const std::map<std::string, MeritFctName> supported_type_names = get_supported_type_list();

  auto cit = supported_type_names.cbegin();
  while (cit != supported_type_names.cend())
  {
    if (type_name == cit->first)
    {
      infeasibility_type_ = cit->second;
      break;
    }
    ++cit;
  }

  if (cit == supported_type_names.cend())
  {
    std::cout << "\n\n=====================================================\n";
    std::cout << "Supported infeasibility type names:\n"
                 "EXPECTED INPUT [= deduced merit function type]\n";
    for (const auto& supported_pair : supported_type_names)
      std::cout << supported_pair.first << " [= " << MeritFuncName2String(supported_pair.second)
                << "]\n";

    FOUR_C_THROW("Unknown type name: \"%s\"", type_name.c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::MeritFunction::Infeasibility::computef(const ::NOX::Abstract::Group& grp) const
{
  if (not grp.isF())
    FOUR_C_THROW(
        "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");

  // cast the nox-group to the constraint group
  const NOX::Nln::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::Nln::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr) FOUR_C_THROW("Dynamic cast to NOX::Nln::Constraint::Group failed!");

  return constr_grp_ptr->get_model_value(Type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::MeritFunction::Infeasibility::computeGradient(
    const ::NOX::Abstract::Group& group, ::NOX::Abstract::Vector& result) const
{
  FOUR_C_THROW("Currently unsupported.");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::MeritFunction::Infeasibility::computeSlope(
    const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    FOUR_C_THROW(
        "The current function value was not computed yet. Please call "
        "computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::Nln::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::Nln::CONSTRAINT::Group*>(&grp);
  if (not constr_grp_ptr) FOUR_C_THROW("Dynamic cast to NOX::Nln::Constraint::Group failed!");

  // compute the slope
  return constr_grp_ptr->get_linearized_model_terms(dir, Type(), linorder_first, lin_wrt_all_dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::MeritFunction::Infeasibility::computeQuadraticModel(
    const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const
{
  FOUR_C_THROW("Currently unsupported.");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::MeritFunction::Infeasibility::computeQuadraticMinimizer(
    const ::NOX::Abstract::Group& grp, ::NOX::Abstract::Vector& result) const
{
  FOUR_C_THROW("Currently unsupported.");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::Nln::MeritFunction::Infeasibility::name() const
{
  return merit_function_name_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::Nln::MeritFunction::MeritFctName NOX::Nln::MeritFunction::Infeasibility::Type() const
{
  return infeasibility_type_;
}

FOUR_C_NAMESPACE_CLOSE
