/*---------------------------------------------------------------------*/
/*! \file

\brief Templated Evaluate file for inter-acinar linker element containing
       the action types for a reduced inter-acinar linker element
       RedInterAcinarDep. The actual implementation of the routines called
       during the possible actions is contained in inter_acinar_dep_impl.cpp


\level 3

*/
/*---------------------------------------------------------------------*/


#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_red_airways_interacinardep_impl.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN



/*---------------------------------------------------------------------*
 |Evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedInterAcinarDep::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::RedInterAcinarDep::ActionType act = RedInterAcinarDep::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedInterAcinarDep::calc_sys_matrix_rhs;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedInterAcinarDep::calc_sys_matrix_rhs_iad;
  else if (action == "get_initial_state")
    act = RedInterAcinarDep::get_initial_state;
  else if (action == "set_bc")
    act = RedInterAcinarDep::set_bc;
  else if (action == "calc_flow_rates")
    act = RedInterAcinarDep::calc_flow_rates;
  else if (action == "calc_elem_volumes")
    act = RedInterAcinarDep::calc_elem_volumes;
  else if (action == "get_coupled_values")
    act = RedInterAcinarDep::get_coupled_values;
  else if (action == "get_junction_volume_mix")
    act = RedInterAcinarDep::get_junction_volume_mix;
  else if (action == "solve_scatra")
    act = RedInterAcinarDep::solve_scatra;
  else if (action == "calc_cfl")
    act = RedInterAcinarDep::calc_cfl;
  else if (action == "eval_nodal_essential_values")
    act = RedInterAcinarDep::eval_nodal_ess_vals;
  else if (action == "solve_blood_air_transport")
    act = RedInterAcinarDep::solve_blood_air_transport;
  else if (action == "update_scatra")
    act = RedInterAcinarDep::update_scatra;
  else if (action == "eval_PO2_from_concentration")
    act = RedInterAcinarDep::eval_PO2_from_concentration;
  else
  {
    char errorout[200];
    sprintf(
        errorout, "Unknown type of action (%s) for inter-acinar linker element", action.c_str());

    FOUR_C_THROW(errorout);
  }

  /*
    Here one must add the steps for evaluating an element
  */
  Teuchos::RCP<CORE::MAT::Material> mat = Material();

  switch (act)
  {
    case calc_sys_matrix_rhs:
    {
      return DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(this)->Evaluate(
          this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(this)->Initial(
          this, params, discretization, lm, elevec3, mat);
    }
    break;
    case set_bc:
    {
      DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(this)->EvaluateTerminalBC(
          this, params, discretization, lm, elevec1, mat);
    }
    break;
    case calc_flow_rates:
    {
    }
    break;
    case calc_elem_volumes:
    {
    }
    break;

    case get_coupled_values:
    {
      DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(this)->GetCoupledValues(
          this, params, discretization, lm, mat);
    }
    break;
    case get_junction_volume_mix:
    {
      // do nothing
    }
    break;
    case solve_scatra:
    {
      // do nothing
    }
    break;
    case calc_cfl:
    {
      // do nothing
    }
    break;
    case solve_blood_air_transport:
    {
      // do nothing
    }
    break;
    case eval_nodal_ess_vals:
    {
      // do nothing
    }
    break;
    case eval_PO2_from_concentration:
    {
      // do nothing
    }
    break;
    case update_scatra:
    {
      // do nothing
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of action for reduced dimensional acinuss");
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::RedInterAcinarDep::Evaluate


int DRT::ELEMENTS::RedInterAcinarDep::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedInterAcinarDep::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | Get optimal gaussrule for discretisation type                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule1D DRT::ELEMENTS::RedInterAcinarDep::getOptimalGaussrule(
    const CORE::FE::CellType& distype)
{
  CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::line2:
      rule = CORE::FE::GaussRule1D::line_2point;
      break;
    case CORE::FE::CellType::line3:
      rule = CORE::FE::GaussRule1D::line_3point;
      break;
    default:
      FOUR_C_THROW("Unknown number of nodes for Gaussrule initialization in inter-acinar linker.");
      break;
  }
  return rule;
}


/*----------------------------------------------------------------------*
 | Check, whether higher order derivatives for shape functions          |
 | (dxdx, dxdy, ...) are necessary|                                     |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedInterAcinarDep::is_higher_order_element(
    const CORE::FE::CellType distype) const
{
  bool hoel = true;
  switch (distype)
  {
    case CORE::FE::CellType::line3:
      hoel = true;
      break;
    case CORE::FE::CellType::line2:
      hoel = false;
      break;
    default:
      FOUR_C_THROW("distype unknown!");
      break;
  }
  return hoel;
}

FOUR_C_NAMESPACE_CLOSE
