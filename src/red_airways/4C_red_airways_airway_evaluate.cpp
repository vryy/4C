/*---------------------------------------------------------------------*/
/*! \file

\brief Templated Evaluate file for airway element containing the action
       types for a reduced airway element RedAirway. The actual
       implementation of the routines called during the possible actions
       is contained in airway_impl.cpp


\level 3

*/
/*---------------------------------------------------------------------*/
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_list.hpp"
#include "4C_red_airways_airway_impl.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 |evaluate the element (public)                            ismail 01/10|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirway::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::RedAirway::ActionType act = RedAirway::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedAirway::calc_sys_matrix_rhs;
  else if (action == "get_initial_state")
    act = RedAirway::get_initial_state;
  else if (action == "set_bc")
    act = RedAirway::set_bc;
  else if (action == "calc_flow_rates")
    act = RedAirway::calc_flow_rates;
  else if (action == "get_coupled_values")
    act = RedAirway::get_coupled_values;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedAirway::calc_sys_matrix_rhs_iad;
  else if (action == "get_junction_volume_mix")
    act = RedAirway::get_junction_volume_mix;
  else if (action == "solve_scatra")
    act = RedAirway::solve_scatra;
  else if (action == "solve_junction_scatra")
    act = RedAirway::solve_junction_scatra;
  else if (action == "calc_cfl")
    act = RedAirway::calc_cfl;
  else if (action == "eval_nodal_essential_values")
    act = RedAirway::eval_nodal_ess_vals;
  else if (action == "solve_blood_air_transport")
    act = RedAirway::solve_blood_air_transport;
  else if (action == "update_scatra")
    act = RedAirway::update_scatra;
  else if (action == "update_elem12_scatra")
    act = RedAirway::update_elem12_scatra;
  else if (action == "eval_PO2_from_concentration")
    act = RedAirway::eval_PO2_from_concentration;
  else if (action == "calc_elem_volumes")
    act = RedAirway::calc_elem_volumes;
  else
  {
    char errorout[200];
    sprintf(errorout, "Unknown type of action (%s) for reduced dimensional airway", action.c_str());

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
      return DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->Evaluate(
          this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->Initial(
          this, params, discretization, lm, elevec1, elevec2, mat);
    }
    break;
    case set_bc:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->EvaluateTerminalBC(
          this, params, discretization, lm, elevec1, mat);
    }
    break;
    case calc_flow_rates:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->CalcFlowRates(
          this, params, discretization, lm, mat);
    }
    break;
    case get_coupled_values:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->GetCoupledValues(
          this, params, discretization, lm, mat);
    }
    break;
    case get_junction_volume_mix:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->get_junction_volume_mix(
          this, params, discretization, elevec1, lm, mat);
    }
    break;
    case solve_scatra:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->SolveScatra(
          this, params, discretization, elevec1, elevec2, lm, mat);
    }
    break;
    case solve_junction_scatra:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->solve_scatra_bifurcations(
          this, params, discretization, elevec1, elevec2, lm, mat);
    }
    break;
    case calc_cfl:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->CalcCFL(
          this, params, discretization, lm, mat);
    }
    break;
    case solve_blood_air_transport:
    {
      // do nothing
    }
    break;
    case eval_nodal_ess_vals:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->eval_nodal_essential_values(
          this, params, discretization, elevec1, elevec2, elevec3, lm, mat);
    }
    break;
    case eval_PO2_from_concentration:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->EvalPO2FromScatra(
          this, params, discretization, lm, mat);
    }
    break;
    case update_scatra:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->UpdateScatra(
          this, params, discretization, lm, mat);
    }
    break;
    case update_elem12_scatra:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->UpdateElem12Scatra(
          this, params, discretization, lm, mat);
    }
    break;
    case calc_elem_volumes:
    {
      DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->CalcElemVolume(
          this, params, discretization, lm, mat);
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of action for reduced dimensional airways");
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::RedAirway::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/10|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirway::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/10|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirway::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | get optimal gaussrule for discretisation type                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule1D DRT::ELEMENTS::RedAirway::getOptimalGaussrule(
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
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
      break;
  }
  return rule;
}


/*----------------------------------------------------------------------*
 | Check, whether higher order derivatives for shape functions          |
 | (dxdx, dxdy, ...) are necessary|                                     |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::is_higher_order_element(const CORE::FE::CellType distype) const
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
