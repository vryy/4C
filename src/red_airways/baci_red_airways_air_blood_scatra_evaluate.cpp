/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Templated Evaluate file for blood scatra element containing the
action types for a reduced blood scatra element RedAirBloodScatra. The actual implementation of the
routines called during the possible actions is contained in red_air_blood_scatra_impl.cpp


\level 3

*/
/*---------------------------------------------------------------------*/


#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_lib_discret.H"
#include "baci_lib_exporter.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_list.H"
#include "baci_mat_newtonianfluid.H"
#include "baci_red_airways_air_blood_scatra_impl.H"
#include "baci_red_airways_elementbase.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_SerialDenseSolver.hpp>

using namespace DRT::UTILS;


/*---------------------------------------------------------------------*
 |evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirBloodScatra::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::RedAirBloodScatra::ActionType act = RedAirBloodScatra::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedAirBloodScatra::calc_sys_matrix_rhs;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedAirBloodScatra::calc_sys_matrix_rhs_iad;
  else if (action == "get_initial_state")
    act = RedAirBloodScatra::get_initial_state;
  else if (action == "set_bc")
    act = RedAirBloodScatra::set_bc;
  else if (action == "calc_flow_rates")
    act = RedAirBloodScatra::calc_flow_rates;
  else if (action == "calc_elem_volumes")
    act = RedAirBloodScatra::calc_elem_volumes;
  else if (action == "get_coupled_values")
    act = RedAirBloodScatra::get_coupled_values;
  else if (action == "get_junction_volume_mix")
    act = RedAirBloodScatra::get_junction_volume_mix;
  else if (action == "solve_scatra")
    act = RedAirBloodScatra::solve_scatra;
  else if (action == "calc_cfl")
    act = RedAirBloodScatra::calc_cfl;
  else if (action == "solve_junction_scatra")
    act = RedAirBloodScatra::solve_junction_scatra;
  else if (action == "eval_nodal_essential_values")
    act = RedAirBloodScatra::eval_nodal_ess_vals;
  else if (action == "solve_blood_air_transport")
    act = RedAirBloodScatra::solve_blood_air_transport;
  else if (action == "update_scatra")
    act = RedAirBloodScatra::update_scatra;
  else if (action == "update_elem12_scatra")
    act = RedAirBloodScatra::update_elem12_scatra;
  else if (action == "eval_PO2_from_concentration")
    act = RedAirBloodScatra::eval_PO2_from_concentration;
  else
  {
    char errorout[200];
    sprintf(errorout, "Unknown type of action (%s) for reduced dimensional acinus", action.c_str());

    dserror(errorout);
  }

  /*
  Here must add the steps for evaluating an element
  */
  Teuchos::RCP<MAT::Material> mat = Material();

  switch (act)
  {
    case calc_sys_matrix_rhs:
    {
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      DRT::ELEMENTS::RedAirBloodScatraImplInterface::Impl(this)->Initial(
          this, params, discretization, lm, mat);
    }
    break;
    case set_bc:
    {
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
    case solve_junction_scatra:
    {
      // do nothing
    }
    break;
    case update_scatra:
    {
      // do nothing
    }
    break;
    case update_elem12_scatra:
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
    case solve_blood_air_transport:
    {
      DRT::ELEMENTS::RedAirBloodScatraImplInterface::Impl(this)->SolveBloodAirTransport(
          this, elevec1, elevec2, elevec3, params, discretization, lm, mat);
    }
    break;
    default:
      dserror("Unkown type of action for reduced dimensional acinuss");
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::RedAirBloodScatra::Evaluate


int DRT::ELEMENTS::RedAirBloodScatra::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirBloodScatra::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  return 0;
}


// get optimal gaussrule for discretization type
CORE::DRT::UTILS::GaussRule1D DRT::ELEMENTS::RedAirBloodScatra::getOptimalGaussrule(
    const CORE::FE::CellType& distype)
{
  CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::line2:
      rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
      break;
    case CORE::FE::CellType::line3:
      rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::ELEMENTS::RedAirBloodScatra::isHigherOrderElement(const CORE::FE::CellType distype) const
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
      dserror("distype unknown!");
  }
  return hoel;
}
