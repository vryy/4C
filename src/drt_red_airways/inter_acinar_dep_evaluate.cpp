/*---------------------------------------------------------------------*/
/*! \file

\brief Templated Evaluate file for inter-acinar linker element containing
       the action types for a reduced inter-acinar linker element
       RedInterAcinarDep. The actual implementation of the routines called
       during the possible actions is contained in inter_acinar_dep_impl.cpp

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/


#include "red_airway.H"
#include "inter_acinar_dep_impl.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"

#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;


/*---------------------------------------------------------------------*
 |Evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedInterAcinarDep::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::RedInterAcinarDep::ActionType act = RedInterAcinarDep::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
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

    dserror(errorout);
  }

  /*
    Here one must add the steps for evaluating an element
  */
  Teuchos::RCP<MAT::Material> mat = Material();

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
      dserror("Unknown type of action for reduced dimensional acinuss");
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::RedInterAcinarDep::Evaluate


int DRT::ELEMENTS::RedInterAcinarDep::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedInterAcinarDep::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | Get optimal gaussrule for discretisation type                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
GaussRule1D DRT::ELEMENTS::RedInterAcinarDep::getOptimalGaussrule(const DiscretizationType& distype)
{
  DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
  switch (distype)
  {
    case line2:
      rule = DRT::UTILS::intrule_line_2point;
      break;
    case line3:
      rule = DRT::UTILS::intrule_line_3point;
      break;
    default:
      dserror("Unknown number of nodes for Gaussrule initialization in inter-acinar linker.");
      break;
  }
  return rule;
}


/*----------------------------------------------------------------------*
 | Check, whether higher order derivatives for shape functions          |
 | (dxdx, dxdy, ...) are necessary|                                     |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedInterAcinarDep::isHigherOrderElement(
    const DRT::Element::DiscretizationType distype) const
{
  bool hoel = true;
  switch (distype)
  {
    case line3:
      hoel = true;
      break;
    case line2:
      hoel = false;
      break;
    default:
      dserror("distype unknown!");
      break;
  }
  return hoel;
}
