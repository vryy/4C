/*---------------------------------------------------------------------*/
/*!

\brief Templated Evaluate file for acinus element containing the action
       types for a reduced acinus element RedAcinus. The actual
       implementation of the routines called during the possible actions
       is contained in acinus_impl.cpp

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/

#include "red_airway.H"
#include "acinus_impl.H"

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
 |evaluate the element (public)                            ismail 09/12|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAcinus::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::RedAcinus::ActionType act = RedAcinus::none;

  // get the action required
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = RedAcinus::calc_sys_matrix_rhs;
  else if (action == "calc_sys_matrix_rhs_iad")
    act = RedAcinus::calc_sys_matrix_rhs_iad;
  else if (action == "get_initial_state")
    act = RedAcinus::get_initial_state;
  else if (action == "set_bc")
    act = RedAcinus::set_bc;
  else if (action == "calc_flow_rates")
    act = RedAcinus::calc_flow_rates;
  else if (action == "calc_elem_volumes")
    act = RedAcinus::calc_elem_volumes;
  else if (action == "get_coupled_values")
    act = RedAcinus::get_coupled_values;
  else if (action == "get_junction_volume_mix")
    act = RedAcinus::get_junction_volume_mix;
  else if (action == "solve_scatra")
    act = RedAcinus::solve_scatra;
  else if (action == "solve_junction_scatra")
    act = RedAcinus::solve_junction_scatra;
  else if (action == "calc_cfl")
    act = RedAcinus::calc_cfl;
  else if (action == "eval_nodal_essential_values")
    act = RedAcinus::eval_nodal_ess_vals;
  else if (action == "solve_blood_air_transport")
    act = RedAcinus::solve_blood_air_transport;
  else if (action == "update_scatra")
    act = RedAcinus::update_scatra;
  else if (action == "update_elem12_scatra")
    act = RedAcinus::update_elem12_scatra;
  else if (action == "eval_PO2_from_concentration")
    act = RedAcinus::eval_PO2_from_concentration;
  else
  {
    char errorout[200];
    sprintf(errorout, "Unknown type of action (%s) for reduced dimensional acinus", action.c_str());

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
      return DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->Evaluate(
          this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case calc_sys_matrix_rhs_iad:
    {
    }
    break;
    case get_initial_state:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->Initial(
          this, params, discretization, lm, mat);
    }
    break;
    case set_bc:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->EvaluateTerminalBC(
          this, params, discretization, lm, elevec1, mat);
    }
    break;
    case calc_flow_rates:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->CalcFlowRates(
          this, params, discretization, lm, mat);
    }
    break;
    case calc_elem_volumes:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->CalcElemVolume(
          this, params, discretization, lm, mat);
    }
    break;
    case get_coupled_values:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->GetCoupledValues(
          this, params, discretization, lm, mat);
    }
    break;
    case get_junction_volume_mix:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->GetJunctionVolumeMix(
          this, params, discretization, elevec1, lm, mat);
    }
    break;
    case solve_scatra:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->SolveScatra(
          this, params, discretization, elevec1, elevec2, lm, mat);
    }
    break;
    case solve_junction_scatra:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->SolveScatraBifurcations(
          this, params, discretization, elevec1, elevec2, lm, mat);
    }
    break;
    case update_scatra:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->UpdateScatra(
          this, params, discretization, lm, mat);
    }
    break;
    case update_elem12_scatra:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->UpdateElem12Scatra(
          this, params, discretization, lm, mat);
    }
    break;
    case calc_cfl:
    {
      // do  nothing
    }
    break;
    case solve_blood_air_transport:
    {
      // do nothing
    }
    break;
    case eval_nodal_ess_vals:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->EvalNodalEssentialValues(
          this, params, discretization, elevec1, elevec2, elevec3, lm, mat);
    }
    break;
    case eval_PO2_from_concentration:
    {
      DRT::ELEMENTS::RedAcinusImplInterface::Impl(this)->EvalPO2FromScatra(
          this, params, discretization, lm, mat);
    }
    break;
    default:
      dserror("Unkown type of action for reduced dimensional acinuss");
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::RedAcinus::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 09/12|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAcinus::EvaluateNeumann(Teuchos::ParameterList& params,
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
int DRT::ELEMENTS::RedAcinus::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | get optimal gaussrule for discretisation type                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
GaussRule1D DRT::ELEMENTS::RedAcinus::getOptimalGaussrule(const DiscretizationType& distype)
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
      dserror("unknown number of nodes for gaussrule initialization");
      break;
  }
  return rule;
}


/*----------------------------------------------------------------------*
 | Check, whether higher order derivatives for shape functions          |
 | (dxdx, dxdy, ...) are necessary|                                     |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::isHigherOrderElement(
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
