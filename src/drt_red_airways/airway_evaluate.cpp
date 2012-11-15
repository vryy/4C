
/*!----------------------------------------------------------------------
\file airway_evaluate.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            Ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/


#include "red_airway.H"
#include "airway_impl.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"

#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;



/*---------------------------------------------------------------------*
 |evaluate the element (public)                            ismail 01/10|
 *---------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirway::Evaluate(ParameterList& params,
                                       DRT::Discretization&      discretization,
                                       vector<int>&              lm,
                                       Epetra_SerialDenseMatrix& elemat1,
                                       Epetra_SerialDenseMatrix& elemat2,
                                       Epetra_SerialDenseVector& elevec1,
                                       Epetra_SerialDenseVector& elevec2,
                                       Epetra_SerialDenseVector& elevec3)
{

  DRT::ELEMENTS::RedAirway::ActionType act = RedAirway::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
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
  else
  {

    char errorout[200];
    sprintf(errorout,"Unknown type of action (%s) for reduced dimensional airway",action.c_str());

    dserror(errorout);
  }

/*
Here must add the steps for evaluating an element
*/
  RCP<MAT::Material> mat = Material();

  switch(act)
  {
  case calc_sys_matrix_rhs:
  {
    return DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->Evaluate(this,
                                                                       params,
                                                                       discretization,
                                                                       lm,
                                                                       elemat1,
                                                                       elemat2,
                                                                       elevec1,
                                                                       elevec2,
                                                                       elevec3,
                                                                       mat);
  }
  break;
  case calc_sys_matrix_rhs_iad:
  {
  }
  break;
  case get_initial_state:
  {
    DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->Initial(this,
                                                               params,
                                                               discretization,
                                                               lm,
                                                               mat);
    
  }
  break;
  case set_bc:
  {
    DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->EvaluateTerminalBC(this,
                                                                          params,
                                                                          discretization,
                                                                          lm,
                                                                          elevec1,
                                                                          mat);
    
  }
  break;
  case calc_flow_rates:
  {
    DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->CalcFlowRates(this,
                                                                     params,
                                                                     discretization,
                                                                     elevec1,
                                                                     elevec2,
                                                                     lm,
                                                                     mat);

  }
  break;
  case get_coupled_values:
  {
    DRT::ELEMENTS::RedAirwayImplInterface::Impl(this)->GetCoupledValues(this,
                                                                        params,
                                                                        discretization,
                                                                        lm,
                                                                        mat);
    
  }
  break;
  default:
    dserror("Unkown type of action for reduced dimensional airways");
  }// end of switch(act)
  
  return 0;
} // end of DRT::ELEMENTS::RedAirway::Evaluate


int DRT::ELEMENTS::RedAirway::EvaluateNeumann(ParameterList& params,
                                              DRT::Discretization& discretization,
                                              DRT::Condition& condition,
                                              vector<int>& lm,
                                              Epetra_SerialDenseVector& elevec1,
                                              Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/10|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RedAirway::EvaluateDirichlet(ParameterList& params,
                                                DRT::Discretization&      discretization,
                                                DRT::Condition&           condition,
                                                vector<int>&              lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


// get optimal gaussrule for discretization type
GaussRule1D DRT::ELEMENTS::RedAirway::getOptimalGaussrule(const DiscretizationType& distype)
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
    }
  return rule;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::ELEMENTS::RedAirway::isHigherOrderElement(
  const DRT::Element::DiscretizationType  distype) const
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
  }
  return hoel;
}

