
/*!----------------------------------------------------------------------
\file artery_evaluate.cpp

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/
#ifdef D_ARTNET


#include "artery.H"
#include "artery_lin_exp.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/matlist.H"

#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;



/*---------------------------------------------------------------------*
 //evaluate the element (public)                            ismail 06/09
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Artery::ActionType act = Artery::none;

  // get the action required
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = Artery::calc_sys_matrix_rhs;
  else if (action == "get_initail_artery_state")
    act = Artery::get_initail_artery_state;
  else if (action == "solve_riemann_problem")
    act = Artery::solve_riemann_problem;
  else if (action == "set_term_bc")
    act = Artery::set_term_bc;
  else if (action == "set_scatra_term_bc")
    act = Artery::set_scatra_term_bc;
  else if (action == "calc_postprocessing_values")
    act = Artery::calc_postpro_vals;
  else if (action == "calc_scatra_sys_matrix_rhs")
    act = Artery::calc_scatra_sys_matrix_rhs;
  else if (action == "calc_scatra_from_scatra_fb")
     act = Artery::calc_scatra_from_scatra_fb;
  else if(action == "evaluate_wf_wb")
     act = Artery::evaluate_wf_wb;
  else if(action == "evaluate_scatra_analytically")
     act = Artery::evaluate_scatra_analytically;
  else
  {

    char errorout[200];
    sprintf(errorout,"Unknown type of action (%s) for 1D_Artery",action.c_str());

    dserror(errorout);
  }

/*
Here must add the steps for evaluating an element
*/
  Teuchos::RCP<MAT::Material> mat = Material();

  switch(act)
  {
    case calc_sys_matrix_rhs:
    {
      return DRT::ELEMENTS::ArteryExpInterface::Expl(this)->Evaluate(this,
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
    case get_initail_artery_state:
    {
      DRT::ELEMENTS::ArteryExpInterface::Expl(this)->Initial(this,
                                                             params,
                                                             discretization,
                                                             lm,
                                                             mat);

    }
    break;
    case set_term_bc:
    {
      DRT::ELEMENTS::ArteryExpInterface::Expl(this)->EvaluateTerminalBC(this,
                                                                        params,
                                                                        discretization,
                                                                        lm,
                                                                        mat);

    }
    break;
    case set_scatra_term_bc:
    {
      DRT::ELEMENTS::ArteryExpInterface::Expl(this)->EvaluateScatraBC(this,
                                                                        params,
                                                                        discretization,
                                                                        lm,
                                                                        mat);

    }
    break;
    case solve_riemann_problem:
    {
      DRT::ELEMENTS::ArteryExpInterface::Expl(this)->SolveRiemann(this,
                                                                  params,
                                                                  discretization,
                                                                  lm,
                                                                  mat);

    }
    break;
    case calc_postpro_vals:
    {
      DRT::ELEMENTS::ArteryExpInterface::Expl(this)->CalcPostprocessingValues(this,
                                                                              params,
                                                                              discretization,
                                                                              lm,
                                                                              mat);

    }
    break;
  case calc_scatra_sys_matrix_rhs:
  {
    return DRT::ELEMENTS::ArteryExpInterface::Expl(this)->ScatraEvaluate(this,
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
  case  calc_scatra_from_scatra_fb:
  {
    DRT::ELEMENTS::ArteryExpInterface::Expl(this)->CalcScatraFromScatraFW(this,
                                                                          params,
                                                                          discretization,
                                                                          lm,
                                                                          mat);
  }
  break;
  case evaluate_wf_wb:
  {
    DRT::ELEMENTS::ArteryExpInterface::Expl(this)->EvaluateWfAndWb(this,
                                                                   params,
                                                                   discretization,
                                                                   lm,
                                                                   mat);
  }
  break;
  case evaluate_scatra_analytically:
  {
    DRT::ELEMENTS::ArteryExpInterface::Expl(this)-> SolveScatraAnalytically(this,
                                                                   params,
                                                                   discretization,
                                                                   lm,
                                                                   mat);
  }
  break;
  default:
      dserror("Unkown type of action for Artery");
  }// end of switch(act)


  return 0;
} // end of DRT::ELEMENTS::Artery::Evaluate


int DRT::ELEMENTS::Artery::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    DRT::Condition& condition,
    std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::EvaluateDirichlet(Teuchos::ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


// get optimal gaussrule for discretization type
GaussRule1D DRT::ELEMENTS::Artery::getOptimalGaussrule(const DiscretizationType& distype)
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
bool DRT::ELEMENTS::Artery::isHigherOrderElement(
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


#endif  // #ifdef D_ARTNET
