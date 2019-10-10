
/*! \file

 \brief evaluate routines for the artery element

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/

#include "artery.H"
#include "artery_ele_action.H"
#include "artery_ele_factory.H"
#include "artery_ele_interface.H"
#include "artery_ele_calc_lin_exp.H"
#include "../drt_inpar/inpar_bio.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/matlist.H"


#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;



/*---------------------------------------------------------------------*
 //evaluate the element (public)                            ismail 06/09
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // check for the action parameter
  const ARTERY::Action action = DRT::INPUT::get<ARTERY::Action>(params, "action");
  /*
  Here must add the steps for evaluating an element
  */
  Teuchos::RCP<MAT::Material> mat = Material();

  switch (action)
  {
    case ARTERY::calc_sys_matrix_rhs:
    {
      return DRT::ELEMENTS::ArtNetFactory::ProvideImpl(Shape(), impltype_, discretization.Name())
          ->Evaluate(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3, mat);
    }
    break;
    case ARTERY::calc_scatra_sys_matrix_rhs:
    {
      return DRT::ELEMENTS::ArtNetFactory::ProvideImpl(Shape(), impltype_, discretization.Name())
          ->ScatraEvaluate(this, params, discretization, la[0].lm_, elemat1, elemat2, elevec1,
              elevec2, elevec3, mat);
      break;
    }
    case ARTERY::get_initial_artery_state:
    case ARTERY::set_term_bc:
    case ARTERY::set_scatra_term_bc:
    case ARTERY::set_scatra_bc:
    case ARTERY::solve_riemann_problem:
    case ARTERY::calc_postpro_vals:
    case ARTERY::calc_scatra_from_scatra_fb:
    case ARTERY::evaluate_wf_wb:
    case ARTERY::evaluate_scatra_analytically:
    case ARTERY::calc_flow_pressurebased:
    {
      return DRT::ELEMENTS::ArtNetFactory::ProvideImpl(Shape(), impltype_, discretization.Name())
          ->EvaluateService(this, action, params, discretization, la, elemat1, elemat2, elevec1,
              elevec2, elevec3, mat);
    }
    break;
    default:
      dserror("Unkown type of action %d for Artery", action);
  }  // end of switch(act)


  return 0;
}  // end of DRT::ELEMENTS::Artery::Evaluate


int DRT::ELEMENTS::Artery::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::EvaluateDirichlet(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
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
  }
  return hoel;
}
