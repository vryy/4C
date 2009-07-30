
/*!----------------------------------------------------------------------
\file artery_evaluate.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            Ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ARTNET
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "artery.H"
#include "artery_lin_exp.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../drt_mat/matlist.H"

#include <blitz/array.h>
#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;



/*---------------------------------------------------------------------*
 //evaluate the element (public)                            ismail 06/09
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Artery::ActionType act = Artery::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_sys_matrix_rhs")
    act = Artery::calc_sys_matrix_rhs;
  else if (action == "get_stiffness_matrix")
    act = Artery::get_stiffness_matrix;
  else if (action == "calc_stiffness_matrix")
    act = Artery::calc_stiffness_matrix;
  else if (action == "calc_charcteristic_variables")
    act = Artery::calc_charcteristic_variables;
  else if (action == "calc_mass_matrix")
    act = Artery::calc_mass_matrix;
  else if (action == "get_density")
    act = Artery::get_density;
  else if (action == "get_viscosity")
    act = Artery::get_viscosity;
  else if (action == "get_mass_matrix")
    act = Artery::get_mass_matrix;
  else if (action == "get_beta")
    act = Artery::get_beta;
  else if (action == "get_sound_speed")
    act = Artery::get_sound_speed;
  else if (action == "get_initail_artery_state")
    act = Artery::get_initail_artery_state;
  else if (action == "set_term_bc")
    act = Artery::set_term_bc;
  else
  {

    char errorout[200];
    sprintf(errorout,"Unknown type of action (%s) for 1D_Artery",action.c_str());

    dserror(errorout);
  }

/*
Here must add the steps for evaluating an element
*/
  RefCountPtr<MAT::Material> mat = Material();

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
    case calc_stiffness_matrix:
    {
    }
    break;
    case get_stiffness_matrix:
    {
      
    }
    break;
    case get_mass_matrix:
    {
    }
    break;
    case calc_mass_matrix:
    {
    }
    break;
    case calc_charcteristic_variables:
    {
    }
    break;
    case get_density:
    {
    }
    break;
    case get_viscosity:
    {
    }
    break;
    case get_beta:
    {
    }
    break;
    case get_sound_speed:
    {
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
                                                                        elemat1,
                                                                        elevec1,
                                                                        mat);

    }
    break;
    default:
      dserror("Unkown type of action for Artery");
  }// end of switch(act)




  return 0;
} // end of DRT::ELEMENTS::Artery::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     ismail 01/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Artery::EvaluateDirichlet(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
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


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                               ismail 07/09|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ArteryRegister::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ARTNET
