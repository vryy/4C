/*!----------------------------------------------------------------------
\file so3_thermo_evaluate.cpp
\brief

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so3_thermo.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include <iterator>

// TODO 2012-08-02 include headers of thermo-materials
//#include "../drt_mat/fluidporo.H"
//#include "../drt_mat/structporo.H"
//#include "../drt_mat/micromaterial.H"
//#include "../drt_mat/robinson.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_globalproblem.H"

// TODO 2012-08-02 check if it works if header is included at the beginning of
// the file, cf. poro at the end
#include "so3_thermo_fwd.hpp"


/*----------------------------------------------------------------------*
 | pre-evaluate the element (public)                         dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::PreEvaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la
  )
{
    return;
}

/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Thermo< so3_ele, distype>::Evaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{
  // start with "none"
  typename So3_Thermo::ActionType act = So3_Thermo::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");

  // TODO adopt actions from hex8_thermo_lin_evaluate
  else if (action=="calc_struct_stifftemp")   act = So3_Thermo::calc_struct_stifftemp;
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // coupling terms in force-vector and stiffness matrix
  case So3_Thermo::calc_struct_stifftemp:
  {
    CouplThrEvaluate(
      params,
      discretization,
      la,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }
  break;
  //==================================================================================
  default:
  {
    //in some cases we need to write/change some data before evaluating
    PreEvaluate(
      params,
      discretization,
      la
      );

    // call the purely structural methods
    so3_ele::Evaluate(
      params,
      discretization,
      la[0].lm_,  // only the first column, i.e. the structural field is passed
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );

    // add the temperature-dependent terms to the structural field
    CouplThrEvaluate(
      params,
      discretization,
      la,  // coupled TSI is considered, i.e. pass the compled location array
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }

  } // action

  return 0;
}  // Evaluate()


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 | here is the action for the coupling to the thermal field             |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::CouplThrEvaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{

  // start with "none"
  ActionType act = none;

  // TODO 2012-07-26 adopt action types
  // get the required action for coupling with the thermal field
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_update_istep")          act = calc_struct_update_istep;
  else if (action=="calc_struct_internalforce")         act = calc_struct_internalforce;
  else if (action=="calc_struct_nlnstiff")              act = calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")          act = calc_struct_nlnstiffmass;
  else if (action=="postprocess_stress")                act = postprocess_stress;
  else dserror("Unknown type of action for So3_Thermo: %s",action.c_str());
  // what should the element do
  switch(act)
  {
  //============================================================================
  // nonlinear stiffness, damping and internal force vector for poroelasticity
  case calc_struct_nlnstiff:
  {}
  break;

  //============================================================================
  // nonlinear stiffness, mass matrix and internal force vector for poroelasticity
  case calc_struct_nlnstiffmass:
  {}
  break;

  //============================================================================
  case calc_struct_update_istep:
  {}
  break;

  //============================================================================
  case postprocess_stress:
  {}
  break;

  //============================================================================
  default:
    dserror("Unknown type of action for So3_Thermo");
  } // action

  return 0;
}


/*----------------------------------------------------------------------*
 | initialise Jacobian                                       dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::InitJacobianMapping()
{
  //const static vector<LINALG::Matrix<numdim_,numnod_> > derivs;// = soh8_derivs();
  LINALG::Matrix<numdim_,numnod_> deriv ;
  LINALG::Matrix<numnod_,numdim_> xrefe;
  for (int i=0; i<numnod_; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  for (int gp=0; gp<numgpt_; ++gp)
  {
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim=0;idim<numdim_;idim++)
    {
       xsi_[gp](idim) = gpcoord[idim];
    }

    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    //invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(deriv,xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0",detJ_[gp]);
  }

  return;
}  // InitJacobianMapping


/*----------------------------------------------------------------------*/
// TODO 2012-08-02 check if header has to be included at the end of the file
//#include "so3_thermo_fwd.hpp"

