/*!----------------------------------------------------------------------
\file fluid3_xfem_line_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3Line::Evaluate(	ParameterList& params,
                                    		DRT::Discretization&      discretization,
                                    		vector<int>&              lm,
                                   	 		Epetra_SerialDenseMatrix& elemat1,
                                    		Epetra_SerialDenseMatrix& elemat2,
                                    		Epetra_SerialDenseVector& elevec1,
                                    		Epetra_SerialDenseVector& elevec2,
                                    		Epetra_SerialDenseVector& elevec3)
{
	DRT::Elements::Fluid3Line::ActionType act = Fluid3Line::none;
	string action = params.get<string>("action","none");
	if (action == "none") dserror("No action supplied");
	else if (action == "calc_Shapefunction")
        act = Fluid3Line::calc_Shapefunction;
    else if (action == "calc_ShapeDeriv1")
        act = Fluid3Line::calc_ShapeDeriv1;
    else if (action == "calc_ShapeDeriv2")
        act = Fluid3Line::calc_ShapeDeriv2;
  	else dserror("Unknown type of action for Fluid3_Line");
  	
  	const DiscretizationType distype = this->Shape();
    switch(act)
    {
    case calc_Shapefunction:
    {
     	shape_function_1D(elevec1,elevec2[0],distype);
      break;
    }
    case calc_ShapeDeriv1:
    {
    	shape_function_1D_deriv1(elemat1,elevec2[0],distype);
     	break;
    }
    case calc_ShapeDeriv2:
    {
		shape_function_1D_deriv2(elemat2,elevec2[0],distype);
      break;
    }
    default:
        dserror("Unknown type of action for Fluid3_Line");
  	} // end of switch(act)
    
  	return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3Line::EvaluateNeumann(	ParameterList& params,
                                           		DRT::Discretization&      discretization,
                                           		DRT::Condition&           condition,
                                           		vector<int>&              lm,
                                           		Epetra_SerialDenseVector& elevec1)
{  
 	dserror("Neumann condition on line not implemented");
	return 0;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
