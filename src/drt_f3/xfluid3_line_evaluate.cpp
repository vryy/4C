/*!----------------------------------------------------------------------
\file xfluid3_line_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3Line::Evaluate(   ParameterList& params,
                                    		DRT::Discretization&      discretization,
                                    		vector<int>&              lm,
                                   	 		Epetra_SerialDenseMatrix& elemat1,
                                    		Epetra_SerialDenseMatrix& elemat2,
                                    		Epetra_SerialDenseVector& elevec1,
                                    		Epetra_SerialDenseVector& elevec2,
                                    		Epetra_SerialDenseVector& elevec3)
{
	DRT::ELEMENTS::XFluid3Line::ActionType act = XFluid3Line::none;
	string action = params.get<string>("action","none");
	if (action == "none") dserror("No action supplied");
  	else dserror("Unknown type of action for XFluid3_Line");
  	
  	const DiscretizationType distype = this->Shape();
    switch(act)
    {
    default:
        dserror("Unknown type of action for Fluid3_Line");
  	} // end of switch(act)
    
  	return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XFluid3Line::EvaluateNeumann(	ParameterList& params,
                                           		DRT::Discretization&      discretization,
                                           		DRT::Condition&           condition,
                                           		vector<int>&              lm,
                                           		Epetra_SerialDenseVector& elevec1)
{  
 	dserror("Neumann condition on line not implemented");
	return 0;
}




#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
