/*!----------------------------------------------------------------------
\file so_disp_line_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDispLine::Evaluate(	ParameterList& params,
                                    		DRT::Discretization&      discretization,
                                    		vector<int>&              lm,
                                   	 		Epetra_SerialDenseMatrix& elemat1,
                                    		Epetra_SerialDenseMatrix& elemat2,
                                    		Epetra_SerialDenseVector& elevec1,
                                    		Epetra_SerialDenseVector& elevec2,
                                    		Epetra_SerialDenseVector& elevec3)
{
	DRT::ELEMENTS::SoDispLine::ActionType act = SoDispLine::none;
	string action = params.get<string>("action","none");
	if (action == "none") dserror("No action supplied");
  	else dserror("Unknown type of action for SoDisp_Line");
  	
    switch(act)
    {
    default:
        dserror("Unknown type of action for SoDisp_Line");
  	} // end of switch(act)
    
  	return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition           gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDispLine::EvaluateNeumann(	ParameterList& params,
                                           		DRT::Discretization&      discretization,
                                           		DRT::Condition&           condition,
                                           		vector<int>&              lm,
                                           		Epetra_SerialDenseVector& elevec1)
{  
 	dserror("Neumann condition on line not implemented");
	return 0;
}




#endif  // #ifdef CCADISCRET
#endif // #ifdef
