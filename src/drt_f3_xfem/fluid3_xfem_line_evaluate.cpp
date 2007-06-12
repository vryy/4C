/*!----------------------------------------------------------------------
\file fluid3_line_evaluate.cpp
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
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3Line::Evaluate(	ParameterList& params,
                                    		DRT::Discretization&      discretization,
                                    		vector<int>&              lm,
                                   	 		Epetra_SerialDenseMatrix& elemat1,
                                    		Epetra_SerialDenseMatrix& elemat2,
                                    		Epetra_SerialDenseVector& elevec1,
                                    		Epetra_SerialDenseVector& elevec2,
                                    		Epetra_SerialDenseVector& elevec3)
{
	DRT::Elements::XFluid3Line::ActionType act = XFluid3Line::none;
	string action = params.get<string>("action","none");
	if (action == "none") dserror("No action supplied");
	else if (action == "calc_ShapefunctDeriv1Deriv2")
  		act = XFluid3Line::calc_ShapefunctDeriv1Deriv2;
  	else dserror("Unknown type of action for XFluid3_Line");
  	
  	switch(act)
  	{
		case calc_ShapefunctDeriv1Deriv2:
      {
      	// functions, deriv1, deriv2 iel, r,          
      	f3_shapefunction_for_line(elevec1,elemat1,elemat2,lm[0],elevec2[0]); 											
      	dserror("Implement shapefunctions for XFluid3_Line");
      }
      default:
        dserror("Unknown type of action for XFluid3_Line");
  	} // end of switch(act)
  	return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3Line::EvaluateNeumann(	ParameterList& params,
                                           			DRT::Discretization&      discretization,
                                           			DRT::Condition&           condition,
                                           			vector<int>&              lm,
                                           			Epetra_SerialDenseVector& elevec1)
{  
 	dserror("Neumann condition on line not implemented");
	return 0;
}


/*
 * Line node numbering: Linear 
 * 
 * 
 */

void DRT::Elements::XFluid3Line::f3_shapefunction_for_line(	
											Epetra_SerialDenseVector&	funct ,
  											Epetra_SerialDenseMatrix& 	deriv1,
  											Epetra_SerialDenseMatrix& 	deriv2,
  											int&                 		iel   ,
  											double&              		r     )
{
    /*------------------------------- selection of polynomial interpolation */
    switch (iel)
    {
        case 2: /* LINEAR shape functions, 1.derivatives, 2.derivatives  ----*/
        {    
            double rp=1.0 + r;
    		double rm=1.0 - r;
    
    		funct[0]=0.5*rm;
    		funct[1]=0.5*rp;
    
    		deriv1(0,0)= -0.5;
     		deriv1(1,0)=  0.5;
     		
     		deriv2(0,0)= 0.0;
     		deriv2(1,0)= 0.0;
     		
    		break;
    	}
    	case 3: /* QUADRATIC shape functions, 1.derivatives, 2.derivatives  */
    	{
            double rp= 1.0 + r;
    		double rm= 1.0 - r;
    		double r2= 1.0 - r*r;

    		funct[0]= 0.5*r*rm;
    		funct[1]= r2;
    		funct[2]= 0.5*r*rp;
 
      	    deriv1(0,0)= r - 0.5;
      	    deriv1(1,0)= 1.0 - 2.0*r;
      	    deriv1(2,0)= r + 0.5;
      	
      	    deriv2(0,0)=  1.0;
            deriv2(1,0)= -2.0;
            deriv2(2,0)=  1.0;     	   
    		break;
    	}
    	default:
    		dserror("distyp unknown\n");
  } /* end switch(iel) */
 
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
