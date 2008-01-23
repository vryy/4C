/*!----------------------------------------------------------------------
\file condif2_line_evaluate.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "condif2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 08/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Condif2Line::Evaluate(       ParameterList&            params,
                                                DRT::Discretization&      discretization,
                                                vector<int>&              lm,
                                                Epetra_SerialDenseMatrix& elemat1,
                                                Epetra_SerialDenseMatrix& elemat2,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseVector& elevec2,
                                                Epetra_SerialDenseVector& elevec3)
{
    DRT::Elements::Condif2Line::ActionType act = Condif2Line::none;
    string action = params.get<string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action == "integrate_Shapefunction")
        act = Condif2Line::integrate_Shapefunction;	
    else dserror("Unknown type of action for Condif2_Surface");
    
    switch(act)
    {
    case integrate_Shapefunction:
    {
        dserror("No implementation yet for convection-diffusion problems");
        break;
    }
    default:
        dserror("Unknown type of action for Condif2Line");
    } // end of switch(act)

    return 0;
    
} // DRT::Elements::Condif2Line::Evaluate


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)        vg 05/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Condif2Line::EvaluateNeumann(
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{  
  dserror("No implementation yet for convection-diffusion problems");
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
