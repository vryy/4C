/*!----------------------------------------------------------------------
\file condif3_surface_evaluate.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element (tri or quad)

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>


*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "condif3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 06/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3Surface::Evaluate(ParameterList&            params,
                                            DRT::Discretization&      discretization,
                                            vector<int>&              lm,
                                            Epetra_SerialDenseMatrix& elemat1,
                                            Epetra_SerialDenseMatrix& elemat2,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseVector& elevec2,
                                            Epetra_SerialDenseVector& elevec3)
{
    DRT::ELEMENTS::Condif3Surface::ActionType act = Condif3Surface::none;
    string action = params.get<string>("action","none");
    if (action == "none") dserror("No action supplied");
    else dserror("Unknown type of action for Condif3_Surface");

    switch(act)
    {
    case Condif3Surface::none:
      dserror("action=none");
    default:
        dserror("Unknown type of action for Condif3_Surface");
    } // end of switch(act)

    return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3Surface::EvaluateNeumann(
                                           ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  dserror("EvaluateNeumann not implemented.");
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
