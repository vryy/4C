/*!----------------------------------------------------------------------
\file fluid3_surface_evaluate.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3.H"
#include "linalg_utils.H"
#include "drt_utils.H"
#include "drt_discret.H"
#include "drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
}
#include "dstrc.H"


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  bauer 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3Surface::EvaluateNeumann(
                                           ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{  
  DSTraceHelper dst("Fluid3Surface::EvaluateNeumann");
  dserror("Surface Neumann condition not yet implemented for Fluid3");

  return 0;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
