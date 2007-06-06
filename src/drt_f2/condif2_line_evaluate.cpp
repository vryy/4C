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
#ifdef TRILINOS_PACKAGE

#include "condif2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"


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
  DSTraceHelper dst("Condif2Line::EvaluateNeumann");
  dserror("Line Neumann condition not yet implemented for Condif2");

  return 0;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
