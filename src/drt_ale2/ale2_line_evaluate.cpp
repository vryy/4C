/*!----------------------------------------------------------------------
\file ale2_line_evaluate.cpp
\brief

<pre>
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET

#include "ale2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Ale2Line::EvaluateNeumann(
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  dserror("Line Neumann condition not implemented for Ale2");

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_ALE2
