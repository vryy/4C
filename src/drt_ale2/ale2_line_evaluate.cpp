//-----------------------------------------------------------------------
/*!
\file ale2_line_evaluate.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE

#include "ale2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2Line::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  dserror("Line Neumann condition not implemented for Ale2");

  return 0;
}


#endif // #ifdef D_ALE2
