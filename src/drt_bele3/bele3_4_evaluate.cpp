/*!----------------------------------------------------------------------
\file bele3_evaluate.cpp
\brief

<pre>
 * Maintainer: Benedikt Schott
 *             schott@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "bele3_4.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/newtonianfluid.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           schott 11/11|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3_4::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                     schott 11/11|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3_4::EvaluateNeumann(Teuchos::ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

