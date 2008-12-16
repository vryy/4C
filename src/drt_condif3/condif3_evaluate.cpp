/*!
\file condif3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "condif3.H"
#include "condif3_utils.H"
#include "condif3_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // for CalculateFlux()
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <Epetra_SerialDenseSolver.h>

#include "../drt_lib/drt_globalproblem.H"


using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 06/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::Evaluate(ParameterList& params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: condif3 and condif2)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
      this,
      params,
      discretization,
      lm,
      elemat1,
      elemat2,
      elevec1,
      elevec2,
      elevec3
      );

} //DRT::ELEMENTS::Condif3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                        gjb 06/08|
 |                                                                      |
 |  The function is just a dummy. For the condif elements, the          |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::EvaluateNeumann(ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
