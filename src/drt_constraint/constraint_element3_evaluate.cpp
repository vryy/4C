/*!----------------------------------------------------------------------
\file constraint_element3_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "constraint_element3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "Epetra_SerialDenseSolver.h"

using namespace DRT::UTILS;



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  return 0;
} // end of DRT::ELEMENTS::ConstraintElement3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      g.bau 07/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid2 elements, the          |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
