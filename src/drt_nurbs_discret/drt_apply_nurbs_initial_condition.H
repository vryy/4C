/*----------------------------------------------------------------------*/
/*!
\file drt_apply_nurbs_initial_condition.H

\brief A service method allowing the application of initial conditions
       for nurbs discretisations.

Since nurbs shape functions are not interpolating, it is not as
straightforward to apply initial conditions to the degrees of freedom.
(dofs are always associated with control points, i.e. the location
associated with the 'node'=control point is not the physical location
and the value at the control point is not the prescribed value at this
position since dofs associated with neighbouring control points influence
the function value as well)


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_APPLY_NURBS_INITIAL_CONDITION_H
#define DRT_APPLY_NURBS_INITIAL_CONDITION_H

#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

namespace DRT
{

namespace NURBS
{
/*----------------------------------------------------------------------*/
/*!
\brief A service method allowing the application of initial conditions
       for nurbs discretisations.

This method provides the following:

Given an initial distribution u_0(x) of initial values (for example by a
spatial function), we compute control point values u_cp such that they
minimize the following least-squares problem:

<pre>
                  ||                                   || 2
                  || +----+                            ||
                  ||  \                                ||
             min  ||   +    N   (x) * u     -  u   (x) ||
                  ||  /      cp       - cp     - 0     ||
             u    || +----+                            ||
             - cp ||   cp                              ||
</pre>


This is equivalent to the solution of the following linear system:

<pre>

         /                         \               /                         \
        |    /                      |             |    /                      |
 +----+ |   |                       |             |   |                       |
  \     |   |                       |    dim      |   |            dim        |
   +    |   | N   (x) * N   (x) dx  | * u     =   |   | N   (x) * u   (x) dx  |
  /     |   |  cp        cp         |    cp       |   |  cp        0          |
 +----+ |   |    j         i        |      j      |   |    j                  |
   cp   |  /                        |             |  /                        |
     j   \                         /               \                         /

        |                           |             |                           |
        +---------------------------+             +---------------------------+

               M(assmatrix)                                   rhs

</pre>


\param dis         (i) the discretisation
\param solver      (i) a solver object for the least-squares problem
\param ndim        (i) problem dimension
\param startfuncno (i) the number of the startfunction defining
                       the initial field (i.e. u_0(x))
\param initialvals (o) the initial field on output (i.e. u_cp)

\date 04/09
*/
  void apply_nurbs_initial_condition(
    DRT::Discretization&        dis         ,
    LINALG::Solver&             solver      ,
    const int                   ndim        ,
    const int                   startfuncno ,
    Teuchos::RCP<Epetra_Vector> initialvals );


} // namespace NURBS

} //namespace DRT

#endif // #ifdef DRT_APPLY_NURBS_INITIAL_CONDITION_H
#endif // #ifdef CCADISCRET
