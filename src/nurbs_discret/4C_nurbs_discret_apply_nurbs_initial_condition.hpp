/*----------------------------------------------------------------------*/
/*! \file

\brief A service method allowing the application of initial conditions
       for nurbs discretisations.

Since nurbs shape functions are not interpolating, it is not as
straightforward to apply initial conditions to the degrees of freedom.
(dofs are always associated with control points, i.e. the location
associated with the 'node'=control point is not the physical location
and the value at the control point is not the prescribed value at this
position since dofs associated with neighbouring control points influence
the function value as well)



\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_NURBS_DISCRET_APPLY_NURBS_INITIAL_CONDITION_HPP
#define FOUR_C_NURBS_DISCRET_APPLY_NURBS_INITIAL_CONDITION_HPP

#include "4C_config.hpp"

#include "4C_utils_function.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;
}
namespace Core::LinAlg
{
  class Solver;
}

namespace Discret
{
  namespace Nurbs
  {
    /*----------------------------------------------------------------------*/
    /*!
    \brief A service method allowing the application of initial conditions
           for nurbs discretisations. Recommended version with separate
           solver allocation

    \param dis            (i) the discretisation
    \param solverparams   (i) a list with solver parameters
    \param start_function (i) a function defining the initial field (i.e. u_0(x))
    \param initialvals    (o) the initial field on output (i.e. u_cp)

    \date 08/11
    */
    void apply_nurbs_initial_condition(Discret::Discretization& dis,
        const Teuchos::ParameterList& solverparams,
        const Core::UTILS::FunctionOfSpaceTime& start_function,
        Teuchos::RCP<Epetra_Vector> initialvals);

  }  // namespace Nurbs

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
