/*----------------------------------------------------------------------*/
/*! \file

\brief Computation of specific solver parameters

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_METHOD_PARAMETERS_H
#define BACI_LINEAR_SOLVER_METHOD_PARAMETERS_H

#include <Epetra_Map.h>
#include <Teuchos_ParameterListAcceptor.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::LINEAR_SOLVER
{
  class Parameters
  {
   public:
    /*!
      \brief Setting parameters related to specific solvers

      This method sets specific solver parameters such as the nullspace vectors,
      coordinates as well as block information and number of degrees of freedom
      into the solver parameter list.
    */
    static void ComputeSolverParameters(
        BACI::DRT::Discretization& dis, Teuchos::ParameterList& solverlist);

    /*!
     * \brief Fix the nullspace to match a new given map
     *
     * The nullspace is looked for in the parameter list. If found, it is assumed that
     * it matches the oldmap. Then it is fixed to match the new map.
     *
     * \param field (in): field name (just used for output)
     * \param oldmap (in): row map of nullspace
     * \param newmap (in): row map of nullspace upon exit
     * \param solveparams (in): parameterlist including nullspace vector
     */
    static void FixNullSpace(std::string field, const Epetra_Map& oldmap, const Epetra_Map& newmap,
        Teuchos::ParameterList& solveparams);
  };
}  // namespace CORE::LINEAR_SOLVER

BACI_NAMESPACE_CLOSE

#endif
