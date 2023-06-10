/*----------------------------------------------------------------------*/
/*! \file

\brief Computation of specific solver parameters

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef LINEAR_SOLVER_METHOD_PARAMETERS_H
#define LINEAR_SOLVER_METHOD_PARAMETERS_H

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
        ::DRT::Discretization& dis, Teuchos::ParameterList& solverlist);
  };
}  // namespace CORE::LINEAR_SOLVER

#endif
