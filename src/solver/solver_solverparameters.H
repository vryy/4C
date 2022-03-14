/*----------------------------------------------------------------------*/
/*! \file

\brief Computation of specific solver parameters

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERPARAMETERS_H_
#define SOLVER_SOLVERPARAMETERS_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_nullspace.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

namespace LINALG
{
  namespace SOLVER
  {
    class Parameters
    {
     public:
      /*!
        \brief Setting parameters related to specific solvers

        This method sets specific solver parameters as the nullspace,
        block information and  number of degrees of freedom.
      */

      static void ComputeSolverParameters(DRT::Discretization& dis,
          Teuchos::ParameterList& solverlist, const int numdf, const int dimns)
      {
        // set nodal block information for solver
        {
          solverlist.set("PDE equations", numdf);
          solverlist.set("null space: dimension", dimns);
          solverlist.set("null space: type", "pre-computed");
          solverlist.set("null space: add default vectors", false);
        }

        // set specific information for the multigrid preconditioner ML / MueLu
        {
          Teuchos::RCP<Epetra_MultiVector> nullspace =
              LINALG::NULLSPACE::ComputeNullSpace(dis, numdf, dimns);

          solverlist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
          solverlist.set("null space: vectors", nullspace->Values());
          solverlist.set<bool>("ML validate parameter list", false);
        }
      }
    };

  }  // namespace SOLVER
}  // namespace LINALG

#endif
