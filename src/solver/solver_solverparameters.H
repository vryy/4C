/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of nullspace for multigrid solvers

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERPARAMETERS_H_
#define SOLVER_SOLVERPARAMETERS_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_nullspace.H"

namespace LINALG
{
  namespace SOLVER
  {
    class Parameters
    {
     public:
      static void ComputeSolverParameters(DRT::Discretization& dis,
          Teuchos::ParameterList& solverlist, const int numdf, const int dimns)
      {
        // set nodal block information for solver
        {

        }

        // set specific information for the multigrid preconditioner ML / MueLu
        {
          solverlist.set("PDE equations", numdf);
          solverlist.set("null space: dimension", dimns);
          solverlist.set("null space: type", "pre-computed");
          solverlist.set("null space: add default vectors", false);

          // allocate dimns times the local length of the rowmap
          const int lrows = dis.DofRowMap(0)->NumMyElements();
          Teuchos::RCP<std::vector<double>> ns =
              Teuchos::rcp(new std::vector<double>(dimns * lrows));
          double* nullsp = &((*ns)[0]);

          solverlist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
          solverlist.set("null space: vectors", nullsp);
          solverlist.set<bool>("ML validate parameter list", false);  // TODO: What does this?

          // compute null space directly. that will call eletypes.
          LINALG::Nullspace::ComputeNullSpace(dis, ns, numdf, dimns);
        }
      }
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif