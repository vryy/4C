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
#include "../drt_lib/drt_elementtype.H"

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

      static void ComputeSolverParameters(
          DRT::Discretization& dis, Teuchos::ParameterList& solverlist)
      {
        Teuchos::RCP<Epetra_Map> nullspaceMap =
            solverlist.get<Teuchos::RCP<Epetra_Map>>("null space: map", Teuchos::null);

        int numdf = 1;
        int dimns = 1;
        int nv = 0;
        int np = 0;

        // set nodal block information for solver
        {
          if (nullspaceMap == Teuchos::null)
          {
            // no map given, just grab the block information on the first element that appears
            DRT::Element* dwele = dis.lRowElement(0);
            dwele->ElementType().NodalBlockInformation(dwele, numdf, dimns, nv, np);
          }
          else
          {
            // if a map is given, grab the block information of the first element in that map
            for (int i = 0; i < dis.NumMyRowNodes(); ++i)
            {
              DRT::Node* actnode = dis.lRowNode(i);
              std::vector<int> dofs = dis.Dof(0, actnode);

              const int localIndex = nullspaceMap->LID(dofs[0]);

              if (localIndex == -1) continue;

              DRT::Element* dwele = dis.lRowElement(localIndex);
              actnode->Elements()[0]->ElementType().NodalBlockInformation(
                  dwele, numdf, dimns, nv, np);
              break;
            }
          }

          // communicate data to procs without row element
          int ldata[4] = {numdf, dimns, nv, np};
          int gdata[4] = {0, 0, 0, 0};
          dis.Comm().MaxAll(&ldata[0], &gdata[0], 4);
          numdf = gdata[0];
          dimns = gdata[1];
          nv = gdata[2];
          np = gdata[3];

          solverlist.set("PDE equations", numdf);
          solverlist.set("null space: dimension", dimns);
          solverlist.set("null space: type", "pre-computed");
          solverlist.set("null space: add default vectors", false);
        }

        // set specific information for the multigrid preconditioner ML / MueLu
        {
          Teuchos::RCP<Epetra_MultiVector> nullspace = Teuchos::null;
          if (nullspaceMap == Teuchos::null)
          {
            // if no map is given, we calculate the nullspace on the map describing the
            // whole discretization
            nullspaceMap = Teuchos::rcp(new Epetra_Map(*dis.DofRowMap()));
            nullspace = LINALG::NULLSPACE::ComputeNullSpace(dis, numdf, dimns, nullspaceMap);
          }
          else
          {
            // if a map is given, we calculate the nullspace on that map
            nullspace = LINALG::NULLSPACE::ComputeNullSpace(dis, numdf, dimns, nullspaceMap);
          }

          solverlist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
          solverlist.set("null space: vectors", nullspace->Values());
          solverlist.set<bool>("ML validate parameter list", false);
        }
      }
    };

  }  // namespace SOLVER
}  // namespace LINALG

#endif
