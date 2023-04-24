/*----------------------------------------------------------------------*/
/*! \file

\brief Computation of specific solver parameters

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERPARAMETERS_H
#define SOLVER_SOLVERPARAMETERS_H

#include "lib_discret.H"
#include "lib_dserror.H"
#include "lib_elementtype.H"

#include "linalg_nullspace.H"
#include "linalg_utils_sparse_algebra_manipulation.H"

namespace LINALG::SOLVER
{
  class Parameters
  {
   public:
    /*!
      \brief Setting parameters related to specific solvers

      This method sets specific solver parameters as the nullspace,
      coordinates as well as block information and number of degrees of freedom.
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

      // set parameter information for solver
      {
        if (nullspaceMap == Teuchos::null and dis.NumMyRowNodes() > 0)
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
        std::array<int, 4> ldata{numdf, dimns, nv, np};
        std::array<int, 4> gdata{0, 0, 0, 0};
        dis.Comm().MaxAll(ldata.data(), gdata.data(), 4);
        numdf = gdata[0];
        dimns = gdata[1];
        nv = gdata[2];
        np = gdata[3];

        // store nullspace information in solver list
        solverlist.set("PDE equations", numdf);
        solverlist.set("null space: dimension", dimns);
        solverlist.set("null space: type", "pre-computed");
        solverlist.set("null space: add default vectors", false);
      }

      // set coordinate information
      {
        Teuchos::RCP<Epetra_MultiVector> coordinates;
        if (nullspaceMap == Teuchos::null)
          coordinates = dis.BuildNodeCoordinates();
        else
          coordinates = dis.BuildNodeCoordinates(nullspaceMap);

        solverlist.set<Teuchos::RCP<Epetra_MultiVector>>("Coordinates", coordinates);
      }

      // set nullspace information
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
}  // namespace LINALG::SOLVER

#endif
