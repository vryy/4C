/*----------------------------------------------------------------------*/
/*! \file

\brief Computation of specific solver parameters

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_method_parameters.hpp"

#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_discret_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

#include <Xpetra_EpetraIntMultiVector.hpp>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

void CORE::LINEAR_SOLVER::Parameters::compute_solver_parameters(
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
      CORE::Elements::Element* dwele = dis.lRowElement(0);
      dwele->ElementType().nodal_block_information(dwele, numdf, dimns, nv, np);
    }
    else
    {
      // if a map is given, grab the block information of the first element in that map
      for (int i = 0; i < dis.NumMyRowNodes(); ++i)
      {
        CORE::Nodes::Node* actnode = dis.lRowNode(i);
        std::vector<int> dofs = dis.Dof(0, actnode);

        const int localIndex = nullspaceMap->LID(dofs[0]);

        if (localIndex == -1) continue;

        CORE::Elements::Element* dwele = dis.lRowElement(localIndex);
        actnode->Elements()[0]->ElementType().nodal_block_information(dwele, numdf, dimns, nv, np);
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
      coordinates = dis.build_node_coordinates();
    else
      coordinates = dis.build_node_coordinates(nullspaceMap);

    solverlist.set<Teuchos::RCP<Epetra_MultiVector>>("Coordinates", coordinates);
  }

  // set nullspace information
  {
    if (nullspaceMap == Teuchos::null)
    {
      // if no map is given, we calculate the nullspace on the map describing the
      // whole discretization
      nullspaceMap = Teuchos::rcp(new Epetra_Map(*dis.dof_row_map()));
    }

    auto nullspace = DRT::ComputeNullSpace(dis, numdf, dimns, nullspaceMap);

    solverlist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
    solverlist.set("null space: vectors", nullspace->Values());
    solverlist.set<bool>("ML validate parameter list", false);
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::Parameters::FixNullSpace(std::string field, const Epetra_Map& oldmap,
    const Epetra_Map& newmap, Teuchos::ParameterList& solveparams)
{
  if (!oldmap.Comm().MyPID()) printf("Fixing %s Nullspace\n", field.c_str());

  // there is no ML or MueLu list, do nothing
  if (!solveparams.isSublist("ML Parameters") && !solveparams.isSublist("MueLu Parameters") &&
      !solveparams.isSublist("MueLu (FSI) Parameters"))
    return;

  // find the ML or MueLu list
  Teuchos::ParameterList* params_ptr = nullptr;
  if (solveparams.isSublist("ML Parameters"))
    params_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
    params_ptr = &(solveparams.sublist("MueLu Parameters"));
  else
    params_ptr = &(solveparams);
  Teuchos::ParameterList& params = *params_ptr;

  const int ndim = params.get("null space: dimension", -1);
  if (ndim == -1) FOUR_C_THROW("List does not contain nullspace dimension");

  Teuchos::RCP<Epetra_MultiVector> nullspace =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
  if (nullspace == Teuchos::null) FOUR_C_THROW("List does not contain nullspace");

  const int nullspaceLength = nullspace->MyLength();
  const int newmapLength = newmap.NumMyElements();

  if (nullspaceLength == newmapLength) return;
  if (nullspaceLength != oldmap.NumMyElements())
    FOUR_C_THROW("Nullspace map of length %d does not match old map length of %d", nullspaceLength,
        oldmap.NumMyElements());
  if (newmapLength > nullspaceLength)
    FOUR_C_THROW("New problem size larger than old - full rebuild of nullspace neccessary");

  Teuchos::RCP<Epetra_MultiVector> nullspaceNew =
      Teuchos::rcp(new Epetra_MultiVector(newmap, ndim, true));

  for (int i = 0; i < ndim; i++)
  {
    Epetra_Vector* nullspaceData = (*nullspace)(i);
    Epetra_Vector* nullspaceDataNew = (*nullspaceNew)(i);
    const int myLength = nullspaceDataNew->MyLength();

    for (int j = 0; j < myLength; j++)
    {
      int gid = newmap.GID(j);
      int olid = oldmap.LID(gid);
      if (olid == -1) continue;
      (*nullspaceDataNew)[j] = (*nullspaceData)[olid];
    }
  }

  params.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspaceNew);
  params.set("null space: vectors", nullspaceNew->Values());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CORE::LINEAR_SOLVER::Parameters::extract_nullspace_from_parameterlist(
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& row_map,
    Teuchos::ParameterList& list)
{
  if (!list.isParameter("null space: dimension"))
    FOUR_C_THROW(
        "CORE::LINEAR_SOLVER::Parameters::extract_nullspace_from_parameterlist: Multigrid "
        "parameter "
        "'null space: dimension' missing  in solver parameter list.");

  const int nullspace_dimension = list.get<int>("null space: dimension");
  if (nullspace_dimension < 1)
    FOUR_C_THROW(
        "CORE::LINEAR_SOLVER::Parameters::extract_nullspace_from_parameterlist: Multigrid "
        "parameter "
        "'null space: dimension' wrong. It has to be > 0.");

  Teuchos::RCP<Epetra_MultiVector> nullspace_data =
      list.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
  if (nullspace_data.is_null())
    FOUR_C_THROW(
        "CORE::LINEAR_SOLVER::Parameters::extract_nullspace_from_parameterlist: Nullspace data is "
        "null.");

  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace =
      Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(nullspace_data));

  nullspace->replaceMap(row_map);

  return nullspace;
}

FOUR_C_NAMESPACE_CLOSE
