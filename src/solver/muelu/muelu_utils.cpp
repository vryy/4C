/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for interaction with MueLu

\level 3
*/
/*----------------------------------------------------------------------*/
#include "muelu_utils.H"

#include "../../drt_lib/drt_dserror.H"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_UseShortNames.hpp>

Teuchos::RCP<MultiVector> LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromMLList(
    const Teuchos::RCP<const Map>& rowMap, Teuchos::ParameterList& mllist)
{
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;

  // Extract info about nullspace dimension
  if (!mllist.isParameter("null space: dimension"))
    dserror("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int nspDimension = mllist.get<int>("null space: dimension");
  if (nspDimension < 1)
    dserror("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  // Create an Xpetra::MultiVector, where the i-th column will then be filled with the i-th
  // nullspace vector
  RCP<MultiVector> nullspace = MultiVectorFactory::Build(rowMap, nspDimension, true);
  RCP<std::vector<double>> nsdata =
      mllist.get<RCP<std::vector<double>>>("nullspace", Teuchos::null);
  for (size_t dim = 0; dim < Teuchos::as<size_t>(nspDimension); ++dim)
  {
    ArrayRCP<Scalar> nspVectorData = nullspace->getDataNonConst(dim);
    const LocalOrdinal myLength = nullspace->getLocalLength();
    for (LocalOrdinal dofLID = 0; dofLID < myLength; ++dofLID)
      nspVectorData[dofLID] = (*nsdata)[dim * myLength + dofLID];
  }

  return nullspace;
}
