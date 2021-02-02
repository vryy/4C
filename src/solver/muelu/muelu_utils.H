/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for interaction with MueLu

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELU_UTILS_H_
#define SOLVER_MUELU_UTILS_H_

#include <Kokkos_DefaultNode.hpp>

#include <MueLu_UseDefaultTypes.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>

namespace LINALG
{
  namespace SOLVER
  {
    namespace MUELU
    {
      namespace UTILS
      {
        /*!
        \brief Extract nullspace from Baci-provided ML list and convert to Xpetra::MultiVector

        \pre The input parameter list \c mllist needs to contain these entries:
        - "null space: dimension" (type: \c int )
        - "nullspace" (type: \c RCP<std::vector<double>> )

        @param[in] rowMap Map to be used to create the null space vectors
        @param[in] mllist Parameter list, where Baci has stored the null space data in ML style
        @return MultiVector with nullspace data
        */
        Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
        ExtractNullspaceFromMLList(
            const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
            Teuchos::ParameterList& mllist);
      }  // namespace UTILS
    }    // namespace MUELU
  }      // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_MUELU_UTILS_H_ */
