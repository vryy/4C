/*---------------------------------------------------------------------*/
/*! \file

\brief L2 projection on nodes

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_FEM_GENERAL_L2_PROJECTION_HPP
#define BACI_DISCRETIZATION_FEM_GENERAL_L2_PROJECTION_HPP

#include "baci_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
}

namespace CORE::FE
{
  /*!
    \brief compute L2 projection of a dof based field onto a node based field in a least
    squares sense.
    WARNING: Make sure to pass down a dofrowmap appropriate for your discretization.

    \return an Epetra_MultiVector based on the discret's node row map containing numvec vectors
            with the projected state

    \author Georg Hammerl
    \date 06/14
   */
  Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(
      Teuchos::RCP<DRT::Discretization> dis,  ///< underlying discretization
      const std::string& statename,           ///< name of state which will be set
      const int& numvec,                      ///< number of entries per node to project
      Teuchos::ParameterList& params,         ///< parameter list that contains the element action
      const Teuchos::ParameterList&
          solverparams);  ///< solver parameters for solving the resulting global system;

  Teuchos::RCP<Epetra_MultiVector> EvaluateAndSolveNodalL2Projection(DRT::Discretization& dis,
      const Epetra_Map& noderowmap, const std::string& statename, const int& numvec,
      Teuchos::ParameterList& params, const Teuchos::ParameterList& solverparams,
      const Epetra_Map* fullnoderowmap = nullptr,
      const std::map<int, int>* slavetomastercolnodesmap = nullptr);

  Teuchos::RCP<Epetra_MultiVector> SolveNodalL2Projection(CORE::LINALG::SparseMatrix& massmatrix,
      Epetra_MultiVector& rhs, const Epetra_Comm& comm, const int& numvec,
      const Teuchos::ParameterList& solverparams, const Epetra_Map& noderowmap,
      const Epetra_Map* fullnoderowmap = nullptr,
      const std::map<int, int>* slavetomastercolnodesmap = nullptr);

}  // namespace CORE::FE


BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_H
