/*---------------------------------------------------------------------*/
/*! \file

\brief L2 projection on nodes

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_L2_PROJECTION_HPP
#define FOUR_C_FEM_GENERAL_L2_PROJECTION_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::FE
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
  Teuchos::RCP<Epetra_MultiVector> compute_nodal_l2_projection(
      Teuchos::RCP<Core::FE::Discretization> dis,  ///< underlying discretization
      const std::string& statename,                ///< name of state which will be set
      const int& numvec,                           ///< number of entries per node to project
      Teuchos::ParameterList& params,  ///< parameter list that contains the element action
      const Teuchos::ParameterList&
          solverparams,  ///< solver parameters for solving the resulting global system;
      const std::function<const Teuchos::ParameterList&(int)>
          get_solver_params  ///< function that returns the solver parameters for the i-th solver
  );

  Teuchos::RCP<Epetra_MultiVector> evaluate_and_solve_nodal_l2_projection(
      Core::FE::Discretization& dis, const Epetra_Map& noderowmap, const std::string& statename,
      const int& numvec, Teuchos::ParameterList& params, const Teuchos::ParameterList& solverparams,
      const std::function<const Teuchos::ParameterList&(int)> get_solver_params,
      const Epetra_Map& fullnoderowmap, const std::map<int, int>& slavetomastercolnodesmap);

  Teuchos::RCP<Epetra_MultiVector> solve_nodal_l2_projection(Core::LinAlg::SparseMatrix& massmatrix,
      Epetra_MultiVector& rhs, const Epetra_Comm& comm, const int& numvec,
      const Teuchos::ParameterList& solverparams,
      const std::function<const Teuchos::ParameterList&(int)> get_solver_params,
      const Epetra_Map& noderowmap, const Epetra_Map& fullnoderowmap,
      const std::map<int, int>& slavetomastercolnodesmap);

}  // namespace Core::FE


FOUR_C_NAMESPACE_CLOSE

#endif
