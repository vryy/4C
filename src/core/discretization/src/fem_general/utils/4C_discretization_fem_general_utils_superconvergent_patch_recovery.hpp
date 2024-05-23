/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_SUPERCONVERGENT_PATCH_RECOVERY_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_SUPERCONVERGENT_PATCH_RECOVERY_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <random>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}


namespace CORE::FE
{
  /*!
    \brief reconstruct nodal values via superconvergent patch recovery

    \return an Epetra_MultiVector based on the discret's node row map containing numvec vectors
            with the reconstruced state
   */
  template <int dim>
  Teuchos::RCP<Epetra_MultiVector> compute_superconvergent_patch_recovery(
      DRT::Discretization& dis,       ///< underlying discretization
      const Epetra_Vector& state,     ///< state vector needed on element level
      const std::string& statename,   ///< name of state which will be set
      const int numvec,               ///< number of entries per node to project
      Teuchos::ParameterList& params  ///< parameter list that contains the element action
  );
}  // namespace CORE::FE


FOUR_C_NAMESPACE_CLOSE

#endif
