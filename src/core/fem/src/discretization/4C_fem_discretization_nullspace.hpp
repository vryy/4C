// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_DISCRETIZATION_NULLSPACE_HPP
#define FOUR_C_FEM_DISCRETIZATION_NULLSPACE_HPP

#include "4C_config.hpp"

#include "4C_linalg_multi_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;

  /*!
  \brief Compute the nullspace of a discretization

  This method looks in the solver parameters whether algebraic multigrid (AMG)
  is used as preconditioner. AMG desires the nullspace of the
  system of equations which is then computed here if it does not already exist
  in the parameter list.

  \param discretization (in): discretization to compute the nullspace for
  \param solveparams (in, out): List of parameters
  \param recompute (in)  : force method to recompute the nullspace
  */
  FOUR_C_API(FOUR_C_CORE)
  void compute_null_space_if_necessary(const Discretization& discretization,
      Teuchos::ParameterList& solveparams, bool recompute = false);

  /*!
   \brief Calculate the nullspace based on a given discretization

  The nullspace is build by looping over all nodes of a discretization and stored
          in the respective variable.

     \param dis (in): discretization
     \param dimns (in): nullspace dimension
     \param map (in): nullspace map
      */
  FOUR_C_API(FOUR_C_CORE)
  std::shared_ptr<Core::LinAlg::MultiVector<double>> compute_null_space(
      const Core::FE::Discretization& dis, const int dimns, const Core::LinAlg::Map& dofmap);
}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
