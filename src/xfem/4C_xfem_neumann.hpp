// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_XFEM_NEUMANN_HPP
#define FOUR_C_XFEM_NEUMANN_HPP


#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>



FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  template <typename T>
  class Vector;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace XFEM
{

  /// evaluate Neumann boundary conditions
  void evaluate_neumann(Teuchos::ParameterList& params,
      std::shared_ptr<Core::FE::Discretization> discret, Core::LinAlg::Vector<double>& systemvector,
      Core::LinAlg::SparseOperator* systemmatrix = nullptr);

  /// evaluate standard Neumann boundary conditions
  void evaluate_neumann_standard(
      std::multimap<std::string, Core::Conditions::Condition*>& condition, const double time,
      bool assemblemat, Teuchos::ParameterList& params, Core::FE::Discretization& discret,
      Core::LinAlg::Vector<double>& systemvector, Core::LinAlg::SparseOperator* systemmatrix);


}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
