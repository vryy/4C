// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_CONSTRAINTENFORCER_FACTORY_HPP
#define FOUR_C_FBI_CONSTRAINTENFORCER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FBIConstraintenforcer;

  /**
   *  \brief Factory that creates the appropriate constraint enforcement strategy
   *
   *  To create a constraint enforcer implementation for FBI, call the static CreateEnforcer
   * function directly! No instance of ConstraintEnforcerFactory has to be created! If you try to
   * call the constructor, you will get an error message, since it is set to be private.
   */
  class ConstraintEnforcerFactory
  {
   private:
    /// constructor
    ConstraintEnforcerFactory() = delete;

   public:
    /**
     *  \brief Creates the appropriate constraint enforcement strategy
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] fsidyn List of FSI Input parameters
     * \param[in] fbidyn List of FBI Input parameters
     *
     * \return FBI constraint enforcement strategy
     */
    static std::shared_ptr<FBIConstraintenforcer> create_enforcer(
        const Teuchos::ParameterList& fsidyn, const Teuchos::ParameterList& fbidyn);
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
