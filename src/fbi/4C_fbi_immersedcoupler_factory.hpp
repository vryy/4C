// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_IMMERSEDCOUPLER_FACTORY_HPP
#define FOUR_C_FBI_IMMERSEDCOUPLER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace FBI
{
  class FBIGeometryCoupler;

  /**
   *  \brief Factory that creates the appropriate geometry coupler strategy
   *
   *  To create an immersed coupling implementation for FBI, call the static create_geometry_coupler
   * function directly! No instance of GeometryCouplerFactory has to be created! If you try to
   * call the constructor, you will get an error message.
   */
  class GeometryCouplerFactory
  {
    /// constructor
    GeometryCouplerFactory() = delete;

   public:
    /**
     *  \brief Creates the appropriate geometry coupler strategy
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] fbidyn List of FBI Input parameters
     *
     * \return FBI geometry coupler strategy
     */
    static std::shared_ptr<FBI::FBIGeometryCoupler> create_geometry_coupler(
        const Teuchos::ParameterList& fbidyn);
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
