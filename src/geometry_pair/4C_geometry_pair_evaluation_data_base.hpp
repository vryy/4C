// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_EVALUATION_DATA_BASE_HPP
#define FOUR_C_GEOMETRY_PAIR_EVALUATION_DATA_BASE_HPP


#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * \brief A base class that all geometry pair evaluation data container have to inherit from.
   */
  class FOUR_C_API(FOUR_C_CORE) GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryEvaluationDataBase() = default;

    /**
     * \brief Destructor.
     */
    virtual ~GeometryEvaluationDataBase() = default;

    /**
     * \brief Clear data in this evaluation container.
     */
    virtual void clear() {};
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
