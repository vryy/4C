// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_FACTORY_HPP
#define FOUR_C_GEOMETRY_PAIR_FACTORY_HPP


#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace Core::Elements
{
  class Element;
}
namespace GEOMETRYPAIR
{
  class GeometryPair;
  class GeometryEvaluationDataBase;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
  /**
   * \brief Create the correct geometry pair for line to volume coupling.
   * @return RCP to created geometry pair.
   */
  template <typename ScalarType, typename Line, typename Volume>
  Teuchos::RCP<GeometryPair> geometry_pair_line_to_volume_factory(
      const Core::Elements::Element* element1, const Core::Elements::Element* element2,
      const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data_ptr);

  /**
   * \brief Create the correct geometry pair for line to surface coupling.
   * @return RCP to created geometry pair.
   */
  template <typename ScalarType, typename Line, typename Surface>
  Teuchos::RCP<GeometryPair> geometry_pair_line_to_surface_factory(
      const Core::Elements::Element* element1, const Core::Elements::Element* element2,
      const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data_ptr);

  /**
   * \brief Create the correct geometry pair for line to surface coupling with FAD scalar types.
   *
   * The default geometry_pair_line_to_surface_factory would be sufficient for this, however, for
   * performance reasons it is better use the wrapped pairs, which are created in this function.
   *
   * @return RCP to created geometry pair.
   */
  template <typename ScalarType, typename Line, typename Surface>
  Teuchos::RCP<GeometryPair> geometry_pair_line_to_surface_factory_fad(
      const Core::Elements::Element* element1, const Core::Elements::Element* element2,
      const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data_ptr);
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
