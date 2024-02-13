/*-----------------------------------------------------------*/
/*! \file

\brief Contains access traits definitions used by ArborX to access
       the baci-specific bounding volume data structure

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRIC_SEARCH_ACCESS_TRAITS_HPP
#define BACI_DISCRETIZATION_GEOMETRIC_SEARCH_ACCESS_TRAITS_HPP

#include "baci_config.hpp"

#include "baci_discretization_geometric_search_bounding_volume.hpp"

#ifdef BACI_WITH_ARBORX
#include <ArborX.hpp>

namespace ArborX
{
  using namespace BACI;

  template <>
  struct AccessTraits<std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>,
      PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector,
        std::size_t i)
    {
      return ArborX::Box{vector[i].second.bounding_volume_};
    }
  };

  template <>
  struct AccessTraits<std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>,
      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector,
        std::size_t i)
    {
      auto const& vector_entry = vector[i];
      return intersects(vector_entry.second.bounding_volume_);
    }
  };
}  // namespace ArborX

#endif

BACI_NAMESPACE_OPEN
BACI_NAMESPACE_CLOSE

#endif
