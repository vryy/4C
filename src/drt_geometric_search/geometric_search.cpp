/*----------------------------------------------------------------------*/
/*! \file

\brief Todo

\level 3

*/


#include "geometric_search.H"

#include "bounding_box.H"


void CollisionSearch(const std::vector<std::pair<int, BoundingBox>>& primitives,
    const std::vector<std::pair<int, BoundingBox>>& predicates, std::vector<int>& indices,
    std::vector<int>& offsets)
{
  indices.clear();
  offsets.clear();
  offsets.push_back(0);

  for (const auto& predicate : predicates)
  {
    int predicate_collisions = 0;
    for (size_t i_primitive = 0; i_primitive < primitives.size(); i_primitive++)
    {
      const auto& primitive = primitives[i_primitive];
      if (BoxesIntersect(primitive.second, predicate.second))
      {
        indices.push_back(i_primitive);
        predicate_collisions += 1;
      }
    }
    offsets.push_back(offsets.back() + predicate_collisions);
  }
}
