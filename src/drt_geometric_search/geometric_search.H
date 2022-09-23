/*----------------------------------------------------------------------*/
/*! \file

\brief TODO

\level 3

*/
// End doxygen header.



#ifndef TODO_H_
#define TODO_H_


#include <vector>

struct BoundingBox;

// Todo: Put this in an namespace

// From arbor documentation:
// Finds all primitives meeting the predicates and record results in {indices, offsets}. indices
// stores the indices of the objects that satisfy the predicates. offsets stores the locations in
// the indices view that start a predicate, that is, predicates(i) is satisfied by
// primitives(indices(j)) for offsets(i) <= j < offsets(i+1). Following the usual convention,
// offsets(n) == indices.size(), where n is the number of queries that were performed and
// indices.size() is the total number of collisions.
void CollisionSearch(const std::vector<std::pair<int, BoundingBox>>& primitives,
    const std::vector<std::pair<int, BoundingBox>>& predicates, std::vector<int>& indices,
    std::vector<int>& offsets);



#endif
