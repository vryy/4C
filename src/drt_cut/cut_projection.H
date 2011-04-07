#ifndef CUT_PROJECTION_H
#define CUT_PROJECTION_H

#include <set>
#include <vector>

#include "cut_point.H"

namespace GEO
{
namespace CUT
{
  class Mesh;
  class Element;
  class VolumeCell;
  class Facet;

class Hex8Projection
{
public:

  Hex8Projection( Mesh & mesh, Element * element, VolumeCell * cell, Point::PointPosition position, const std::set<Facet*> & facets, int axis, double r );

private:

  std::vector<Facet*> inner_facets_;
  std::vector<Point*> inner_points_;
  std::vector<Point*> projected_points_;
};

}
}

#endif
