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
  class IntegrationCell;

class Hex8Projection
{
public:

  static void HorizontalCut( Mesh & mesh,
                             Element * element,
                             VolumeCell * cell,
                             Point::PointPosition position,
                             const std::set<Facet*> & facets,
                             std::set<IntegrationCell*> & integrationcells,
                             int axis,
                             double r );

  static bool EdgeCut( Mesh & mesh,
                       Element * element,
                       VolumeCell * cell,
                       Point::PointPosition position,
                       const std::set<Facet*> & facets,
                       std::set<IntegrationCell*> & integrationcells,
                       int cutside1,
                       int cutside2,
                       int upside,
                       int downside );

  static void CreateTetMesh( Mesh & mesh, VolumeCell * cell, Point::PointPosition position, const std::set<Facet*> & myfacets );

private:

  static Facet * FindFacet( Element * element, const std::set<Facet*> & facets, int sideid );

  static void FindNeighborFacets( const std::vector<Point*> & facet_points,
                                  const std::set<Facet*> & facets,
                                  Facet * us,
                                  Facet * ds,
                                  std::set<Facet*> & myfacets );

};

}
}

#endif
