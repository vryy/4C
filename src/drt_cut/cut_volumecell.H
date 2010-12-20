#ifndef CUT_VOLUMECELL_H
#define CUT_VOLUMECELL_H

#include <map>
#include <vector>
#include <set>

namespace GEO
{
namespace CUT
{
  class Point;
  class Facet;
  class Element;
  class IntegrationCell;
  class Mesh;

class VolumeCell
{
public:

  VolumeCell( const std::set<Facet*> & facets,
              const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
              Element * element );

  bool Contains( Point * p );

  void Neighbors( Point * p,
                  const std::set<VolumeCell*> & cells,
                  std::set<VolumeCell*> & connected,
                  std::set<Element*> & elements );

  Element * ParentElement() { return element_; }

  void CreateIntegrationCells( Mesh & mesh );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  void ConnectNodalDOFSets();

private:

  Element * element_;

  std::set<Facet*> facets_;

  std::set<IntegrationCell*> integrationcells_;

  std::vector<int> nodaldofset_;
};

}
}

#endif
