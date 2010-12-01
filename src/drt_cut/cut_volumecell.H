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
  class LinearElement;
  class IntegrationCell;

class VolumeCell
{
public:

  VolumeCell( const std::set<Facet*> & facets,
              const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
              LinearElement * element );

private:

  LinearElement * element_;

  std::set<Facet*> facets_;
};

}
}

#endif
