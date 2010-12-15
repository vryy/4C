#ifndef CELL_CELL_H
#define CELL_CELL_H

#include <map>
#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>

#include "cell_facetset.H"
#include "cell_volume.H"

namespace GEO
{
  namespace CUT
  {
    class Point;
    class Facet;
    class LinearElement;
  }

namespace CELL
{

class Cell
{
public:

  Cell( GEO::CUT::LinearElement * element );

  void AddVolume( const std::set<GEO::CUT::Facet*> & collected_facets,
                  const std::map<std::pair<GEO::CUT::Point*, GEO::CUT::Point*>, std::set<GEO::CUT::Facet*> > & volume_lines );

  Point* NewPoint( const GEO::CUT::Point* p );

  Line* NewLine( Point* p1, Point* p2 );

  Facet* NewFacet( GEO::CUT::Facet * f, const std::vector<Point*> & ps );

  void NewFacets( GEO::CUT::Facet * f, std::set<Facet *> & nf );

  Volume* NewVolume();

private:

  GEO::CUT::LinearElement * element_;

  std::map<const GEO::CUT::Point*, Teuchos::RCP<Point> > points_;
  std::vector<Teuchos::RCP<Line> >     lines_;
  std::map<const GEO::CUT::Facet*, std::vector<Teuchos::RCP<Facet> > > facets_;
  std::vector<Teuchos::RCP<Volume> >   volumes_;

};

}
}

#endif
