#ifndef CELL_FACETSET_H
#define CELL_FACETSET_H

#include "cell_facet.H"

namespace GEO
{
namespace CELL
{

class FacetSet
{
public:

  FacetSet( int i ) : id_( i ) {}

  void AddFacets( const std::set<Facet *> & nf );

private:

  int id_;
  std::set<Facet*> facets_;
};

}
}

#endif
