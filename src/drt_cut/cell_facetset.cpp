
#include "cell_facetset.H"

void GEO::CELL::FacetSet::AddFacets( const std::set<Facet *> & nf )
{
  std::copy( nf.begin(), nf.end(), std::inserter( facets_, facets_.begin() ) );
}
