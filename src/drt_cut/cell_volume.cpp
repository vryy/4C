
#include "cell_volume.H"
#include "cell_facetset.H"

GEO::CELL::FacetSet* GEO::CELL::Volume::NewFacetSet( int sid )
{
  std::map<int, Teuchos::RCP<FacetSet> >::iterator i = facetsets_.find( sid );
  if ( i!=facetsets_.end() )
  {
    return &*i->second;
  }
  FacetSet * fs = new FacetSet( sid );
  facetsets_[sid] = Teuchos::rcp( fs );
  return fs;
}
