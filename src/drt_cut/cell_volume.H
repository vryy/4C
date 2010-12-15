#ifndef CELL_VOLUME_H
#define CELL_VOLUME_H

#include <map>

#include <Teuchos_RCP.hpp>

#include "cell_facet.H"

namespace GEO
{
namespace CELL
{
  class FacetSet;

class Volume
{
public:

  FacetSet* NewFacetSet( int sid );

private:

  std::map<int, Teuchos::RCP<FacetSet> > facetsets_;

};

}
}

#endif
