
#include "cut_integrationcell.H"

GEO::CUT::Hex8MainIntegrationCell::Hex8MainIntegrationCell( ConcreteElement<DRT::Element::hex8> * e )
  : element_( e )
{
  // Get all cut points.
  // Which sides are touched?

  std::set<Point*> cut_points;
  e->GetCutPoints( cut_points );
}

GEO::CUT::Tet4MainIntegrationCell::Tet4MainIntegrationCell( ConcreteElement<DRT::Element::tet4> * e )
  : element_( e )
{
}

GEO::CUT::Wedge6MainIntegrationCell::Wedge6MainIntegrationCell( ConcreteElement<DRT::Element::wedge6> * e )
  : element_( e )
{
}

GEO::CUT::Pyramid5MainIntegrationCell::Pyramid5MainIntegrationCell( ConcreteElement<DRT::Element::pyramid5> * e )
  : element_( e )
{
}
