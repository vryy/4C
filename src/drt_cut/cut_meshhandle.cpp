
#include "cut_meshhandle.H"


GEO::CUT::SideHandle * GEO::CUT::MeshHandle::CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::quad4:
  case DRT::Element::tri3:
  {
    std::map<int, LinearSideHandle>::iterator i = linearsides_.find( sid );
    if ( i!=linearsides_.end() )
    {
      return &i->second;
    }

    Side * s = mesh_.CreateSide( sid, nids, distype );
    LinearSideHandle & lsh = linearsides_[sid];
    lsh = LinearSideHandle( s );
    return &lsh;
  }
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  case DRT::Element::tri6:
  {
    std::map<int, Teuchos::RCP<QuadraticSideHandle> >::iterator i = quadraticsides_.find( sid );
    if ( i!=quadraticsides_.end() )
    {
      return &*i->second;
    }

    QuadraticSideHandle * qsh = NULL;
    switch ( distype )
    {
    case DRT::Element::quad8:
    {
      qsh = new Quad8SideHandle( mesh_, sid, nids );
      break;
    }
    case DRT::Element::quad9:
    {
      qsh = new Quad9SideHandle( mesh_, sid, nids );
      break;
    }
    case DRT::Element::tri6:
    {
      qsh = new Tri6SideHandle( mesh_, sid, nids );
      break;
    }
    default:
      throw std::runtime_error( "unsupported distype" );
    }
    quadraticsides_[sid] = Teuchos::rcp( qsh );
    return qsh;
  }
  default:
    throw std::runtime_error( "unsupported distype" );
  }
}

GEO::CUT::ElementHandle * GEO::CUT::MeshHandle::CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::hex8:
  case DRT::Element::tet4:
  case DRT::Element::pyramid5:
  case DRT::Element::wedge6:
  {
    std::map<int, LinearElementHandle>::iterator i = linearelements_.find( eid );
    if ( i!=linearelements_.end() )
    {
      return &i->second;
    }

    Element * e = mesh_.CreateElement( eid, nids, distype );
    LinearElementHandle & leh = linearelements_[eid];
    leh = LinearElementHandle( e );
    return &leh;
  }
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
  {
    std::map<int, Teuchos::RCP<QuadraticElementHandle> >::iterator i = quadraticelements_.find( eid );
    if ( i!=quadraticelements_.end() )
    {
      return &*i->second;
    }

    QuadraticElementHandle * qeh = NULL;
    switch ( distype )
    {
    case DRT::Element::hex20:
    {
      qeh = new Hex20ElementHandle( mesh_, eid, nids );
      break;
    }
    case DRT::Element::hex27:
    {
      qeh = new Hex27ElementHandle( mesh_, eid, nids );
      break;
    }
    case DRT::Element::tet10:
    {
      qeh = new Tet10ElementHandle( mesh_, eid, nids );
      break;
    }
    case DRT::Element::wedge15:
    default:
      throw std::runtime_error( "unsupported distype" );
    }
    quadraticelements_[eid] = Teuchos::rcp( qeh );
    return qeh;
  }
  default:
    throw std::runtime_error( "unsupported distype" );
  }
}

GEO::CUT::SideHandle * GEO::CUT::MeshHandle::GetSide( int eid )
{
  std::map<int, LinearSideHandle>::iterator i = linearsides_.find( eid );
  if ( i!=linearsides_.end() )
  {
    return &i->second;
  }
  std::map<int, Teuchos::RCP<QuadraticSideHandle> >::iterator j = quadraticsides_.find( eid );
  if ( j!=quadraticsides_.end() )
  {
    return &*j->second;
  }
  return NULL;
}

GEO::CUT::ElementHandle * GEO::CUT::MeshHandle::GetElement( int eid )
{
  std::map<int, LinearElementHandle>::iterator i = linearelements_.find( eid );
  if ( i!=linearelements_.end() )
  {
    return &i->second;
  }
  std::map<int, Teuchos::RCP<QuadraticElementHandle> >::iterator j = quadraticelements_.find( eid );
  if ( j!=quadraticelements_.end() )
  {
    return &*j->second;
  }
  return NULL;
}

