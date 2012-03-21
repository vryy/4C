/*!-----------------------------------------------------------------------------------------------*
\file cut_meshhandle.cpp

\brief handle that holds the mesh specific information

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "cut_meshhandle.H"


GEO::CUT::SideHandle * GEO::CUT::MeshHandle::CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
#ifdef DRT_CUT_DUMPCREATION
  std::cout << "CreateSide( " << sid << ", ";
  std::copy( nids.begin(), nids.end(), std::ostream_iterator<int>( std::cout, ", " ) );
  std::cout << distype << " );\n";
#endif
  switch ( distype )
  {
//  case DRT::Element::quad4:
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
  case DRT::Element::quad4:
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
    case DRT::Element::quad4:
    {
      qsh = new Quad4SideHandle( mesh_, sid, nids );
      break;
    }
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
#ifdef DRT_CUT_DUMPCREATION
  std::cout << "CreateElement( " << eid << ", ";
  std::copy( nids.begin(), nids.end(), std::ostream_iterator<int>( std::cout, ", " ) );
  std::cout << distype << " );\n";
#endif
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


GEO::CUT::Node * GEO::CUT::MeshHandle::GetNode( int nid ) const
{
  return mesh_.GetNode( nid );
}

GEO::CUT::SideHandle * GEO::CUT::MeshHandle::GetSide( int eid ) const
{
  std::map<int, LinearSideHandle>::const_iterator i = linearsides_.find( eid );
  if ( i!=linearsides_.end() )
  {
    return const_cast<LinearSideHandle*>( &i->second );
  }
  std::map<int, Teuchos::RCP<QuadraticSideHandle> >::const_iterator j = quadraticsides_.find( eid );
  if ( j!=quadraticsides_.end() )
  {
    return &*j->second;
  }
  return NULL;
}

GEO::CUT::ElementHandle * GEO::CUT::MeshHandle::GetElement( int eid ) const
{
  std::map<int, LinearElementHandle>::const_iterator i = linearelements_.find( eid );
  if ( i!=linearelements_.end() )
  {
    return const_cast<LinearElementHandle*>( &i->second );
  }
  std::map<int, Teuchos::RCP<QuadraticElementHandle> >::const_iterator j = quadraticelements_.find( eid );
  if ( j!=quadraticelements_.end() )
  {
    return &*j->second;
  }
  return NULL;
}
