
#include "cut_pointcycle.H"
#include "cut_point_impl.H"
#include "cut_line.H"
#include "cut_mesh.H"

GEO::CUT::InnerCutPoints::InnerCutPoints( const std::map<Point*, std::vector<Point*>::const_iterator*> & point_positions,
                                          const std::vector<Point*> & points )
{
  std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i1 = point_positions.find( points.front() );
  std::map<Point*, std::vector<Point*>::const_iterator*>::const_iterator i2 = point_positions.find( points.back() );

  if ( i1==point_positions.end() or i2==point_positions.end() )
  {
    throw std::runtime_error( "inner cut list start point not in cut point list" );
  }

  points_.reserve( points.size() );
  if ( *i1->second <= *i2->second )
  {
    points_.assign( points.begin(), points.end() );
  }
  else
  {
    points_.assign( points.rbegin(), points.rend() );
  }
}

bool GEO::CUT::InnerCutPoints::operator==( const std::vector<Point*> & other ) const
{
  if ( points_.size()!=other.size() )
    return false;

  if ( points_.front()==other.front() )
  {
    return std::mismatch( points_.begin(), points_.end(), other.begin() ).first==points_.end();
  }
  else if ( points_.front()==other.back() )
  {
    return std::mismatch( points_.begin(), points_.end(), other.rbegin() ).first==points_.end();
  }
  else
  {
    return false;
  }
}

bool GEO::CUT::PointCycle::Split( const std::vector<Point*> & cut_points )
{
  if ( lhs_!=Teuchos::null /* or rhs_!=Teuchos::null */ )
  {
    bool lhs_split = lhs_->MyCut( cut_points );
    bool rhs_split = rhs_->MyCut( cut_points );

    if ( lhs_split and rhs_split )
    {
      throw std::runtime_error( "cannot split both lhs and rhs" );
    }
    if ( lhs_split )
    {
      return lhs_->Split( cut_points );
    }
    else if ( rhs_split )
    {
      return rhs_->Split( cut_points );
    }
    else
    {
      throw std::runtime_error( "split of neither side" );
    }
  }
  else
  {
    if ( NeedsSplit( cut_points ) )
    {
      std::vector<Point*>::iterator i1 = std::find( facet_points_.begin(), facet_points_.end(), cut_points.front() );
      std::vector<Point*>::iterator i2 = std::find( facet_points_.begin(), facet_points_.end(), cut_points.back() );

      if ( i1==facet_points_.end() or i2==facet_points_.end() )
      {
        throw std::runtime_error( "inner cut point not in facet point list" );
      }

      std::vector<Point*> s1;
      std::vector<Point*> s2;

      if ( i1 < i2 )
      {
        std::copy( facet_points_.begin(), i1, std::back_inserter( s1 ) );
        std::copy( cut_points.begin(), cut_points.end(), std::back_inserter( s1 ) );
        std::copy( i2+1, facet_points_.end(), std::back_inserter( s1 ) );

        std::copy( i1+1, i2, std::back_inserter( s2 ) );
        std::copy( cut_points.rbegin(), cut_points.rend(), std::back_inserter( s2 ) );
      }
      else if ( i1 > i2 )
      {
        std::copy( facet_points_.begin(), i2, std::back_inserter( s1 ) );
        std::copy( cut_points.rbegin(), cut_points.rend(), std::back_inserter( s1 ) );
        std::copy( i1+1, facet_points_.end(), std::back_inserter( s1 ) );

        std::copy( i2+1, i1, std::back_inserter( s2 ) );
        std::copy( cut_points.begin(), cut_points.end(), std::back_inserter( s2 ) );
      }
      else
      {
        std::copy( facet_points_.begin(), i1, std::back_inserter( s1 ) );
        std::copy( cut_points.begin(), cut_points.end(), std::back_inserter( s1 ) );
        std::copy( i1+1, facet_points_.end(), std::back_inserter( s1 ) );

        std::copy( cut_points.begin()+1, cut_points.end(), std::back_inserter( s2 ) );
      }

      lhs_ = Teuchos::rcp( new PointCycle( s1 ) );
      rhs_ = Teuchos::rcp( new PointCycle( s2 ) );
    }

    return true;
  }
}

bool GEO::CUT::PointCycle::NeedsSplit( const std::vector<Point*> & cut_points )
{
  for ( std::vector<Point*>::const_iterator i=cut_points.begin();
        i!=cut_points.end();
        ++i )
  {
    if ( std::find( facet_points_.begin(), facet_points_.end(), *i )!=facet_points_.end() )
    {
      return true;
    }
  }
  return false;
}

bool GEO::CUT::PointCycle::MyCut( const std::vector<Point*> & cut_points )
{
  std::vector<Point*>::iterator i1 = std::find( facet_points_.begin(), facet_points_.end(), cut_points.front() );
  std::vector<Point*>::iterator i2 = std::find( facet_points_.begin(), facet_points_.end(), cut_points.back() );

  return i1!=facet_points_.end() and i2!=facet_points_.end();
}

void GEO::CUT::PointCycle::CreateFacets( Mesh & mesh, Side * side, std::vector<Facet*> & facets )
{
  if ( lhs_!=Teuchos::null /* or rhs_!=Teuchos::null */ )
  {
    lhs_->CreateFacets( mesh, side, facets );
    rhs_->CreateFacets( mesh, side, facets );
  }
  else
  {
    if ( facet_points_.size() > 2 )
      facets.push_back( mesh.NewFacet( facet_points_, 0, side ) );
  }
}

GEO::CUT::PointCycleList::PointCycleList( Element * element,
                                          Side * side,
                                          const std::vector<Point*> & facet_points,
                                          const std::set<Point*, PointPidLess> & cut_points )
  : pcl_( facet_points )
{
  std::vector<std::vector<Point*>::const_iterator> ips;
  for ( std::set<Point*, PointPidLess>::const_iterator i=cut_points.begin();
        i!=cut_points.end();
        ++i )
  {
    Point * p = *i;
    std::vector<Point*>::const_iterator j = find( facet_points.begin(), facet_points.end(), p );
    if ( j==facet_points.end() )
    {
      throw std::runtime_error( "cut point not in facet point list" );
    }
    ips.push_back( j );
  }

  std::sort( ips.begin(), ips.end() );

  std::map<Point*, std::vector<Point*>::const_iterator*> point_positions;
  std::vector<InnerCutPoints> ipvec;

  for ( std::vector<std::vector<Point*>::const_iterator>::iterator i=ips.begin();
        i!=ips.end();
        ++i )
  {
    Point* start = **i;
    point_positions[start] = &*i;
  }

  for ( std::vector<std::vector<Point*>::const_iterator>::iterator i=ips.begin();
        i!=ips.end();
        ++i )
  {
    Point* start = **i;

    std::set<Line*> cut_lines;
    start->CutLines( side, cut_lines );

    for ( std::set<Line*>::iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
    {
      Line * line = *i;
      std::vector<Point*> ip;
      ip.push_back( start );
      Point* point;
      for ( point = line->OtherPoint( start );
            cut_points.count( point )==0;
            point = line->OtherPoint( point ) )
      {
        ip.push_back( point );
        line = point->CutLine( line, side, element );
        if ( line==NULL )
        {
          throw std::runtime_error( "no cut line found on side" );
        }
      }
      ip.push_back( point );

      if ( std::find( ipvec.begin(), ipvec.end(), ip )==ipvec.end() )
      {
        ipvec.push_back( InnerCutPoints( point_positions, ip ) );
      }
    }
  }

  std::sort( ipvec.begin(), ipvec.end(), InnerCutPointsLess( point_positions ) );

  for ( std::vector<InnerCutPoints>::iterator i=ipvec.begin(); i!=ipvec.end(); ++i )
  {
    InnerCutPoints & icp = *i;
    pcl_.Split( icp );
  }
}

void GEO::CUT::PointCycleList::CreateFacets( Mesh & mesh, Side * side, std::vector<Facet*> & facets )
{
  pcl_.CreateFacets( mesh, side, facets );
}
