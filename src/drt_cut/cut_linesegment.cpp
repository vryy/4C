
#include "cut_linesegment.H"
#include "cut_line.H"
#include "cut_point.H"
#include "cut_point_impl.H"
#include "cut_edge.H"
#include "cut_mesh.H"
#include "cut_element.H"

GEO::CUT::LineSegment::LineSegment( Mesh & mesh, Element * element, Side * side, std::set<Line*> & cut_lines, bool inner )
{
  Line * line = *cut_lines.begin();
  Point * begin = line->BeginPoint();

  facet_lines_.push_back( line );
  facet_points_.push_back( begin );

  for ( Point * p = line->EndPoint(); p!=begin; )
  {
    facet_points_.push_back( p );
    line = p->CutLine( line, side, element );
    if ( line==NULL )
    {
      // we have an open path. search backward.
      p = begin;
      line = begin->CommonLine( facet_points_[1] );
      if ( line==NULL )
      {
        throw std::runtime_error( "line disappeared" );
      }
      for ( ;; )
      {
        line = p->CutLine( line, side, element );
        if ( line==NULL )
        {
          break;
        }
        facet_lines_.insert( facet_lines_.begin(), line );
        p = line->OtherPoint( p );
        facet_points_.insert( facet_points_.begin(), p );
      }
      break;
    }
    facet_lines_.push_back( line );
    p = line->OtherPoint( p );
  }

  // remove used lines (no cutting cuts!)
  for ( std::vector<Line*>::iterator i=facet_lines_.begin();
        i!=facet_lines_.end();
        ++i )
  {
    Line* l = *i;
    cut_lines.erase( l );
  }

  closed_ = line != NULL;

  if ( inner and not closed_ )
  {
    ClosedOnEdge( mesh, element, side );
  }
}

bool GEO::CUT::LineSegment::ClosedOnEdge( Mesh & mesh, Element * element, Side * side )
{
  Point * begin = facet_points_.back();
  Point * end   = facet_points_.front();

  Line * begin_line = facet_lines_.back();
  Line * end_line   = facet_lines_.front();

  // If the gap between the end of the point list and the begin of the point
  // list can be closed, this is a closed graph.
  if ( ClosedOnEdge( mesh, element, side, begin, end, begin_line, end_line ) )
  {
    closed_ = true;
    return true;
  }
  else
  {
    return false;
  }
}

bool GEO::CUT::LineSegment::ClosedOnEdge( Mesh & mesh, Element * element, Side * side,
                                          Point * begin, Point * end,
                                          Line * begin_line, Line * end_line )
{
  std::vector<Edge*> begin_edges;
  std::vector<Edge*> end_edges;

  begin->CutEdge( side, begin_line, begin_edges );
  end->CutEdge( side, end_line, end_edges );

  if ( begin_edges.size()==0 or end_edges.size()==0 )
  {
    if ( begin_line==end_line )
    {
      // if this is already one line that is on an edge, and we do not find
      // new edges, we have a touch and not a cut
      return false;
    }

    // If both points belong to the same element side, we can assume a
    // straight line between those points. This is a level set case with two
    // cuts at the same side.
    if ( begin->CutSide( side, end ) != NULL )
    {
      CloseGap( mesh, element, side, end );
      return true;
    }

    throw std::runtime_error( "no edges here!" );
  }

  if ( ( not begin->NodalPoint( side->Nodes() ) and begin_edges.size()>1 ) or
       ( not end  ->NodalPoint( side->Nodes() ) and   end_edges.size()>1 ) )
  {
    throw std::runtime_error( "edge must be unique" );
  }

  for ( std::vector<Edge*>::iterator i=begin_edges.begin(); i!=begin_edges.end(); ++i )
  {
    Edge * begin_edge = *i;

    for ( std::vector<Edge*>::iterator i=end_edges.begin(); i!=end_edges.end(); ++i )
    {
      Edge * end_edge = *i;

      if ( begin_edge==end_edge )
      {
        // Both points are on the same edge. Close the loop.

        std::vector<Point*> line;
        begin_edge->CutPointsBetween( begin, end, line );
        AppendLines( mesh, element, side, line );

        CloseGap( mesh, element, side, end );

        return true;
      }
      else
      {
#if 1
        // look for a common node within the element

        // If there are more than two nodes within the element, we do not find it.

        Point * begin_point = begin_edge->NodeInElement( element, begin );
        Point * end_point = end_edge->NodeInElement( element, end );

        if ( begin_point!=NULL and end_point!=NULL )
        {
          if ( begin_point==end_point )
          {
            std::vector<Point*> begin_line;
            begin_edge->CutPointsBetween( begin, end_point, begin_line );
            begin_line.push_back( end_point );

            AppendLines( mesh, element, side, begin_line );

            std::vector<Point*> end_line;
            end_edge->CutPointsBetween( end_point, end, end_line );

            AppendLines( mesh, element, side, end_line );

            CloseGap( mesh, element, side, end );

            return true;
          }
          else
          {
            Edge * middle_edge = side->FindEdge( begin_point, end_point );
            if ( middle_edge!=NULL )
            {
              std::vector<Point*> begin_line;
              begin_edge->CutPointsBetween( begin, begin_point, begin_line );
              begin_line.push_back( begin_point );

              AppendLines( mesh, element, side, begin_line );

              std::vector<Point*> middle_line;
              middle_edge->CutPointsBetween( begin_point, end_point, middle_line );
              middle_line.push_back( end_point );

              AppendLines( mesh, element, side, middle_line );

              std::vector<Point*> end_line;
              end_edge->CutPointsBetween( end_point, end, end_line );

              AppendLines( mesh, element, side, end_line );

              return true;
            }
          }
        }
#endif
      }
    }
  }
  return false;
}

void GEO::CUT::LineSegment::AppendLines( Mesh & mesh,
                                         Element * element,
                                         Side * side,
                                         const std::vector<Point*> & points )
{
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    CloseGap( mesh, element, side, *i );
    facet_points_.push_back( *i );
  }
}

void GEO::CUT::LineSegment::AppendLinesReversed( Mesh & mesh,
                                                 Element * element,
                                                 Side * side,
                                                 const std::vector<Point*> & points )
{
  for ( std::vector<Point*>::const_reverse_iterator i=points.rbegin(); i!=points.rend(); ++i )
  {
    CloseGap( mesh, element, side, *i );
    facet_points_.push_back( *i );
  }
}

void GEO::CUT::LineSegment::InsertLines( Mesh & mesh,
                                         Element * element,
                                         Side * side,
                                         const std::vector<Point*> & points )
{
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    facet_lines_.insert( facet_lines_.begin(),
                         mesh.NewLine( facet_points_.front(), *i, side, NULL, element ) );
    facet_points_.insert( facet_points_.begin(), *i );
  }
}

void GEO::CUT::LineSegment::InsertLinesReversed( Mesh & mesh,
                                                 Element * element,
                                                 Side * side,
                                                 const std::vector<Point*> & points )
{
  for ( std::vector<Point*>::const_reverse_iterator i=points.rbegin(); i!=points.rend(); ++i )
  {
    facet_lines_.insert( facet_lines_.begin(),
                         mesh.NewLine( facet_points_.front(), *i, side, NULL, element ) );
    facet_points_.insert( facet_points_.begin(), *i );
  }
}

void GEO::CUT::LineSegment::CloseGap( Mesh & mesh, Element * element, Side * side, Point * end )
{
  facet_lines_.push_back( mesh.NewLine( facet_points_.back(), end, side, NULL, element ) );
}

bool GEO::CUT::LineSegment::Combine( Mesh & mesh, Element * element, Side * side, LineSegment & other )
{
  if ( IsClosed() or other.IsClosed() )
  {
    return false;
  }

  if ( ClosedOnEdge( mesh, element, side, EndPoint(), other.BeginPoint(), EndLine(), other.BeginLine() ) )
  {
    const std::vector<Point*> & points = other.Points();
    AppendLines( mesh, element, side, points );
    ClosedOnEdge( mesh, element, side );
    return true;
  }
  else if ( ClosedOnEdge( mesh, element, side, EndPoint(), other.EndPoint(), EndLine(), other.EndLine() ) )
  {
    const std::vector<Point*> & points = other.Points();
    AppendLinesReversed( mesh, element, side, points );
    ClosedOnEdge( mesh, element, side );
    return true;
  }
  else if ( ClosedOnEdge( mesh, element, side, BeginPoint(), other.BeginPoint(), BeginLine(), other.BeginLine() ) )
  {
    const std::vector<Point*> & points = other.Points();
    InsertLines( mesh, element, side, points );
    ClosedOnEdge( mesh, element, side );
    return true;
  }
  else if ( ClosedOnEdge( mesh, element, side, BeginPoint(), other.EndPoint(), BeginLine(), other.EndLine() ) )
  {
    const std::vector<Point*> & points = other.Points();
    InsertLinesReversed( mesh, element, side, points );
    ClosedOnEdge( mesh, element, side );
    return true;
  }

  return false;
}

GEO::CUT::Side* GEO::CUT::LineSegment::OnSide( Element * element )
{
  std::set<Side*> sides( element->Sides().begin(), element->Sides().end() );

  for ( std::vector<Point*>::iterator i=facet_points_.begin(); i!=facet_points_.end(); ++i )
  {
    Point * p = *i;
    p->Intersection( sides );
  }
  if ( sides.size()>1 )
  {
    throw std::runtime_error( "can touch exactly one element side" );
  }
  else if ( sides.size()==1 )
  {
    return *sides.begin();
  }
  return NULL;
}

void GEO::CUT::LineSegmentList::Create( Mesh & mesh, Element * element, Side * side, bool inner )
{
  std::set<Line*> lines;

  const std::vector<Line*> & cut_lines = side->CutLines();
  for ( std::vector<Line*>::const_iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( element ) )
    {
      lines.insert( l );
    }
  }

  Create( mesh, element, side, lines, inner );

  if ( inner )
  {
    for ( unsigned i=0; i<segments_.size(); ++i )
    {
      if ( segments_[i] != Teuchos::null )
      {
        LineSegment & s1 = *segments_[i];
        for ( unsigned j=i+1; j<segments_.size(); ++j )
        {
          if ( segments_[j] != Teuchos::null )
          {
            LineSegment & s2 = *segments_[j];
            if ( s1.Combine( mesh, element, side, s2 ) )
            {
              segments_[j] = Teuchos::null;
            }
          }
        }
      }
    }
  }
}

void GEO::CUT::LineSegmentList::Create( Mesh & mesh, Element * element, Side * side, std::set<Line*> & lines, bool inner )
{
  while ( lines.size() )
  {
    segments_.push_back( Teuchos::rcp( new LineSegment( mesh, element, side, lines, inner ) ) );
  }
}
