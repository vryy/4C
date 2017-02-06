/*!-----------------------------------------------------------------------------------------------*
 \file cut_edge.cpp

 \brief class representing a geometrical edge

 <pre>
\level 3
\maintainer Andy Wirtz
 wirtz@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15270
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_levelsetside.H"

#include "../drt_lib/drt_globalproblem.H"

#include <string>
#include <stack>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::Edge::Create(
    DRT::Element::DiscretizationType edgetype,
    const std::vector<Node*> & nodes)
{
  EdgeFactory factory;
  return factory.CreateEdge( edgetype, nodes );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::Edge::Create(
    unsigned shardskey,
    const std::vector<Node*> & nodes)
{
  return Edge::Create( DRT::ShardsKeyToDisType( shardskey ), nodes );
}

/*-----------------------------------------------------------------------------*
 *  Find points at which this edge which is in "side" cuts the "other"
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Edge::FindCutPoints( Mesh & mesh,
                                    Element * element,
                                    Side & side,
                                    Side & other,
                                    int recursion )
{
  bool cut = false;
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( &other ) ) //cut points may already obtained when some other side was considered
    {
      cut = true;
      p->AddElement( element );
    }
  }
#if 1
  if ( cut and recursion > 0 )
  {
    return true;
  }
#endif

  // test for the cut of edge and side first time in this mesh

  PointSet cut_points;
  other.Cut( mesh, *this, cut_points );
  for ( PointSet::iterator i=cut_points.begin();
        i!=cut_points.end();
        ++i )
  {
    Point * p = *i;

    p->AddEdge( this );
    p->AddElement( element );

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->AddSide( &side );
    p->AddSide( &other );
    AddPoint( p );
  }

  return cut_points.size()>0;
}

/*----------------------------------------------------------------------------*
 * Cut points falling on this edge that are common to the two given sides are
 * extracted
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::GetCutPoints( Element * element,
                                   Side & side,
                                   Side & other,
                                   PointSet & cuts )
{
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( &other ) and p->IsCut( element ) )
    {
      cuts.insert( p );
    }
  }
}

/*----------------------------------------------------------------------------*
 * Cut points falling on this edge that are common to the given edge are
 * extracted
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::GetCutPoints( Edge * other, PointSet & cuts )
{
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( other ) )
    {
      cuts.insert( p );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::AddPoint( Point* cut_point )
{
  // make sure the position of the point on this edge is known
  cut_point->t( this );

#if 0
  std::vector<Point*>::iterator j = std::lower_bound( cut_points_.begin(), cut_points_.end(), cut_point, PointPositionLess( this ) );
  if ( j==cut_points_.end() or *j!=cut_point )
  {
    cut_points_.push_back( cut_point );
    std::sort( cut_points_.begin(), cut_points_.end(), PointPositionLess( this ) );
  }
#else
  cut_points_.insert( cut_point );
#endif

#if 0
#ifdef DEBUGCUTLIBRARY
  PointSet cp;
  std::copy( cut_points_.begin(), cut_points_.end(), std::inserter( cp, cp.begin() ) );
  if ( cut_points_.size() != cp.size() )
    throw std::runtime_error( "broken cut points" );
#endif
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPoint( Node* edge_start, Node* edge_end,
    std::vector<Point*> & edge_points )
{
  Point * bp = BeginNode()->point();
  Point * ep = EndNode()  ->point();
  if ( bp==edge_start->point() and ep==edge_end->point() )
  {
    edge_points.clear();
    edge_points.assign( cut_points_.begin(), cut_points_.end() );
  }
  else if ( ep==edge_start->point() and bp==edge_end->point() )
  {
    edge_points.clear();
    edge_points.assign( cut_points_.rbegin(), cut_points_.rend() );
  }
  else
  {
    throw std::runtime_error( "not a valid edge" );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPoints( Side * side, PointSet & cut_points )
{
  IMPL::SideCutFilter filter( side );
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    plain_line_set cut_lines;
    p->CutLines( filter, cut_lines );
    if ( cut_lines.size()>0 )
    {
      cut_points.insert( p );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsBetween( Point* begin, Point* end,
    std::vector<Point*> & line )
{
//   PointPositionLess::iterator bi = cut_points_.find( begin );
//   PointPositionLess::iterator ei = cut_points_.find( end );
  PointPositionSet::iterator bi = std::lower_bound( cut_points_.begin(),
      cut_points_.end(), begin, PointPositionLess( this ) );
  PointPositionSet::iterator ei = std::lower_bound( cut_points_.begin(),
      cut_points_.end(), end  , PointPositionLess( this ) );

  if ( *bi != begin )
    bi = cut_points_.end();

  if ( *ei != end )
    ei = cut_points_.end();

  double bt = begin->t( this );
  double et = end->t( this );

  if ( bt < et )
  {
    if ( bi!=cut_points_.end() )
    {
      ++bi;
    }
    //std::copy( bi, ei, std::back_inserter( line ) );
    line.insert( line.end(), bi, ei );
  }
  else if ( bt > et )
  {
    if ( ei!=cut_points_.end() )
    {
      ++ei;
    }
    //std::copy( ei, bi, std::back_inserter( line ) );
    line.insert( line.end(), ei, bi );
  }
  else
  {
    if ( begin!=end )
    {
      throw std::runtime_error( "different points at the same place" );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsIncluding( Point* begin, Point* end,
    std::vector<Point*> & line )
{
  PointPositionSet::iterator bi = std::lower_bound( cut_points_.begin(),
      cut_points_.end(), begin, PointPositionLess( this ) );
  PointPositionSet::iterator ei = std::lower_bound( cut_points_.begin(),
      cut_points_.end(), end  , PointPositionLess( this ) );

  if ( *bi != begin )
    throw std::runtime_error( "begin point not on edge" );

  if ( *ei != end )
    throw std::runtime_error( "end point not on edge" );

  double bt = begin->t( this );
  double et = end->t( this );

  if ( bt < et )
  {
    ++ei;
    //std::copy( bi, ei, std::back_inserter( line ) );
    line.insert( line.end(), bi, ei );
  }
  else if ( bt > et )
  {
    ++bi;
    //std::copy( ei, bi, std::inserter( line, line.begin() ) );
    line.insert( line.begin(), ei, bi );
  }
  else
  {
    if ( begin!=end )
    {
      throw std::runtime_error( "different points at the same place" );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsInside( Element * element, std::vector<Point*> & line )
{
  Point * first = NULL;
  Point * last = NULL;
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( element ) )
    {
      if ( first == NULL )
      {
        first = last = p;
      }
      else
      {
        last = p;
      }
    }
  }
  if ( first!=NULL and first!=last )
  {
    CutPointsIncluding( first, last, line );
#if 0
    // Rectify numerical problems that occationally occur. There might be
    // middle points that do not know their position on an element side.
    if ( line.size() > 2 )
    {
      plain_side_set common;
      first->CommonSide( last, common );
      if ( common.size() > 0 )
      {
        for ( std::vector<Point*>::iterator i=line.begin(); i!=line.end(); ++i )
        {
          Point * p = *i;
          for ( plain_side_set::iterator i=common.begin(); i!=common.end(); ++i )
          {
            Side * s = *i;
            p->AddSide( s );
          }
        }
      }
    }
#endif
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Edge::IsCut( Side * side )
{
  // cutpoints contains end-points and internal cut-points
  for ( PointPositionSet::iterator i=cut_points_.begin(); i!=cut_points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( side ) )
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::Edge::NodeInElement( Element * element, Point * other )
{
  Point * p = BeginNode()->point();
  if ( p!=other and p->IsCut( element ) )
  {
    return p;
  }
  p = EndNode()->point();
  if ( p!=other and p->IsCut( element ) )
  {
    return p;
  }
  return NULL;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::RectifyCutNumerics()
{
  if ( cut_points_.size() > 2 )
  {
    // Rectify numerical problems that occasionally occur. There might be middle
    // points that do not know their position on an element side.
    //
    // This assumes linear sides. There might be a problem with quad4
    // sides. Those are actually not supported.

    std::map<Side*, PointPositionSet::iterator> sidecuts;

    for ( PointPositionSet::iterator pi=cut_points_.begin(); pi!=cut_points_.end(); ++pi )
    {
      Point * p = *pi;
      const plain_side_set & cutsides = p->CutSides();
      for ( plain_side_set::const_iterator i=cutsides.begin(); i!=cutsides.end(); ++i )
      {
        Side * s = *i;
        std::map<Side*, PointPositionSet::iterator>::iterator j = sidecuts.find( s );
        if ( j!=sidecuts.end() )
        {
          //if ( std::distance( j->second, pi ) > 1 )
          {
            PointPositionSet::iterator next = j->second;
            next++;

            for ( PointPositionSet::const_iterator i=next; i!=pi; ++i )
            {
              Point * p = *i;
              p->AddSide( s );
            }
          }
          j->second = pi;
        }
        else
        {
          sidecuts[s] = pi;
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*
 *  Gives this edge a selfcutposition and spreads the positional
 *  information                                                 wirtz 05/13
 *------------------------------------------------------------------------*/
void GEO::CUT::Edge::SelfCutPosition( Point::PointPosition pos )
{

#ifdef DEBUGCUTLIBRARY
  if( (selfcutposition_ == Point::inside and pos == Point::outside) or
      (selfcutposition_ == Point::outside and pos == Point::inside) )
  {
    dserror("Are you sure that you want to change the edge-position from inside to outside or vice versa?");
  }
#endif

  if( selfcutposition_ == Point::undecided)
  {

    if ( selfcutposition_ != pos )
    {
      selfcutposition_ = pos;
      if ( pos==Point::outside or pos==Point::inside )
      {
        for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
        {
          Node * n = *i;
          if ( n->SelfCutPosition()==Point::undecided )
          {
            n->SelfCutPosition( pos );
          }
        }
        for ( plain_side_set::iterator i=sides_.begin(); i!=sides_.end(); ++i )
        {
          Side * s = *i;
          s->GetSelfCutPosition( pos );
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*
 *  Changes the selfcutposition of this edge and spreads the positional
 *  information                                                 wirtz 07/16
 *------------------------------------------------------------------------*/
void GEO::CUT::Edge::ChangeSelfCutPosition( Point::PointPosition pos )
{

  if ( selfcutposition_ != pos )
  {
    selfcutposition_ = pos;
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node * n = *i;
      n->ChangeSelfCutPosition( pos );
    }
    for ( plain_side_set::iterator i=sides_.begin(); i!=sides_.end(); ++i )
    {
      Side * s = *i;
      s->ChangeSelfCutPosition( pos );
    }
  }

}

/*------------------------------------------------------------------------*
 *  Replaces the node "nod" of the edge with given node "replwith"
 *                                                              sudhakar 09/13
 *------------------------------------------------------------------------*/
void GEO::CUT::Edge::replaceNode( Node* nod, Node* replwith )
{
  for( unsigned i=0; i < nodes_.size(); i++ )
  {
    Node* orig = nodes_[i];

    if( orig->Id() == nod->Id() )
    {
      nodes_[i] = replwith;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probDim,
           DRT::Element::DiscretizationType edgeType,
           unsigned dimEdge ,
           unsigned numNodesEdge >
bool GEO::CUT::ConcreteEdge<probDim,edgeType,dimEdge,numNodesEdge>::Cut(
    Mesh & mesh, Side & side, PointSet & cuts )
{
  Teuchos::RCP<GEO::CUT::IntersectionBase> inter_ptr =
      IntersectionPtr(side.Shape());

  inter_ptr->Init(&mesh,this,&side,false,false,false);
  return inter_ptr->Intersect( cuts );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probDim,
           DRT::Element::DiscretizationType edgeType,
           unsigned dimEdge ,
           unsigned numNodesEdge >
bool GEO::CUT::ConcreteEdge<probDim,edgeType,dimEdge,numNodesEdge>::ComputeCut(
    Mesh * mesh,
    Edge * other,
    PointSet * cut_points,
    double & tolerance )
{
  // nodal edge coordinates
  LINALG::Matrix<probDim,numNodesEdge> xyze_this;
  Coordinates( xyze_this.A() );

  // nodal pseudo side coordinates
  /* Note: We insert the this edge as side. In this way we can use the
   * standard intersection object. */
  Epetra_SerialDenseMatrix xyze_other;
  xyze_other.Shape( other->ProbDim(), other->NumNodes() );
  other->Coordinates( xyze_other );

  Teuchos::RCP<GEO::CUT::IntersectionBase> inter_ptr =
      IntersectionPtr( other->Shape() );
  inter_ptr->Init( xyze_other, xyze_this, false, false, false );

  // perform the actual edge - edge intersection
  const enum IntersectionStatus istatus =
      inter_ptr->ComputeEdgeSideIntersection( tolerance );

  switch ( istatus )
  {
    case intersect_single_cut_point:
    {
      double * x_ptr = inter_ptr->FinalPoint();
      const double & pos = inter_ptr->LocalCoordinates()[ 0 ];

      Point * p = NULL;
      p = Point::NewPoint( *mesh, x_ptr, pos, other, NULL, tolerance );
      p->AddEdge( this );

      if ( cut_points )
        cut_points->insert( p );

      break;
    }
    case intersect_multiple_cut_points:
    {
      std::vector<LINALG::Matrix<probDim,1> > xyz_cuts;
      std::vector<LINALG::Matrix<dimEdge,1> > r_cuts;
      inter_ptr->FinalPoints( xyz_cuts );
      inter_ptr->LocalSideCoordinates( r_cuts );

      dsassert( xyz_cuts.size() == r_cuts.size(), "Size mismatch!" );

      Point * p = NULL;
      for ( unsigned i=0; i<xyz_cuts.size(); ++i )
      {
        p = Point::NewPoint( *mesh, xyz_cuts[ i ].A(), r_cuts[ i ]( 0 ), other, NULL, tolerance );
        p->AddEdge( this );

       if ( cut_points )
          cut_points->insert( p );
      }

      break;
    }
    case intersect_no_cut_point:
    {
      /* do nothing */
      break;
    }
    default:
      dserror( "Unsupported intersection status! ( status = %d )", istatus );
      exit(EXIT_FAILURE);
  }
  return ( static_cast<bool>( istatus ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probDim,
           DRT::Element::DiscretizationType edgeType,
           unsigned dimEdge ,
           unsigned numNodesEdge >
Teuchos::RCP<GEO::CUT::IntersectionBase> GEO::CUT::ConcreteEdge<probDim,edgeType,
    dimEdge,numNodesEdge>::IntersectionPtr(
    const DRT::Element::DiscretizationType & sidetype) const
{
  return GEO::CUT::IntersectionBase::Create( edgeType, sidetype );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::EdgeFactory::CreateEdge(
        const DRT::Element::DiscretizationType & edgetype,
        const std::vector<Node*> & nodes) const
{
  Teuchos::RCP<Edge> cedge_ptr = Teuchos::null;
  const int probdim = DRT::Problem::Instance()->NDim();
  switch (edgetype)
  {
    case DRT::Element::line2:
    {
      cedge_ptr = Teuchos::rcp( CreateConcreteEdge<DRT::Element::line2>( nodes, probdim ) );
      break;
    }
    default:
    {
      dserror("Unsupported edge type! ( %d | %s )",
          edgetype, DRT::DistypeToString( edgetype).c_str() );
      break;
    }
  }
  return cedge_ptr;
}

template class GEO::CUT::ConcreteEdge<2,DRT::Element::line2>;
template class GEO::CUT::ConcreteEdge<3,DRT::Element::line2>;
