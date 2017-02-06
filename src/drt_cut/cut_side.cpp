/*!-----------------------------------------------------------------------------------------------*
 \file cut_side.cpp

 \brief class representing a geometrical side

 <pre>
\level 3
\maintainer Andy Wirtz
 wirtz@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15270
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "cut_side.H"
#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_pointgraph.H"
#include "cut_levelsetside.H"

#include "../drt_lib/drt_globalproblem.H"

#include <string>
#include <stack>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side::Side(int sid, const std::vector<Node*> & nodes,
    const std::vector<Edge*> & edges)
    : sid_(sid),
      nodes_(nodes),
      edges_(edges)
{
  for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end();
      ++i)
  {
    Edge * e = *i;
    e->Register(this);
  }
  selfcutposition_ = Point::undecided;

  if (sid > -1)
  {
    boundingvolume_ = Teuchos::rcp( BoundingBox::Create(*this) );
  }
}

/*----------------------------------------------------------------------------*
      Returns the edge of this side with given begin and end points
 *----------------------------------------------------------------------------*/
GEO::CUT::Edge * GEO::CUT::Side::FindEdge( Point * begin, Point * end )
{
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( e->Matches( begin, end ) )
    {
      return e;
    }
  }
  return NULL;
}
/*----------------------------------------------------------------------------*
 * Calculate the points at which the other side intersects with this considered
 * side
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindCutPoints( Mesh & mesh, Element * element, Side & other,
    int recursion )
{
  bool cut = false;
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    if ( e->FindCutPoints( mesh, element, *this, other, recursion ) )
    {
      cut = true;
    }
  }
  return cut;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindCutLines( Mesh & mesh, Element * element, Side & other )
{
  /* check whether cut lines are already created (still need to create lines new
   * in case we create AmbiguousCutLines!) */
  for ( std::vector<Line*>::iterator i=cut_lines_.begin(); i!=cut_lines_.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( this, &other ) )
    {
      if ( not l->IsCut( element ) )
        dserror("Line is cut by both sides but not by the element, check this "
            "situation as it is not expected!");
        //l->AddElement( element );
      other.AddLine( l );
    }
  }

  // creating cut lines between the given sides for the first time
  PointSet cuts;
  GetCutPoints( element, other, cuts );

  switch ( cuts.size() )
  {
    // ------------------------------------------------------------------------
    // no point --> no line!
    // ------------------------------------------------------------------------
    case 0:
    {
      return false;
    }
    // ------------------------------------------------------------------------
    // just touching point, there shouldn't be any line!
    // ------------------------------------------------------------------------
    /* This is also the default case for element.Dim() == 1. But in this case
     * the construction of cut_lines_ seems unnecessary, since the CutPoints
     * should be sufficient.
     *
     * (see also IMPL::PointGraph1D::AddCutLinesToGraph)     hiermeier 11/16 */
    case 1:
    {
      return false;
    }
    // ------------------------------------------------------------------------
    // The normal case. A straight cut.
    // ------------------------------------------------------------------------
    case 2:
    {
      std::vector<Point*> c;
      c.reserve( 2 );
      c.assign( cuts.begin(), cuts.end() );
      mesh.NewLine( c[0], c[1], this, &other, element );
      return true;
    }
    default:
    {
      return other.FindAmbiguousCutLines( mesh, element, *this, cuts );
    }
  }

  dserror("How did you get here?");
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::AllOnNodes( const PointSet & points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( PointSet::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::GetCutPoints( Element * element, Side & other, PointSet & cuts )
{
  //Get all Cut Points intersecting this side and an edge of the other side
  {
    const std::vector<Edge*> & edges = Edges();
    for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      e->GetCutPoints( element, *this, other, cuts );
    }
  }
  //Get all Cut Points intersecting the other side and an edge of this side
  {
    const std::vector<Edge*> & edges = other.Edges();
    for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      e->GetCutPoints( element, other, *this, cuts );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::AddPoint( Point * cut_point )
{
  cut_points_.insert( cut_point );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::AddLine( Line* cut_line )
{
  if ( std::find( cut_lines_.begin(), cut_lines_.end(), cut_line )==cut_lines_.end() )
  {
    cut_lines_.push_back( cut_line );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Facet * GEO::CUT::Side::FindFacet( const std::vector<Point*> & facet_points )
{
  for ( std::vector<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->Equals( facet_points ) )
    {
      return f;
    }
  }
  return NULL;
}

/*----------------------------------------------------------------------------*/
// Find Cut Lines for two Cut Sides, which have more than two common cut points!
//(This happens if the cutsides are in the same plane )              ager 08/15
/*----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindTouchingCutLines( Mesh & mesh, Element * element,
    Side & side, const PointSet & cut )
{
  // More that two cut points show a touch.
  //
  //1// If all nodes are caught and nothing else, the cut surface has hit this
  // surface exactly. No need to cut anything. However, the surface might be
  // required for integration.
  {
  //find if this side is completely inside the other side
    const std::vector<Node*> & nodes = Nodes();
    if ( cut.size()==nodes.size() and AllOnNodes( cut ) )
    {
      for ( unsigned i=0; i<nodes.size(); ++i )
      {
        unsigned j = ( i+1 ) % nodes.size();
        mesh.NewLine( nodes[i]->point(), nodes[j]->point(), this, &side, element );
      }
      return true;
    }
  }
  {
    //find if the other side is completely inside this side
    const std::vector<Node*> & nodes_o = side.Nodes();
    if ( cut.size()==nodes_o.size() and side.AllOnNodes( cut ) )
    {
      for ( unsigned i=0; i<nodes_o.size(); ++i )
      {
        unsigned j = ( i+1 ) % nodes_o.size();
        mesh.NewLine( nodes_o[i]->point(), nodes_o[j]->point(), this, &side, element );
      }
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 * Find Cut Lines for two Cut Sides specially based on a discretization,
 * which have more than two common cut points!
 *   (This happens if the cut sides are in the same plane or due to numerical
 *   tolerances!)                                                    ager 08/15
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindAmbiguousCutLines( Mesh & mesh, Element * element,
    Side & side, const PointSet & cut )
{
  /* More than two cut points shows a touch.
   *
   * (1) If all nodes are caught and nothing else, the cut surface has hit this
   *     surface exactly. No need to cut anything. However, the surface might be
   *     required for integration. */
  if ( FindTouchingCutLines(mesh,element,side,cut) )
    return true;

  /* (2) Not all nodes are caught but some cut points lie on an edge of the cut sides,
   *     try to connect all points which lie on a line and finally connect all other
   *     missing cut points!!! */
  {
    /* this might not be really high performance, but it is just handling of a special
     * case an therefore will not occur too often
     * --> Speed doesn't matter, but robustness does !!! */
    int created_lines = 0;
    std::vector<int> p1lines;
    std::vector<int> p2lines;
    PointSet c;
 //   std::vector<Point*> c;
    // loop over all edges of the side side
    for ( unsigned ledge = 0; ledge < side.Edges().size(); ++ledge )
    {
      c.clear();
      for ( unsigned lcpoint = 0; lcpoint < cut.size(); ++lcpoint )
      {
        //find all common points of this edge and cut
        for ( unsigned lpoint = 0; lpoint < side.Edges()[ledge]->CutPoints().size(); ++lpoint )
        {
     //     std::cout << "check point " << side.Edges()[ledge]->CutPoints()[lpoint]->Id() << " and " << cut[lcpoint]->Id() << std::endl;
     //     std::cout << "with idx " << lpoint << " and " << lcpoint << std::endl;
          if (side.Edges()[ledge]->CutPoints()[lpoint]->Id() == cut[lcpoint]->Id())
            c.insert(cut[lcpoint]); //
        }
      }
      switch ( c.size() )
      {
        case 0:
        case 1:
          break;
        case 2:
        {
          mesh.NewLine( c[0], c[1], this, &side, element );
          created_lines++;
          p1lines.push_back(c[0]->Id());
          p2lines.push_back(c[1]->Id());
          break;
        }
        default:
        {
          FindAmbiguousCutLinesonEdge( mesh, element, side, c );
        }
      }
    }
    // loop over all edges of the side side
    for ( unsigned ledge = 0; ledge < Edges().size(); ++ledge )
    {
      c.clear();
      for ( unsigned lcpoint = 0; lcpoint < cut.size(); ++lcpoint )
      {
        //find all common points of this edge and cut
        for ( unsigned lpoint = 0; lpoint < Edges()[ledge]->CutPoints().size(); ++lpoint)
        {
          if (Edges()[ledge]->CutPoints()[lpoint]->Id() == cut[lcpoint]->Id())
            {
              c.insert(cut[lcpoint]); //
            }
        }
      }
      switch (c.size())
      {
        case 0:
        case 1:
          break;
        case 2:
        {
          mesh.NewLine( c[0], c[1], this, &side, element );
          created_lines++;
          p1lines.push_back(c[0]->Id());
          p2lines.push_back(c[1]->Id());
          break;
        }
        default:
        {
          FindAmbiguousCutLinesonEdge( mesh, element, side, c );
        }
      }
    }

    /* (3) Now connect all missing points together --> the idea is that all
     * points basically are on one line with just a little deviation to this
     * line, therefore we construct a vector which points in direction of this
     * line and sort all cut points according to this line! Finally we connect
     * these points by lines!
     * (REMARK: If this is not enough to construct lines, we run into an dserror
     * and have to rethink if there are other special cases possible!) */
    if (p1lines.size() < cut.size() -1)
    {
      FindAmbiguousCutLinesonEdge( mesh, element, side, cut );
      return true;
    }
    else //connected already all points :-)!
    {
      return true;
    }
  }
  dserror("FindAmbiguousCutLines failed!!!");
  return false;
}

/*----------------------------------------------------------------------------*
 * Create lines based on cut points arising from an intersection, which should
 * basically be a straight line! Due to numeric precision and finally merging
 * of Points this can happen! The approach is to create a vector, which points
 * in direction of the line and then sort all points along this vector, and
 * connect the points in this sorted order!                          ager 08/15
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindAmbiguousCutLinesonEdge( Mesh & mesh, Element * element,
    Side & side, const PointSet & cut )
{
  double dist2 = 0.0;
  int idxp1 = 0;
  int idxp2 = 0;
  //1// Find Poins with maximum distance in between (Brute Force)
  for (uint i = 0; i < cut.size(); ++i)
  {
    for (uint j = 0; j < cut.size(); ++j)
    {
      if (i!=j)
      {
        double tmpdist2 = 0.0;
        for (uint dim = 0; dim < 3; ++dim)
          tmpdist2 += (cut[i]->X()[dim]-cut[j]->X()[dim])*(cut[i]->X()[dim]-cut[j]->X()[dim]);
        if (dist2 < tmpdist2)
        {
          dist2 = tmpdist2;
          idxp1 = i;
          idxp2 = j;
        }
      }
    }
  }

  // (2) create vector between these two points
  std::vector<double> linevec;
  for (uint dim = 0; dim < 3; ++dim)
    linevec.push_back(cut[idxp1]->X()[dim]-cut[idxp2]->X()[dim]);

  /* (3) evaluate distance from each involved cut point to the first point,
   * projected onto the line vector */
  std::vector<double> projpointdist;
  for (uint cpoint = 0; cpoint < cut.size(); ++cpoint)
  {
    projpointdist.push_back(0.0);
    for (uint dim = 0; dim < 3; ++dim)
      projpointdist[cpoint] += (cut[idxp1]->X()[dim]-cut[cpoint]->X()[dim])*linevec[dim];
  }

  // (4) sort the points accordingly to the projected distance from step //3//

  std::vector<int> sortedcutpoints;
  for (uint cpoint = 0; cpoint < cut.size(); ++cpoint)
    sortedcutpoints.push_back(cpoint);

  bool sorted = false;
  while (!sorted)
  {
    sorted = true;
    for (uint cpoint = 0; cpoint < cut.size()-1; ++cpoint)
    {
      if (projpointdist[sortedcutpoints[cpoint]] > projpointdist[sortedcutpoints[cpoint+1]])
      {
        int tmpidx = sortedcutpoints[cpoint+1];
        sortedcutpoints[cpoint+1] = sortedcutpoints[cpoint];
        sortedcutpoints[cpoint] = tmpidx;
        sorted = false;
      }
    }
  }

  // (5) Now create lines between these points!
  for (uint cpoint = 0; cpoint < cut.size()-1; ++cpoint)
  {
    mesh.NewLine( cut[sortedcutpoints[cpoint]], cut[sortedcutpoints[cpoint+1]],
        this, &side, element );
  }
  return true;
}

void GEO::CUT::Side::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetBoundaryCells( bcells );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeOwnedSideFacets( Mesh & mesh, Element * element,
    plain_facet_set & facets )
{
  // --- facet already created from another element ---------------------------
  // fill the output variable and return
  if ( facets_.size() !=0 )
  {
    std::copy( facets_.begin(), facets_.end(), std::inserter( facets, facets.begin() ) );
    return;
  }

  // create the pointgraph for the element_side
  Teuchos::RCP<IMPL::PointGraph> pg = Teuchos::rcp( IMPL::PointGraph::Create(
      mesh, element, this, IMPL::PointGraph::element_side,
      IMPL::PointGraph::all_lines ) );

  for ( IMPL::PointGraph::facet_iterator i=pg->fbegin(); i!=pg->fend(); ++i )
  {
    const Cycle & points = *i;

    Facet * f = mesh.NewFacet( points(), this, IsCutSide() );
    if ( f == NULL )
      run_time_error( "failed to create owned facet" );
    facets_.push_back( f );
  }

  for ( IMPL::PointGraph::hole_iterator i=pg->hbegin(); i!=pg->hend(); ++i )
  {
    const std::vector<Cycle> & hole = *i;

    // If we have a hole and multiple cuts we have to test which facet the
    // hole belongs to. Not supported now.
    unsigned int facetid = 0;
    if ( facets_.size()!=1 )
    {
      for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
      {
        Facet * facet = * i;
        if ( HoleOfFacet( *facet, hole ) )
        {
          break;
        }
        facetid++;
      }
      if ( facetid == facets_.size() )
      {
        run_time_error( "failed to find the facet of the hole" );
      }
    }

    for ( std::vector<Cycle>::const_iterator i=hole.begin(); i!=hole.end(); ++i )
    {
      const Cycle & points = *i;

      Facet * h = mesh.NewFacet( points(), this, false );
      facets_[facetid]->AddHole( h );
    }
  }

  std::copy( facets_.begin(), facets_.end(), std::inserter( facets, facets.begin() ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeInternalFacets( Mesh & mesh, Element * element,
    plain_facet_set & facets )
{
  // create the pointgraph for the cut_side
  Teuchos::RCP<IMPL::PointGraph> pg = Teuchos::rcp( IMPL::PointGraph::Create(
      mesh, element, this, IMPL::PointGraph::cut_side,
      IMPL::PointGraph::all_lines ) );

  for ( IMPL::PointGraph::facet_iterator i=pg->fbegin(); i!=pg->fend(); ++i )
  {
    const Cycle & points = *i;
    MakeInternalFacets( mesh, element, points, facets );
  }
  for ( IMPL::PointGraph::hole_iterator i=pg->hbegin(); i!=pg->hend(); ++i )
  {
    const std::vector<Cycle> & hole = *i;
    for ( std::vector<Cycle>::const_iterator i=hole.begin(); i!=hole.end(); ++i )
    {
      const Cycle & points = *i;
      MakeInternalFacets( mesh, element, points, facets );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeInternalFacets( Mesh & mesh, Element * element,
    const Cycle & points, plain_facet_set & facets )
{
  /* ignore cycles with all points on one and the same edge
   * ( ToDo seems meaningful for 3-D, check it for the other cases ) */
  // ignore cycles with points outside the current element
  if ( ( not points.IsValid() and element->Dim() == 3 ) or
         not points.IsCut( element ) )
    return;

  Side * s = NULL;

  plain_side_set sides( element->Sides().begin(), element->Sides().end() );
  points.Intersection( sides );

  if ( sides.size()>1 )
  {
    std::stringstream str;
    str << "can touch exactly one element side: "
        << points
        << "found sides:\n";
    std::copy( sides.begin(), sides.end(), std::ostream_iterator<Side*>( str, "\n" ) );
    run_time_error( str.str() );
  }
  else if ( sides.size()==1 )
  {
    s = *sides.begin();
  }
  //else s==NULL :-)

  if ( s!=NULL )
  {
    Facet * f = s->FindFacet( points() );
    if ( f!=NULL )
    {
      f->ExchangeSide( this, true );
      facets.insert( f );
      facets_.push_back( f );
    }
    /* this case means that side pointers show, that the cut_facet lies on
     * the side of the element, but we cannot find this facet??? */
    else
    {
      // we throw an error at this point since this is not allowed -- hiermeier 11/16
      {
        std::cout << "\nOn cut side " << this;
        std::cout << "\nFound element side " << s;
        std::cout << "\nIs cut = " << ( s->IsCut() ? "TRUE" : "FALSE" );

        std::cout << "\n --- Mesh " << &mesh << " ---\n";
        for ( std::map<int, Teuchos::RCP<Node> >::const_iterator cit = mesh.Nodes().begin();
            cit != mesh.Nodes().end(); ++cit )
          std::cout << "Point" << cit->second->point() << "\n";

        std::cout << "\n --- Element " << element << " ---\n";
        for ( std::vector<Point*>::const_iterator cit = element->Points().begin();
              cit != element->Points().end(); ++cit )
          std::cout << "Point " << (*cit) << "\n";

        std::cout << "\n --- Side " << s << " ---\n";
        s->Print();

        std::cout << "\n\n --- PointCycle ---" << std::endl;
        for ( std::vector<Point*>::const_iterator cit = points().begin();
              cit != points().end(); ++cit )
          std::cout << "Point " << (*cit) << "\n";
        std::cout << std::endl;

        run_time_error( "The point cycle was found on the side, but the facet was not found "
            "on the side. This shouldn't happen, since the facet has to be owned by the side "
            "in this case! Properly there is some hanging node in your cut mesh. Take a closer "
            "look." );
      }

      // multiple facets on one cut side within one element
      Facet * f = mesh.NewFacet( points(), this, true );
      facets.insert( f );
      facets_.push_back( f );
  //    dserror("must have matching facet on side");
    }
  }
  else
  {
    // insert new internal facet
    Facet * f = mesh.NewFacet( points(), this, true );
    facets.insert( f );
    facets_.push_back( f );
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnSide( const PointSet & points )
{
  if ( nodes_.size()==points.size() )
  {
    for ( std::vector<Node*>::iterator i=nodes_.begin();
          i!=nodes_.end();
          ++i )
    {
      Node * n = *i;
      if ( points.count( n->point() )==0 )
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnEdge( Point * point )
{
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( point->IsCut( e ) )
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnEdge( Line * line )
{
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( line->OnEdge( e ) )
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HaveCommonNode( Side & side )
{
  const std::vector<Node*> & other_nodes = side.Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * e = *i;
    if ( std::find( other_nodes.begin(), other_nodes.end(), e )!=other_nodes.end() )
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HaveCommonEdge( Side & side )
{
  const std::vector<Edge*> & other_edges = side.Edges();
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( std::find( other_edges.begin(), other_edges.end(), e )!=other_edges.end() )
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Element * GEO::CUT::Side::CommonElement( Side * other )
{
  plain_element_set intersection;
  std::set_intersection( elements_.begin(), elements_.end(),
                         other->elements_.begin(), other->elements_.end(),
                         std::inserter( intersection, intersection.begin() ) );
  switch ( intersection.size() )
  {
  case 0:
    return NULL;
  case 1:
    return *intersection.begin();
  default:
    throw std::runtime_error( "sides with more than one element in common" );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::Print()
{
  std::cout << "[ ";
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    e->Print();
    std::cout << " ; ";
    if ( i+1 != edges_.end() )
      std::cout << "\n  ";
  }
  std::cout << " ]";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned GEO::CUT::Side::UncutFacetNumberPerSide() const
{
  if ( elements_.size() < 1 )
    dserror( "We need at least one registered element to be able to tell the "
        "default facet number of this side." );

  const unsigned ele_dim = elements_[ 0 ]->Dim();
  switch ( ele_dim )
  {
    // 3-D elements: in the uncut case we have one facet per side
    case 3:
      return 1;
    /* 2-D elements: in the uncut case we have four facets per side ( lines
     *               surrounding the element / side ) */
    case 2:
      return 4;
    /* 1-D elements: in the uncut case we have two facets per side ( corner
     *               points of the element / side / edge ) */
    case 1:
      return 2;
    default:
      dserror( "Unsupported parent element dimension! (ele->Dim() = %d )",
         ele_dim );
      exit( EXIT_FAILURE );
  }
  // can never be reached
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::IsCut()
{
  if ( facets_.size() > UncutFacetNumberPerSide() )
  {
    return true;
  }
  if ( facets_[0]->OnCutSide() )
  {
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream & operator<<( std::ostream & stream, GEO::CUT::Side & s )
{
  stream << "side: {";
  const std::vector<GEO::CUT::Node*> & nodes = s.Nodes();
  for ( std::vector<GEO::CUT::Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    GEO::CUT::Node * n = *i;
    n->point()->Print( stream );
    stream << ",";
  }
  stream << "}";
  return stream;
}

/*----------------------------------------------------------------------------*
 *  Gets the selfcutposition and spreads the positional information wirtz 05/13
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::GetSelfCutPosition( Point::PointPosition position )
{
  if ( selfcutposition_ != position )
  {
    selfcutposition_ = position;

    for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
    {
      Edge * e = *i;
      Point::PointPosition ep = e->SelfCutPosition();
      if ( ep==Point::undecided )
      {
        e->SelfCutPosition( position );
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *  Changes the selfcutposition of this cutside and spreads the positional
 *  information                                                    wirtz 07/16
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::ChangeSelfCutPosition( Point::PointPosition position )
{
  if ( selfcutposition_ != position )
  {
    selfcutposition_ = position;

    for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
    {
      Edge * e = *i;
      e->ChangeSelfCutPosition( position );
    }
  }
}

/*----------------------------------------------------------------------------*
 * returns true if the hole is inside the facet                     wirtz 05/13
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HoleOfFacet( Facet & facet, const std::vector<Cycle> & hole )
{

  int intersectioncount = 0;
  bool intersectioninpoint = true;
  std::vector<Point*> facetpoints = facet.Points();
  int facetsize = facetpoints.size();
  std::vector<LINALG::Matrix<3,1> > facetpointslocalcoord;
  facetpointslocalcoord.reserve( facetsize );
  for ( std::vector<Point*>::iterator i=facetpoints.begin(); i!=facetpoints.end(); ++i )
  {
    Point * facetpoint = * i;
    LINALG::Matrix<3,1> pointcoord;
    facetpoint->Coordinates( pointcoord.A() );
    LINALG::Matrix<3,1> pointlocalcoord;
    LocalCoordinates( pointcoord, pointlocalcoord, false );
    facetpointslocalcoord.push_back( pointlocalcoord );
  }
  LINALG::Matrix<3,1> holepointcoord;
  LINALG::Matrix<3,1> holepointlocalcoord;
  hole[0]()[0]->Coordinates( holepointcoord.A() );
  LocalCoordinates( holepointcoord, holepointlocalcoord, false );
  double epsilon = 0;
  while ( intersectioninpoint )
  {
    intersectioninpoint = false;
    for ( std::vector<LINALG::Matrix<3,1> >::iterator
        i=facetpointslocalcoord.begin(); i!=facetpointslocalcoord.end(); ++i )
    {
      LINALG::Matrix<3,1> facetpoint1 = * i;
      LINALG::Matrix<3,1> facetpoint2;
      if ( i+1 != facetpointslocalcoord.end() )
      {
        facetpoint2 = * (i+1);
      }
      else
      {
        facetpoint2 = * (i+1-facetsize);
      }
      double A = facetpoint1(0) - facetpoint2(0);
      double B = facetpoint1(1) - facetpoint2(1);
      double C = facetpoint1(0) + facetpoint2(0) - 2*holepointlocalcoord(0) - 2;
      double D = facetpoint1(1) + facetpoint2(1) - 2*holepointlocalcoord(1) - epsilon;
      double N = 2*B - epsilon*A;
      if ( abs(N) > REFERENCETOL )
      {
        double eta = ( B*C - A*D )/N;
        double xsi = ( 2*D - epsilon*C )/N;
        if ( eta < 1 and eta > -1 and xsi < 1 and xsi > -1)
        {
          intersectioncount++;
          double xlocalcoord = holepointlocalcoord(0) + 1 + eta;
          double ylocalcoord = (2*holepointlocalcoord(1) + epsilon + epsilon*eta)/2;
          if ( (abs(xlocalcoord - facetpoint1(0)) < REFERENCETOL and
                abs(ylocalcoord - facetpoint1(1)) < REFERENCETOL ) or
               (abs(xlocalcoord - facetpoint2(0)) < REFERENCETOL and
                abs(ylocalcoord - facetpoint2(1)) < REFERENCETOL ) )
          {
            intersectioninpoint = true;
            intersectioncount = 0;
            epsilon += REFERENCETOL;
            break;
          }
        }
      }
    }
  }
  if ( intersectioncount%2 == 0 )
  {
    return false;
  }
  else
  {
    return true;
  }

}

/*--------------------------------------------------------------------------*
 * replace the Node "nod" with the new node "replwith"         sudhakar 10/13
 * Modify also the edge informations correspondingly
 *--------------------------------------------------------------------------*/
void GEO::CUT::Side::replaceNodes( Node* nod, Node* replwith )
{
  bool replaced = false;

  for( unsigned i=0; i < nodes_.size(); i++ )
  {
    Node* orig = nodes_[i];

    if( orig->Id() == nod->Id() )
    {
      nodes_[i] = replwith;
      replaced = true;
    }
  }

  if( not replaced )
    return;

  // also modify the corresponding edge information
  for( std::vector<Edge*>::iterator it = edges_.begin(); it != edges_.end(); it++ )
  {
    Edge* ed = *it;
    ed->replaceNode( nod, replwith );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType sidetype,
           unsigned numNodesSide,
           unsigned dim>
bool GEO::CUT::ConcreteSide<probdim,sidetype,numNodesSide,dim>::IsCloserSide(
    const LINALG::Matrix<probdim,1>& startpoint_xyz, GEO::CUT::Side* other, bool& is_closer)
{
  /* shoot a ray starting from the startpoint through the midpoint of this side
   * and find an intersection point with the other side */
  LINALG::Matrix<probdim,1> ray_point_xyz(true);
  // as second point on the ray we define the midpoint of this side
//  this->SideCenter(ray_point_xyz);


  /* choose a point inside the side such that the angle between the normal of
   * the other side and the vector between start-point and ray-point is next to
   * 0 or 180 degree to avoid orthogonal ray and normal vector therefore, for
   * almost parallel sides and a start-point next to the common point of the sides
   * the side-center might lead to a ray which is almost parallel to the other side */

  LINALG::Matrix< dim, numNodesSide  > corner_coords_rst( true );
  this->LocalCornerCoordinates( corner_coords_rst.A() );

  /* shrink/perturb the local coordinates around the center point with a given
   * tolerance to obtain points which are next to the corner points however slightly
   * inside */
  LINALG::Matrix<dim,1> rst_center = DRT::UTILS::getLocalCenterPosition<dim>(this->Shape());

  //-----------------------------
  // get perturbed coordinates
  //-----------------------------
  LINALG::Matrix<dim,numNodesSide> inner_corner_coords_rst( true );

  // 1. transform such that coordinates center is located in the element center
  for( unsigned i=0; i<numNodesSide; ++i )
  {
    for( unsigned j=0; j<dim; ++j )
      inner_corner_coords_rst(j,i) = corner_coords_rst(j,i) - rst_center(j);
  }

  /* 2. shrink the side coordinates (hard-coded tolerance w.r.t. local coordinates of
   *    element, this should be fine) */
  const double TOL = 1e-003;
  const double scalefac = 1.0-TOL;
  inner_corner_coords_rst.Scale(scalefac);

  // 3. transform the element back
  for( unsigned i=0; i<numNodesSide; ++i )
  {
    for( unsigned j=0; j<dim; ++j )
      inner_corner_coords_rst(j,i) += rst_center(j);
  }

  //-----------------------------
  /* choose the center point point or a perturbed inner corner point
   * such that the ray between startpoint and this point is as perpendicular
   * as possible to guarantee well-conditioned systems for finding ray-cut points */
  //-----------------------------

  LINALG::Matrix<probdim,1> xyz(true);
  LINALG::Matrix<probdim,1> ray_dir(true);

  // get normal of other side at its center
  LINALG::Matrix<probdim,1> t1,t2,n(true);
  other->BasisAtCenter(t1,t2,n);


  // start with center point
  this->SideCenter(xyz); // as second point on the ray we define the midpoint of this side
  ray_dir.Update(1.0, xyz, -1.0, startpoint_xyz, 0.0);
  ray_point_xyz.Update(1.0, xyz, 0.0);

  double cosine = ray_dir.Dot(n) / ray_dir.Norm2(); // n is normalized


  // loop corner nodes
  for( unsigned i= 0; i< numNodesSide; i++)
  {
    // get ray vector (endpoint-startpoint)
    LINALG::Matrix<dim,1> rst_inner_corner(&inner_corner_coords_rst(0,i),true);
    PointAt( rst_inner_corner, xyz);
    ray_dir.Update(1.0, xyz, -1.0, startpoint_xyz, 0.0);

    double cosine_tmp = ray_dir.Dot(n) / ray_dir.Norm2(); // n is normalized

    /* maximize the absolute value of the cosine to choose the ray which
     * is as perpendicular to the side as possible */
    if(fabs(cosine_tmp) > fabs(cosine))
    {
      ray_point_xyz.Update(1.0, xyz, 0.0);
      cosine = cosine_tmp;
    }
  }
  //-----------------------------
  /* shoot the ray and find a cutpoint with the other side's plane or curved
   * surface space */
  //-----------------------------
  LINALG::Matrix<dim,1> rs(true);
  double line_xi = 0.0;

  bool cut_found = other->RayCut( startpoint_xyz, ray_point_xyz, rs, line_xi);

  if(!cut_found) return false;
  else
  {
    // The main decision if the side lies closer to the start-point than the other side

//    std::cout << "line_xi " << line_xi << std::endl;

    if(line_xi > 1.0+REFERENCETOL)
    {
      // the first side is closer to the start point than the second side
      is_closer = true;
      return true;
    }
    else if(fabs(line_xi - 1.0) <= REFERENCETOL)
    {
      // the found intersection point on the other is the midpoint of the fist side
      // in that case both sides lie within one plane
      // this case is catched in SameNormal afterwards

      // std::cout << "line_xi " << line_xi << std::endl;
      // std::cout << "check if both sides lie in one plane " << std::endl;

      return false;
    }
    else if(line_xi < 1.0-REFERENCETOL and line_xi > -1.0+REFERENCETOL)
    {
      // the most safe check (accept the other side as the nearest side)
      is_closer = false;
      return true;
    }
    else if(fabs(line_xi + 1.0) <= REFERENCETOL )
    {
      /* the intersection point is the same as the start-point of the ray
       * the other side contains the start-point and the cut-point shared
       * with the original side in that case the line between the start-point
       * and the cut-point lies in the second side then the side is orthogonal
       * to the line this side should be removed in */
      std::cout << "line_xi " << line_xi << std::endl;
      std::cout << "start-point: " << startpoint_xyz << std::endl;
      std::cout << "side orthogonal ? " << std::endl; other->Print();

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");

      return false;
    }
    else if(line_xi < -1.0-REFERENCETOL)
    {
      // both sides lead to the same result
      is_closer = true; // false would be also okay
      return true;
    }
    else
    {
      // undermined range of local coordinates!

      std::cout << "line_xi " << line_xi << std::endl;
      std::cout << "cut point found, but the local line coordinates along the "
          "ray-tracing line lies in undefined region" << std::endl;

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");
    }

    return false;
  }

  return false; // return not successful
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType sidetype,
           unsigned numNodesSide,
           unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim,sidetype,numNodesSide,dim>::WithinSide(
    const LINALG::Matrix<probdim, 1> & xyz,
    LINALG::Matrix<dim, 1> & rs, double & dist)
{
  Teuchos::RCP<Position> pos =
      PositionFactory::BuildPosition<probdim,sidetype>( *this, xyz );
  bool success = pos->IsGivenPointWithinElement();
  if ( not success )
  {
//      throw std::runtime_error( "ComputeDistance w.r.t side not successful" );
    rs = 0;
    dist = 9999;     // set large value
    return false;
  }
  pos->LocalCoordinates(rs);
  dist = pos->Distance();

  return pos->WithinLimits( false );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType sidetype,
           unsigned numNodesSide,
           unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim,sidetype,numNodesSide,dim>::LocalCoordinates(
    const LINALG::Matrix<probdim, 1> & xyz,
    LINALG::Matrix<probdim, 1> & rsd, bool allow_dist)
{
  Teuchos::RCP<Position> pos = PositionFactory::BuildPosition<probdim,sidetype>( *this, xyz );
  bool success = pos->Compute();
  LINALG::Matrix<dim,1> rs( true );
  pos->LocalCoordinates( rs );
  // copy the position
  std::copy( rs.A(), rs.A() + dim, &rsd( 0 ) );
  // copy the distance
  switch ( dim )
  {
    /* Side dimension is by a factor of one smaller than the problem dimension,
     * i.e. surface element in 3-D or line element in 2-D. */
    case probdim-1 :
    {
      rsd( dim ) = pos->Distance();
      break;
    }
    /* Side dimension is by a factor of two smaller than the problem dimension,
     * i.e. line element in 3-D. */
    case probdim-2 :
    {
      dserror( "Unsupported dim / probdim combination. I think this can't happen "
          "for side elements. If I'm wrong, just ask me. -- hiermeier" );
      exit( EXIT_FAILURE );
    }
    /* For a 1-D and a 2-D problem, the side dimension is equal to the problem
     * dimension */
    default :
    {
      /* do nothing, there is no distance in the 2-D non-embedded standard case */
      break;
    }
  }

  return success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType sidetype,
           unsigned numNodesSide,
           unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim,sidetype,numNodesSide,dim>::RayCut(
    const LINALG::Matrix<probdim, 1> & p1_xyz,
    const LINALG::Matrix<probdim, 1> & p2_xyz,
    LINALG::Matrix<dim, 1> & rs,
    double & line_xi)
{
  LINALG::Matrix<probdim,numNodesSide> xyze_surface(true);
  this->Coordinates(xyze_surface);

  LINALG::Matrix<probdim,2> xyze_line(true);
  for (unsigned i=0; i<probdim; ++i)
  {
    xyze_line(i,0) = p1_xyz(i);
    xyze_line(i,1) = p2_xyz(i);
  }
  /* The dim+1 entry corresponds to the parameter space coordinate of the
   * 1-D line element. */
  LINALG::Matrix<dim+1,1> xsi(true);

  // do not check for within-limits during the Newton-scheme, since the cut-point is
  // allowed to be not within the side and line
  bool checklimits = false;
  GEO::CUT::KERNEL::ComputeIntersection<probdim,DRT::Element::line2, sidetype>
      ci( xsi, checklimits );
//  GEO::CUT::KERNEL::DebugComputeIntersection<probdim,DRT::Element::line2, sidetype>
//      ci( xsi, checklimits );

  // successful line-side intersection
  if ( ci( xyze_surface, xyze_line ) )
  {
    std::copy(xsi.A(),xsi.A()+dim,&rs(0));
    line_xi = xsi(dim);

    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::SideFactory::CreateSide(
    DRT::Element::DiscretizationType sidetype,
    int sid,
    const std::vector<Node*> & nodes,
    const std::vector<Edge*> & edges) const
{
  Side * s = NULL;
  const int probdim = DRT::Problem::Instance()->NDim();
  switch (sidetype)
  {
    case DRT::Element::line2 :
    {
      s = CreateConcreteSide<DRT::Element::line2>( sid, nodes, edges, probdim );
      break;
    }
    case DRT::Element::tri3 :
      s = CreateConcreteSide<DRT::Element::tri3>( sid, nodes, edges, probdim );
      break;
    case DRT::Element::quad4 :
      s = CreateConcreteSide<DRT::Element::quad4>( sid, nodes, edges, probdim );
      break;
    default:
    {
      dserror("Unsupported side type! ( %d | %s )",
          sidetype,DRT::DistypeToString(sidetype).c_str());
      break;
    }
  }
  return s;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Side::CreateLevelSetSide( const int& sid )
{
  Side * lvs_side_ptr = NULL;

  const int probdim = DRT::Problem::Instance()->NDim();
  switch (probdim)
  {
    case 2:
      lvs_side_ptr = new LevelSetSide<2>(sid);
      break;
    case 3:
      lvs_side_ptr = new LevelSetSide<3>(sid);
      break;
    default:
      dserror("Unsupported problem dimension! (probdim=%d)",
          probdim);
      break;
  }
  return lvs_side_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Side::Create(
    const DRT::Element::DiscretizationType & sidetype,
    const int & sid,
    const std::vector<Node*> & nodes,
    const std::vector<Edge*> & edges)
{
  SideFactory factory;
  return factory.CreateSide(sidetype,sid,nodes,edges);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::Side::Create(
      const unsigned & shardskey,
      const int & sid,
      const std::vector<Node*> & nodes,
      const std::vector<Edge*> & edges)
{
  return Create(DRT::ShardsKeyToDisType(shardskey),sid,nodes,edges);
}


template class GEO::CUT::ConcreteSide<2,DRT::Element::line2>;
template class GEO::CUT::ConcreteSide<3,DRT::Element::line2>;
template class GEO::CUT::ConcreteSide<3,DRT::Element::quad4>;
template class GEO::CUT::ConcreteSide<3,DRT::Element::tri3>;
