#include "cut_intersection.H"
#include "cut_position.H"
#include "cut_tetmesh.H"
#include "cut_options.H"
#include "cut_integrationcellcreator.H"
#include "cut_volumecellgenerator.H"
#include "cut_facetgraph.H"

#include <string>
#include <stack>

/*--------------------------------------------------------------------*
 *            cut this element with given cut_side
 *--------------------------------------------------------------------*/
bool GEO::CUT::Element::Cut( Mesh & mesh, Side & side, int recursion )
{
  bool cut = false;

  // find nodal points inside the element
  const std::vector<Node*> side_nodes = side.Nodes();

  for ( std::vector<Node*>::const_iterator i=side_nodes.begin(); i!=side_nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();

    if ( not p->IsCut( this ) ) // point does not know if this element is a cut_element_
    {
      if ( PointInside( p ) ) // if point is inside the element
      {
        p->AddElement( this ); // add element to cut_element_-list of this point
        cut = true;
      }
    }
    else // point inside this element, already determined by another side
    {
      cut = true;
    }
  }

  // all the other cut points lie on sides of the element (s is an element side, side is the cutter side)
  const std::vector<Side*> & sides = Sides();
  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    if ( FindCutPoints( mesh, *s, side, recursion ) )
    {
      cut = true;
    }
  }

  // insert this side into cut_faces_
  if ( cut )
  {
    cut_faces_.insert( &side );
    return true;
  }
  else
  {
    return false;
  }
}

void GEO::CUT::Element::MakeCutLines( Mesh & mesh, Creator & creator )
{
  for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
  {
    Side & side = **i;

    bool cut = false;

    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
      if ( FindCutLines( mesh, *s, side ) )
      {
        cut = true;
      }
    }

    // find lines inside the element
    const std::vector<Edge*> & side_edges = side.Edges();
    for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
    {
      Edge * e = *i;
      std::vector<Point*> line;
      e->CutPointsInside( this, line );
      mesh.NewLinesBetween( line, &side, NULL, this );
    }

    if ( cut )
    {
      // create any remaining cut lines
      //side.CreateMissingLines( creator, this );
    }
  }
}

bool GEO::CUT::Element::FindCutPoints( Mesh & mesh, Side & side, Side & other, int recursion )
{
  bool cut = side.FindCutPoints( mesh, this, other, recursion );
  bool reverse_cut = other.FindCutPoints( mesh, this, side, recursion );
  return cut or reverse_cut;
}

bool GEO::CUT::Element::FindCutLines( Mesh & mesh, Side & side, Side & other )
{
  bool cut = side.FindCutLines( mesh, this, other );
  bool reverse_cut = other.FindCutLines( mesh, this, side );
  return cut or reverse_cut;
}

void GEO::CUT::Element::MakeFacets( Mesh & mesh )
{
  if ( facets_.size()==0 )
  {
    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      side.MakeOwnedSideFacets( mesh, this, facets_ );
    }
    for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side & cut_side = **i;
      cut_side.MakeInternalFacets( mesh, this, facets_ );
    }
  }
}

void GEO::CUT::Element::FindNodePositions()
{
  LINALG::Matrix<3,1> xyz;
  LINALG::Matrix<3,1> rst;

  const std::vector<Node*> & nodes = Nodes();

  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();

    if ( pos==Point::undecided )
    {
      bool done = false;
      const plain_facet_set & facets = p->Facets();

      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;

        double smallest_dist = 0.0;

        for ( plain_side_set::const_iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
        {
          Side * s = *i;

          // Only take a side that belongs to one of this point's facets and
          // shares a cut edge with this point. If there are multiple cut
          // sides within the element (facets), only the close one will always
          // give the right direction.
          if ( f->IsCutSide( s ) and p->CommonCutEdge( s )!=NULL )
          {
            if ( p->IsCut( s ) )
            {
              p->Position( Point::oncutsurface );
            }
            else
            {
              p->Coordinates( xyz.A() );
              // the local coordinates here are wrong! (local coordinates of a point lying not within the side!?)
              // but the third component gives the smallest distance of the point to the lines (line-segments) of this side
              // (not to the inner of this side!!!, needs only linear operations and no newton)
              // gives a distance != 0.0 only if distance to all lines could be determined (which means the projection of the point would lie
              // within the side)
              // if the distance could not be determined, then the local coordinates are (0,0,0) and also smaller than MINIMALTOL
              s->LocalCoordinates( xyz, rst );
              double d = rst( 2, 0 );
              if ( fabs( d ) > MINIMALTOL )
              {
                if( fabs(smallest_dist) > MINIMALTOL ) // distance for node, this facet and another side already set, distance already set by another side and this facet
                {
                    if ( (d > 0) and (fabs( d ) < fabs( smallest_dist )) ) // new smaller distance found for the same facet with another side
                    {
#ifdef DEBUG
                        if( pos == Point::inside) cout << "!!! position of node " << n->Id()  << " has changed from inside to outside" << endl;
#endif
                        // set new position
                        pos = Point::outside;

                        // set new smallest distance
                        smallest_dist = d;
                    }
                    else if((d < 0) and (fabs( d ) < fabs( smallest_dist ))) //new smaller distance found for the same facet with another side
                    {
#ifdef DEBUG
                        if( pos == Point::outside) cout << "!!! position of node " << n->Id()  << " has changed from outside to inside" << endl;
#endif

                        // set new position
                        pos = Point::inside;

                        // set new smallest distance
                        smallest_dist = d;
                    }
                }
                else // standard case : distance set for the first time (smallest_dist currently 0.0)
                {
                    if ( (d > 0) )
                    {
                        pos = Point::outside;
                        smallest_dist = d;
                    }
                    else // d<0
                    {
                    	pos = Point::inside;
                        smallest_dist = d;
                    }
                }
              }
              else // d=0 or distance smaller than MINIMALTOL
              {
                // within the cut plane but not cut by the side
                break;
              }
            }
            done = true;
//            break; // do not finish the loop over sides (cut_faces_)
          }
        }
        if ( done )
        {
          // set the final found position
          if(pos != Point::undecided ) p->Position(pos);

          break;
        }
      } // end for facets
      if ( p->Position()==Point::undecided )
      {
        // Still undecided! No facets with cut side attached! Will be set in a
        // minute.
      }
    } // end if undecided
    else if ( pos==Point::outside or pos==Point::inside )
    {
      // The nodal position is already known. Set it to my facets. If the
      // facets are already set, this will not have much effect anyway. But on
      // multiple cuts we avoid unset facets this way.
      const plain_facet_set & facets = p->Facets();
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        f->Position( pos );
      }
    }

  } // loop nodes

}

// Uli's original version
//
//void GEO::CUT::Element::FindNodePositions()
//{
//  LINALG::Matrix<3,1> xyz;
//  LINALG::Matrix<3,1> rst;
//
//  const std::vector<Node*> & nodes = Nodes();
//  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
//  {
//    Node * n = *i;
//    Point * p = n->point();
//    Point::PointPosition pos = p->Position();
//    if ( pos==Point::undecided )
//    {
//      bool done = false;
//      const plain_facet_set & facets = p->Facets();
//      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
//      {
//        Facet * f = *i;
//        for ( plain_side_set::const_iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
//        {
//          Side * s = *i;
//
//          // Only take a side that belongs to one of this points facets and
//          // shares a cut edge with this point. If there are multiple cut
//          // sides within the element (facets), only the close one will always
//          // give the right direction.
//          if ( f->IsCutSide( s ) and p->CommonCutEdge( s )!=NULL )
//          {
//            if ( p->IsCut( s ) )
//            {
//              p->Position( Point::oncutsurface );
//            }
//            else
//            {
//              p->Coordinates( xyz.A() );
//              s->LocalCoordinates( xyz, rst );
//              double d = rst( 2, 0 );
//              if ( fabs( d ) > MINIMALTOL )
//              {
//                if ( d > 0 )
//                {
//                  p->Position( Point::outside );
//                }
//                else
//                {
//                  p->Position( Point::inside );
//                }
//              }
//              else
//              {
//                // within the cut plane but not cut by the side
//                break;
//              }
//            }
//            done = true;
//            break;
//          }
//        }
//        if ( done )
//          break;
//      }
//      if ( p->Position()==Point::undecided )
//      {
//        // Still undecided! No facets with cut side attached! Will be set in a
//        // minute.
//      }
//    }
//    else if ( pos==Point::outside or pos==Point::inside )
//    {
//      // The nodal position is already known. Set it to my facets. If the
//      // facets are already set, this will not have much effect anyway. But on
//      // multiple cuts we avoid unset facets this way.
//      const plain_facet_set & facets = p->Facets();
//      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
//      {
//        Facet * f = *i;
//        f->Position( pos );
//      }
//    }
//  }
//}


bool GEO::CUT::Element::IsCut()
{
  if ( cut_faces_.size()>0 )
  {
    return true;
  }
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.IsCut() )
    {
      return true;
    }
  }
  return false;
}

bool GEO::CUT::Element::OnSide( Facet * f )
{
  if ( not f->HasHoles() )
  {
    return OnSide( f->Points() );
  }
  return false;
}

bool GEO::CUT::Element::OnSide( const std::vector<Point*> & facet_points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Point*>::const_iterator i=facet_points.begin();
        i!=facet_points.end();
        ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }

  PointSet points;
  std::copy( facet_points.begin(), facet_points.end(),
             std::inserter( points, points.begin() ) );

  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.OnSide( points ) )
    {
      return true;
    }
  }

  return false;
}


void GEO::CUT::Element::GetIntegrationCells( plain_integrationcell_set & cells )
{
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = *i;
    vc->GetIntegrationCells( cells );
  }
}

void GEO::CUT::Element::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( cut_faces_.count( f->ParentSide() )!= 0 )
    {
      f->GetBoundaryCells( bcells );
    }
  }
}

void GEO::CUT::Element::GetCutPoints( PointSet & cut_points )
{
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side * side = *i;

    for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side * other = *i;
      side->GetCutPoints( this, *other, cut_points );
    }
  }
}

void GEO::CUT::Element::CreateIntegrationCells( Mesh & mesh, int count, bool levelset )
{
  if ( not active_ )
    return;

#ifdef DEBUGCUTLIBRARY
  {
    int volume_count = 0;
    for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
    {
      VolumeCell * vc = *i;

      std::stringstream str;
      str << "volume-" << count << "-" << volume_count << ".plot";
      std::ofstream file( str.str().c_str() );
      vc->Print( file );
      volume_count += 1;
    }
  }
#endif

  if ( cells_.size()==1 )
  {
    VolumeCell * vc = *cells_.begin();
    if ( IntegrationCellCreator::CreateCell( mesh, Shape(), vc ) )
    {
      CalculateVolumeOfCellsTessellation();
      return;
    }
  }

  if ( mesh.CreateOptions().SimpleShapes() )
  {
    if ( IntegrationCellCreator::CreateCells( mesh, this, cells_ ) )
    {
      CalculateVolumeOfCellsTessellation();
      return;
    }
  }

  PointSet cut_points;

  // There are never holes in a cut facet. Furthermore, cut facets are
  // always convex, as all elements and sides are convex. Thus, we are free
  // to triangulate all cut facets. This needs to be done, so repeated cuts
  // work in the right way.

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->OnCutSide() and f->HasHoles() )
      throw std::runtime_error( "no holes in cut facet possible" );
    //f->GetAllPoints( mesh, cut_points, f->OnCutSide() );
#if 1
    f->GetAllPoints( mesh, cut_points, levelset and f->OnCutSide() );
#else
    f->GetAllPoints( mesh, cut_points, false );
#endif
  }

  std::vector<Point*> points;
  points.reserve( cut_points.size() );
  points.assign( cut_points.begin(), cut_points.end() );

  // sort points that go into qhull to obtain the same result independent of
  // pointer values (compiler flags, code structure, memory usage, ...)
  std::sort( points.begin(), points.end(), PointPidLess() );

#if 0
  {
    LINALG::Matrix<3,1> xyz;
    LINALG::Matrix<3,1> rst;
    for ( std::vector<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
    {
      Point * p = *i;
      std::copy( p->X(), p->X()+3, &xyz( 0, 0 ) );
      LocalCoordinates( xyz, rst );
      std::cout << "rst[" << p->Id() << "] = ("
                << std::setprecision( 10 )
                << rst( 0, 0 ) << ","
                << rst( 1, 0 ) << ","
                << rst( 2, 0 ) << ")\n";
    }
    throw std::runtime_error( "debug output done" );
  }
#endif

  TetMesh tetmesh( points, facets_, false );
  tetmesh.CreateElementTets( mesh, this, cells_, cut_faces_, count, levelset );

  CalculateVolumeOfCellsTessellation();
}

void GEO::CUT::Element::RemoveEmptyVolumeCells()
{
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); )
  {
    VolumeCell * vc = *i;
    if ( vc->Empty() )
    {
      vc->Disconnect();
      set_erase( cells_, i );
    }
    else
    {
      ++i;
    }
  }
}

void GEO::CUT::Element::MakeVolumeCells( Mesh & mesh )
{
#if 0
#ifdef DEBUGCUTLIBRARY
  DumpFacets();
#endif
#endif

#if 0
  VolumeCellGenerator vcg( sides_, facets_ );
  vcg.CreateVolumeCells( mesh, this, cells_ );
#else
  FacetGraph fg( sides_, facets_ );
  fg.CreateVolumeCells( mesh, this, cells_ );
#endif
}

bool GEO::CUT::ConcreteElement<DRT::Element::tet4>::PointInside( Point* p )
{
  Position<DRT::Element::tet4> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::hex8>::PointInside( Point* p )
{
  Position<DRT::Element::hex8> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::wedge6>::PointInside( Point* p )
{
  Position<DRT::Element::wedge6> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::PointInside( Point* p )
{
  Position<DRT::Element::pyramid5> pos( *this, *p );
  return pos.Compute();
}


void GEO::CUT::ConcreteElement<DRT::Element::tet4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::tet4> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::hex8>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex8> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::wedge6>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::wedge6> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::pyramid5> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

int GEO::CUT::Element::NumGaussPoints( DRT::Element::DiscretizationType shape )
{
  int numgp = 0;
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = *i;
    numgp += vc->NumGaussPoints( shape );
  }
  return numgp;
}

void GEO::CUT::Element::DebugDump()
{
  std::cout << "Problem in element " << Id() << " of shape " << Shape() << ":\n";
  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    //std::cout << n->LSV();
    n->Plot( std::cout );
  }
  std::cout << "\n";
  const plain_side_set & cutsides = CutSides();
  for ( plain_side_set::const_iterator i=cutsides.begin(); i!=cutsides.end(); ++i )
  {
    Side * s = *i;
    //s->Print();
    const std::vector<Node*> & side_nodes = s->Nodes();
    for ( std::vector<Node*>::const_iterator i=side_nodes.begin(); i!=side_nodes.end(); ++i )
    {
      Node * n = *i;
      n->Plot( std::cout );
    }
    std::cout << "\n";
  }
}

void GEO::CUT::Element::GnuplotDump()
{
  std::stringstream str;
  str << "element" << Id() << ".plot";
  std::ofstream file( str.str().c_str() );

  plain_edge_set all_edges;

  const std::vector<Side*> & sides = Sides();
  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    const std::vector<Edge*> & edges = s->Edges();

    std::copy( edges.begin(), edges.end(), std::inserter( all_edges, all_edges.begin() ) );
  }

  for ( plain_edge_set::iterator i=all_edges.begin(); i!=all_edges.end(); ++i )
  {
    Edge * e = *i;
    e->BeginNode()->point()->Plot( file );
    e->EndNode()  ->point()->Plot( file );
    file << "\n\n";
  }
}

void GEO::CUT::Element::DumpFacets()
{
  std::stringstream str;
  str << "facets" << Id() << ".plot";
  std::string name = str.str();

  std::cout << "write '" << name << "'\n";
  std::ofstream file( name.c_str() );

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Print( file );
  }
}

/*-----------------------------------------------------------------*
 * Calculate volume of all volumecells when Tessellation is used
 *-----------------------------------------------------------------*/
void GEO::CUT::Element::CalculateVolumeOfCellsTessellation()
{
  const plain_volumecell_set& volcells = VolumeCells();
  for(plain_volumecell_set::const_iterator i=volcells.begin();i!=volcells.end();i++)
  {
    VolumeCell *vc1 = *i;
    plain_integrationcell_set ics;
    vc1->GetIntegrationCells(ics);

    double volume = 0;
    for ( plain_integrationcell_set::iterator j=ics.begin(); j!=ics.end(); ++j )
    {
      IntegrationCell * ic = *j;
      volume += ic->Volume();
    }

    vc1->SetVolume(volume);
  }
}

/*------------------------------------------------------------------------------------------------------------------*
      The Gauss rules for each cut element is constructed by performing moment fitting for each volumecells.
           Unless specified moment fitting is performed only for cells placed in the fluid region
*-------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::Element::MomentFitGaussWeights( Mesh & mesh, bool include_inner, std::string Bcellgausstype )
{
  if ( not active_ )
    return;

  //When the cut side touches the element the shape of the element is retained
  /*if(cells_.size()==1)
  {
	  VolumeCell * vc = *cells_.begin();
      if ( IntegrationCellCreator::CreateCell( mesh, Shape(), vc ) )
      {
           return;
      }
  }*/

 /* if ( mesh.CreateOptions().SimpleShapes() )
    {
      if ( IntegrationCellCreator::CreateCells( mesh, this, cells_ ) )
      {
        return;
      }
    }*/

  for(plain_volumecell_set::iterator i=cells_.begin();
                           i!=cells_.end();i++)
  {
    VolumeCell *cell1 = *i;
    cell1->MomentFitGaussWeights(this, mesh, include_inner, Bcellgausstype);
  }
}

/*------------------------------------------------------------------------------------------------------------------*
   The Gauss rules for each cut element is constructed by triangulating the facets and applying divergence theorem
            Unless specified moment fitting is performed only for cells placed in the fluid region
*-------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::Element::DirectDivergenceGaussRule( Mesh & mesh, bool include_inner, std::string Bcellgausstype )
{
  if ( not active_ )
    return;

  for(plain_volumecell_set::iterator i=cells_.begin();
                           i!=cells_.end();i++)
  {
    VolumeCell *cell1 = *i;
    cell1->DirectDivergenceGaussRule(this, mesh, include_inner, Bcellgausstype);
  }
}
