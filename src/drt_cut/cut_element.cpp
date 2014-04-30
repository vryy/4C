#include "cut_intersection.H"
#include "cut_position.H"
#include "cut_tetmesh.H"
#include "cut_options.H"
#include "cut_integrationcellcreator.H"
#include "cut_volumecellgenerator.H"
#include "cut_facetgraph.H"

//#include "cut_createVolumeCell.H"

#include <string>
#include <stack>

#include <Teuchos_TimeMonitor.hpp>

/*--------------------------------------------------------------------*
 * struct for comparison of position of sides using ray-tracing techniques
 * shoot a ray starting from startpoint through the midpoint of one side
 * find the intersection point with the second side
 * dependent on the local coordinates along the ray, decide which
 * side lies in front of the other one
 *--------------------------------------------------------------------*/
struct nextSideAlongRay {
  nextSideAlongRay(GEO::CUT::Point* startpoint, GEO::CUT::Point* cutpoint) : startpoint_(startpoint),cutpoint_(cutpoint) {
    startpoint_->Coordinates(startpoint_xyz_.A());
    cutpoint_->Coordinates(cutpoint_xyz_.A());
  };


  /*--------------------------------------------------------------------*
   * check if both sides have the same normal vector
   *--------------------------------------------------------------------*/
  bool SameNormal(GEO::CUT::Side* s1, GEO::CUT::Side* s2, const LINALG::Matrix<3,1> & cutpoint_xyz)
  {

    LINALG::Matrix<3,1> rst(true);
    LINALG::Matrix<2,1> rs(true);

    //-------------
    // first side
    s1->LocalCoordinates(cutpoint_xyz, rst, false);

    rs(0)=rst(0);
    rs(1)=rst(1);

    LINALG::Matrix<3,1> normal_1(true);
    s1->Normal(rs, normal_1);

    //-------------
    // second side
    s2->LocalCoordinates(cutpoint_xyz, rst, false);

    rs(0)=rst(0);
    rs(1)=rst(1);

    LINALG::Matrix<3,1> normal_2(true);
    s2->Normal(rs, normal_2);

    //-------------
    if(normal_1.Dot(normal_2) > 1-TOLERANCE) return true;

    return false;
  }

  /*--------------------------------------------------------------------*
   * comparator function to sort two sides, which side lies in front of the other along the ray
   *--------------------------------------------------------------------*/
  bool operator()(GEO::CUT::Side* s1, GEO::CUT::Side* s2)
  {

    // REMARK:
    // shoot a ray through the first side s1 starting from startpoint and find intersection with side s2
    // if not successful shoot a second ray through side s2 and find intersection with s1
    // if not successful check if the sides are parallel
    bool is_closer = false;

    if(s1->IsCloserSide(startpoint_xyz_, s2, is_closer))
    {
      if(is_closer) return true;
      else return false;
    }
    else if(s2->IsCloserSide(startpoint_xyz_, s1, is_closer))
    {
      if(!is_closer) return true;
      else return false;
    }
    else if(SameNormal(s1, s2, cutpoint_xyz_)) // check if both sides are parallel to each other, then both sides lead the same position
    {
      return true;
    }
    else
    {
      //TODO: check if we can relax this case (Parallelogramm, both Positions would be the same!?)
      // try to return true or false, both sides should lead to the same position, sorting not necessary

      //return true;

      throw std::runtime_error("ray-tracing-based comparisons to find the nearest side along the ray failed for the first time!");
    }

    return false;
  }


   GEO::CUT::Point* startpoint_;
   GEO::CUT::Point* cutpoint_;

   LINALG::Matrix<3,1> startpoint_xyz_;
   LINALG::Matrix<3,1> cutpoint_xyz_;
};

 /*-----------------------------------------------------------------------------------*
  *  For this shadow element, set corner nodes of parent Quad element      sudhakar 11/13
  *-----------------------------------------------------------------------------------*/
  void GEO::CUT::Element::setQuadCorners( Mesh & mesh, const std::vector<int> & nodeids )
  {
    if( not isShadow_ )
      dserror( "You can't set Quad-corners for non-shadow element\n" );

    for( unsigned i=0; i<nodeids.size(); i++ )
    {
      Node * n1 = mesh.GetNode( nodeids[i] );
      quadCorners_.push_back( n1 );
    }
  }

  /*----------------------------------------------------------------------------------------------------*
   *  Get corner nodes of parent Quad element from which this shadow element is derived       sudhakar 11/13
   *----------------------------------------------------------------------------------------------------*/
  std::vector<GEO::CUT::Node*> GEO::CUT::Element::getQuadCorners()
  {
    if( (not isShadow_) or quadCorners_.size() == 0 )
      dserror("what?! you want Quadratic element corners for linear element?\n");
    return quadCorners_;
  }

/*--------------------------------------------------------------------*
 *            cut this element with given cut_side
 *--------------------------------------------------------------------*/
bool GEO::CUT::Element::Cut( Mesh & mesh, Side & cut_side, int recursion )
{
  bool cut = false;

  // find nodal points inside the element
  const std::vector<Node*> side_nodes = cut_side.Nodes();

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
    else // point cuts this element, already determined by another side
    {
      cut = true;
    }
  }

  // all the other cut points lie on sides of the element (s is an element side, cut_side is the cutter side)
  const std::vector<Side*> & sides = Sides();
  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    if ( FindCutPoints( mesh, *s, cut_side, recursion ) )
    {
      cut = true;
    }
  }

  // insert this side into cut_faces_
  if ( cut )
  {
    cut_faces_.insert( &cut_side );
    return true;
  }
  else
  {
    return false;
  }
}

/*---------------------------------------------------------------------------*
 * After all cut points are found, create cut lines for this element by
 * connecting appropriate cut points
 *---------------------------------------------------------------------------*/
void GEO::CUT::Element::MakeCutLines( Mesh & mesh, Creator & creator )
{
  for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
  {
    Side & cut_side = **i;

    bool cut = false;

    const std::vector<Side*> & sides = Sides();
    // create cut lines over each side of background element
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
      if ( FindCutLines( mesh, *s, cut_side ) )
      {
        cut = true;
      }
    }

    // find lines inside the element
    const std::vector<Edge*> & side_edges = cut_side.Edges();
    for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
    {
      Edge * e = *i;
      std::vector<Point*> line;
      e->CutPointsInside( this, line );
      mesh.NewLinesBetween( line, &cut_side, NULL, this );
    }

    if ( cut )
    {
      // create any remaining cut lines
      //side.CreateMissingLines( creator, this );
    }
  }
}

/*------------------------------------------------------------------------------------------*
 * Find cut points between a background element side and a cut side
 * Cut points are stored correspondingly
 *------------------------------------------------------------------------------------------*/
bool GEO::CUT::Element::FindCutPoints( Mesh & mesh, Side & ele_side, Side & cut_side, int recursion )
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- FindCutPoints" );


  bool cut = ele_side.FindCutPoints( mesh, this, cut_side, recursion );         // edges of element side cuts through cut side
  bool reverse_cut = cut_side.FindCutPoints( mesh, this, ele_side, recursion ); // edges of cut side cuts through element side
  return cut or reverse_cut;
}

/*------------------------------------------------------------------------------------------*
 *     Returns true if cut lines exist between the cut points produced by the two sides
 *------------------------------------------------------------------------------------------*/
bool GEO::CUT::Element::FindCutLines( Mesh & mesh, Side & ele_side, Side & cut_side )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh --- FindCutLines" );


  bool cut = ele_side.FindCutLines( mesh, this, cut_side );
  bool reverse_cut = cut_side.FindCutLines( mesh, this, ele_side );
  return cut or reverse_cut;
}

/*------------------------------------------------------------------------------------------*
 * Create facets
 *------------------------------------------------------------------------------------------*/
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


/*------------------------------------------------------------------------------------------*
 *     Determine the inside/outside/oncutsurface position for the element's nodes
 *------------------------------------------------------------------------------------------*/
void GEO::CUT::Element::FindNodePositions()
{

  // DEBUG flag for FindNodePositions
  // compute positions for nodes again, also if already set by other nodes, facets, vcs (safety check)
 //  #define check_for_all_nodes

#if(1)
  //----------------------------------------------------------------------------------------
  // new implementation based on cosine between normal vector on cut side and line-vec between point and cut-point
  //----------------------------------------------------------------------------------------

  const std::vector<Node*> & nodes = Nodes();

  // determine positions for all the element's nodes
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();

#ifdef check_for_all_nodes
    std::cout << "Position for node " << n->Id() << std::endl;
    // do the computation in all cases
#else
    if ( pos==Point::undecided )
#endif
    {
      bool done = false;

      // a) this line lies on the cut surface, then p has to be on cut surface for a least one cut-side
      // b) the line connects two points, both lying on different cut sides, then p has to be on cut surface for a least one cut-side

      // check if the node lies on a cut-surface
      for ( plain_side_set::const_iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
      {
        Side * s = *i;

        // check if the point lies on one of the element's cut sides
        if ( p->IsCut( s ) )
        {
          p->Position( Point::oncutsurface );

          done = true;
          break;
        }
      }

      if(done) continue; // next node

      // c) search for a line connection between the point p and a cut side in this element
      //
      // is there a facet's (!) line between the point p and a cut-point on the side s ?
      // if there is a line, then no further point lies between p and the cut side
      // this line goes either through the outside or inside region

      const plain_facet_set & facets = p->Facets();

      // loop all the facets sharing this node
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;

        // loop all the cut-faces stored for this element
        // (includes cut-faces that only touch the element at points, edges, sides or parts of them)
        for ( plain_side_set::const_iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
        {
          Side * s = *i;

          // is there a common point between facet and side?
          // and belongs the facet to the element (otherwise we enter the neighboring element via the facet)?
          // - however we include cut-sides of neighboring elements that only touches the facet,
          // - the facet however has to an element's facet
          if( f->IsCutSide( s ) and IsFacet(f) )
          {
            // for inside-outside decision there must be a direct line connection between the point and the cut-side
            // check for a common facet's line between a side's cut point and point p

            std::map<std::pair<Point*, Point*>, plain_facet_set > lines;
            f->GetLines(lines); // all facet's lines, each line sorted by P1->Id() < P2->Id()

            for(std::map<std::pair<Point*,Point*>, plain_facet_set>::iterator line_it = lines.begin();
                line_it != lines.end();
                line_it++)
            {
              std::pair<Point*,Point*> line = line_it->first;

              Point* cutpoint = NULL;

              // find the right facet's line and the which endpoint is the cut-point
              if(line.first->Id() == p->Id() and line.second->IsCut(s) )
              {
                cutpoint = line.second;
              }
              else if(line.second->Id() == p->Id() and line.first->IsCut(s))
              {
                cutpoint = line.first;
              }
              else
              {
                // this line is not a line between the point and the cut-side
                // continue with next line
                continue;
              }

              //---------------------------------------------------
              // call the main routine to compute the position based on the angle between
              // the line-vec (p-c) and an appropriate cut-side
              done = ComputePosition(p, cutpoint, f, s);
              //---------------------------------------------------

              if(done) break;
            } // end lines

            if(done) break;
          }
        } // end cutsides

        if(done) break;
      } // loop facets
      if ( p->Position()==Point::undecided )
      {
        // Still undecided! No facets with cut side attached! Will be set in a
        // minute.
      }

      if(done) continue;
    } // end if undecided

#ifdef check_for_all_nodes
    if ( pos==Point::outside or pos==Point::inside )
#else
      else if ( pos==Point::outside or pos==Point::inside )
#endif
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
    } // end if outside or inside

  } // loop nodes

#else

  //----------------------------------------------------------------------------------------
  // improved old implementation, causes wrong positions especially when line, sides are cut more than one
  //----------------------------------------------------------------------------------------

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
                        if( pos == Point::inside) std::cout << "!!! position of node " << n->Id()  << " has changed from inside to outside" << std::endl;
#endif
                        // set new position
                        pos = Point::outside;

                        // set new smallest distance
                        smallest_dist = d;
                    }
                    else if((d < 0) and (fabs( d ) < fabs( smallest_dist ))) //new smaller distance found for the same facet with another side
                    {
#ifdef DEBUG
                        if( pos == Point::outside) std::cout << "!!! position of node " << n->Id()  << " has changed from outside to inside" << std::endl;
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
#endif
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


/*------------------------------------------------------------------------------------------*
 *  main routine to compute the position based on the angle between the line-vec (p-c) and an appropriate cut-side
 *------------------------------------------------------------------------------------------*/
bool GEO::CUT::Element::ComputePosition(Point * p,                        // the point for that the position has to be computed
                                        Point * cutpoint,                 // the point on cut side which is connected to p via a facets line)
                                        Facet * f,                        // the facet via which p and cutpoint are connected
                                        Side* s                           // the current cut side, the cutpoint lies on
)
{

  // REMARK: the following inside/outside position is based on comparisons of the line vec between point and the cut-point
  // and the normal vector w.r.t cut-side "angle-comparison"
  // in case the cut-side is not unique we have to determine at least one cut-side which can be used for the "angle-criterion"
  // such that the right position is guaranteed
  //
  // in case that the cut-point lies on an edge between two different cut-sides or even on a node between several cut-sides,
  // then we have to find the side (maybe not unique, but at least one of these) that defines the right position
  // based on the "angle-criterion"

  //---------------------------
  // find the element's volume-cell the cut-side and the line has to be adjacent to
  const plain_volumecell_set & facet_cells = f->Cells();
  plain_volumecell_set adjacent_cells;

  for(plain_volumecell_set::const_iterator f_cells_it= facet_cells.begin(); f_cells_it!=facet_cells.end(); f_cells_it++)
  {
    if(cells_.count(*f_cells_it)) adjacent_cells.insert(*f_cells_it); // insert this cell
  }

  if(adjacent_cells.size()> 1)
  {
    std::cout << "Warning: there is not a unique element's volumecell, number="  << adjacent_cells.size() << ", the line and facet is adjacent to-> Check this" << std::endl;
    throw std::runtime_error("Warning: there is not a unique element's volumecell");
  }
  else if(adjacent_cells.size()==0)
  {
    std::cout << "facet cells" << facet_cells.size() << std::endl;
    std::cout << "Warning: there is no adjacent volumecell, the line and facet is adjacent to-> Check this" << std::endl;
    throw std::runtime_error("Warning: there is no adjacent volumecell");
  }

  // that's the adjacent volume-cell
  VolumeCell* vc = adjacent_cells[0];


  //---------------------------
  // get the element's cut-sides adjacent to this cut-point and adjacent to the same volume-cell
  const plain_side_set & cut_sides = this->CutSides();
  plain_side_set point_cut_sides;

  //  plain_side_set e_cut_sides;
  for(plain_side_set::const_iterator side_it = cut_sides.begin(); side_it!=cut_sides.end(); side_it++)
  {
    // is the cut-point a cut-point of this element's cut_side?
    // and is this side a cut-side adjacent to this volume-cell ?
    // and remove sides, if normal vector is orthogonal to side and cut-point lies on edge, since then the angle criterion does not work
    if(     cutpoint->IsCut(*side_it)
        and vc->IsCut(*side_it)
        and !IsOrthogonalSide(*side_it, p, cutpoint) )
    {
      // the angle-criterion has to be checked for this side
      point_cut_sides.insert(*side_it);
    }
  }

  //std::cout << "how many cut_sides found? " << point_cut_sides.size() << std::endl;


  if(point_cut_sides.size()==0)
  {
    // no right cut_side found! -> Either another node can compute the position or hope for distributed positions or hope for parallel communication
    return false;
  }

  //------------------------------------------------------------------------
  // sort the sides and do the check for the first one!
  // the sorting is based on ray-tracing techniques:
  // shoot a ray starting from point p through the midpoint of one of the two sides and find another intersection point
  // the local coordinates along this ray determines the order of the sides
  //------------------------------------------------------------------------
  if( point_cut_sides.size() > 1 )
    std::sort(point_cut_sides.begin(), point_cut_sides.end(), nextSideAlongRay(p, cutpoint));


  //------------------------------------------------------------------------
  // determine the inside/outside position w.r.t the chosen cut-side
  // in case of the right side the "angle-criterion" leads to the right decision (position)
  Side* cut_side = point_cut_sides[0];

  bool successful = PositionByAngle(p, cutpoint, cut_side);
  //------------------------------------------------------------------------


  //if(successful) std::cout << "set position to " << p->Position() << std::endl;
  //else std::cout << "not successful" << std::endl;

  return successful;
}


/*------------------------------------------------------------------------------------------*
 *  determine the position of point p based on the angle between the line (p-c) and the side's normal vector, return if successful
 *------------------------------------------------------------------------------------------*/
bool GEO::CUT::Element::PositionByAngle(Point* p, Point* cutpoint, Side* s )
{

  LINALG::Matrix<3,1> xyz(true);
  LINALG::Matrix<3,1> cut_point_xyz(true);

  p->Coordinates( xyz.A() );
  cutpoint->Coordinates( cut_point_xyz.A() );

  //------------------------------------------------------------------------
  // determine the inside/outside position w.r.t the chosen cut-side
  // in case of the right side the "angle-criterion" leads to the right decision (position)

  LINALG::Matrix<2,1> rs(true);      // local coordinates of the cut-point w.r.t side
  double dist = 0.0;                 // just used for within-side-check, has to be about 0.0
  double Tol_2D = 1e-09;             // choose a very rough tolerance to allow "bad" cuts as in cut-tests, just to catch non-sense
  bool within_side = s->WithinSide(cut_point_xyz, rs, dist, Tol_2D);

  if(!within_side)
  {
    std::cout << "Side: " << std::endl; s->Print();
    std::cout << "Point: " << std::endl; p->Print(std::cout);
    std::cout << "local coordinates " << rs << " dist " << dist << std::endl;
    throw std::runtime_error("cut-point does not lie on side! That's wrong, because it is a side's cut-point!");
  }


  LINALG::Matrix<3,1> normal(true);
  s->Normal(rs, normal); // outward pointing normal at cut-point

  LINALG::Matrix<3,1> line_vec(true);
  line_vec.Update(1.0, xyz ,-1.0, cut_point_xyz); // vector representing the line between p and the cut-point

  // check the cosine between normal and line_vec
  double n_norm = normal.Norm2();
  double l_norm = line_vec.Norm2();
  if(n_norm < MINIMALTOL or l_norm < MINIMALTOL)
  {
    dserror(" the norm of line_vec or n_norm is smaller than %d, should these points be one point in pointpool?, lnorm=%d, nnorm=%d", MINIMALTOL, l_norm, n_norm);
  }

  // cosine between the line-vector and the normal vector
  double cosine = normal.Dot(line_vec);
  cosine /= (n_norm * l_norm);

  if( cosine > 0.0)
  {
    p->Position(Point::outside);
    // std::cout << " set position to outside" << std::endl;
    return true;
  }
  else if ( cosine < 0.0 )
  {
    p->Position(Point::inside);
    // std::cout << " set position to inside" << std::endl;
    return true;
  }
  else
  {
    // Still undecided!
    // There must be another side with cosine != 0.0
    return false;
  }
  //------------------------------------------------------------------------
  return false;
}


/*------------------------------------------------------------------------------------------*
 *  check if the side's normal vector is orthogonal to the line between p and the cutpoint
 *------------------------------------------------------------------------------------------*/
bool GEO::CUT::Element::IsOrthogonalSide(Side* s, Point* p, Point* cutpoint)
{
  if(s->OnEdge(cutpoint)) // check if the point lies on at least one edge of the side, otherwise it cannot be orthogonal
  {
    LINALG::Matrix<3,1> line(true);
    LINALG::Matrix<3,1> p_xyz(true);
    LINALG::Matrix<3,1> cut_point_xyz(true);

    p->Coordinates(p_xyz.A());
    cutpoint->Coordinates(cut_point_xyz.A());
    line.Update(1.0, p_xyz, -1.0, cut_point_xyz);

    double line_norm = line.Norm2();

    if(line_norm > TOLERANCE )
    {
      line.Scale(1./line_norm);
    }
    else
    {
      std::cout << "point: " << p_xyz << std::endl;
      std::cout << "cutpoint: " << cut_point_xyz << std::endl;
      dserror("the line has nearly zero length: %d", line_norm);
    }

    if(s->Shape() != DRT::Element::tri3)
    {
      std::cout << "HERE !tri3 cutsides are used!!!" << std::endl;
//      throw std::runtime_error("expect only tri3 cutsides!");
    }

    // tri3/quad4 element center
    LINALG::Matrix<2,1> rs(true);

    if(s->Shape()==DRT::Element::tri3)
    {
      rs = DRT::UTILS::getLocalCenterPosition<2>(DRT::Element::tri3);
    }
    else if(s->Shape()==DRT::Element::quad4)
    {
      rs = DRT::UTILS::getLocalCenterPosition<2>(DRT::Element::quad4);
    }
    else throw std::runtime_error("unsupported side-shape");


    LINALG::Matrix<3,1> normal(true);
    s->Normal(rs, normal);

    // check for angle=+-90 between line and normal
    if( fabs(normal.Dot(line)) < (0.0 + TOLERANCE) ) return true;

  }

  return false;
}


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


/*------------------------------------------------------------------------------------------*
 * Get cutpoints of this element
 *------------------------------------------------------------------------------------------*/
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

void GEO::CUT::Element::CreateIntegrationCells( Mesh & mesh, int count, bool levelset, bool tetcellsonly )
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

  if(not tetcellsonly)
  {
    if ( cells_.size()==1 ) // in case there is only one volumecell, check if a simple shaped integration cell is possible
    {
      VolumeCell * vc = *cells_.begin();
      if ( IntegrationCellCreator::CreateCell( mesh, Shape(), vc ) )
      {
        CalculateVolumeOfCellsTessellation();
        return; // return if this was possible
      }
    }

    if ( mesh.CreateOptions().SimpleShapes() ) // try to create only simple-shaped integration cells for all! volumecells
    {
      if ( IntegrationCellCreator::CreateCells( mesh, this, cells_ ) )
      {
        CalculateVolumeOfCellsTessellation();
        return; // return if this was possible
      }
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

  // standard subtetrahedralization starts here, there also the boundary cells will be created

  TetMesh tetmesh( points, facets_, false );
  tetmesh.CreateElementTets( mesh, this, cells_, cut_faces_, count, levelset, tetcellsonly );

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

/*------------------------------------------------------------------------------------------*
 * Create volumecells
 *------------------------------------------------------------------------------------------*/
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

#if 1
  FacetGraph fg( sides_, facets_ );
  fg.CreateVolumeCells( mesh, this, cells_ );
#endif

#if 0
  volumecellCreator vcc( facets_ );
  vcc.createCells( mesh, this, cells_ );
#endif

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

/*-------------------------------------------------------------------------------------------------------------*
 * Find local coodinates of given point w.r to the parent Quad element														sudhakar 11/13
 *-------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::Element::LocalCoordinatesQuad( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  if( not isShadow_ )
    dserror("This is not a shadow elemenet\n");

  switch( getQuadShape() )
  {
  case DRT::Element::hex20:
  {
    Position<DRT::Element::hex20> pos( quadCorners_, xyz );
    bool success = pos.Compute();
    if( success ){}
    rst = pos.LocalCoordinates();
    break;
  }
  case DRT::Element::hex27:
  {
    Position<DRT::Element::hex27> pos( quadCorners_, xyz );
    bool success = pos.Compute();
    if( success ){}
    rst = pos.LocalCoordinates();
    break;
  }
  case DRT::Element::tet10:
  {
    Position<DRT::Element::tet10> pos( quadCorners_, xyz );
    bool success = pos.Compute();
    if( success ){}
    rst = pos.LocalCoordinates();
    break;
  }
  default:
  {
    dserror("not implemented yet\n");
    break;
  }
  }
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
Integrate pre-defined functions over each volumecell created from this element when using Tessellation
*-------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::Element::integrateSpecificFunctionsTessellation()
{
  for(plain_volumecell_set::iterator i=cells_.begin();
                                     i!=cells_.end();i++)
  {
    VolumeCell *cell1 = *i;
    cell1->integrateSpecificFunctionsTessellation();
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
