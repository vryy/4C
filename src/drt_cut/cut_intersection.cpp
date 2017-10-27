/*----------------------------------------------------------------------*/
/*!
\file cut_intersection.cpp

\brief here the intersection of a (plane) surface with a line is performed

\level 2

\maintainer Christoph Ager, Michael Hiermeier
*/
/*----------------------------------------------------------------------*/

#include "cut_intersection.H"
#include "cut_boundingbox.H"
#include "cut_side.H"
#include "cut_position.H"
#include "cut_utils.H"

#include "../drt_lib/drt_globalproblem.H"

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::IntersectionBase> GEO::CUT::IntersectionBase::Create(
    const DRT::Element::DiscretizationType & edge_type,
    const DRT::Element::DiscretizationType & side_type )
{
  const IntersectionFactory factory;
  return factory.CreateIntersection( edge_type, side_type );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
void GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::SetCoordinates()
{
  GetEdge().Coordinates( xyze_lineElement_ );
  GetSide().Coordinates( xyze_surfaceElement_ );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::ComputeCut( Edge * sedge, Edge * eedge,
         PointSet & ee_cut_points, double & tolerance )
{
  return eedge->ComputeCut( GetMeshPtr(), sedge, & ee_cut_points, tolerance );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
void GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::TestSideEdges( Point * p,
    const LINALG::Matrix<dimedge+dimside,1> & xsi,
    std::vector<Edge*> & edges )
{
  if ( AtEdge( xsi ) )
  {
    const LINALG::Matrix<dimside,1> rs( xsi.A(), true );
    GetSide().EdgeAt( rs, edges );
    for ( std::vector<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      p->AddEdge( e );
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckBoundingBoxOverlap()
{
  Teuchos::RCP<BoundingBox> sbb = Teuchos::rcp( BoundingBox::Create( GetSide() ) );
  Teuchos::RCP<BoundingBox> ebb = Teuchos::rcp( BoundingBox::Create( GetEdge() ) );

  return CheckBoundingBoxOverlap( *sbb, *ebb );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckBoundingBoxOverlap(
         BoundingBox & ebb, BoundingBox & sbb ) const
{
  return ( not sbb.Within( 1.0, ebb ) );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckParallelism(
         std::vector< LINALG::Matrix<dimside, 1> > & side_rs_intersect,
         std::vector< LINALG::Matrix<dimedge, 1> > & edge_r_intersect,
         double & tolerance )
{
  switch ( dimside )
  {
    case 1:
    {
      return CheckCollinearity( side_rs_intersect,
          edge_r_intersect, tolerance );
    }
    case 2:
    {
      return CheckParallelismBetweenSideAndEdge(
          side_rs_intersect, edge_r_intersect, tolerance );
    }
    default:
    {
      dserror( "The given side element type is currently unsupported! \n"
               "( dim = %d | sideType = %s ", dimside,
               DRT::DistypeToString( sidetype ).c_str() );
      exit( EXIT_FAILURE );
    }
  }
  // this cannot be reached
  exit( EXIT_FAILURE );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckParallelismBetweenSideAndEdge(
         std::vector< LINALG::Matrix<dimside, 1> > & side_rs_intersect,
         std::vector< LINALG::Matrix<dimedge, 1> > & edge_r_intersect,
        double & tolerance )
{
  side_rs_intersect.clear();
  edge_r_intersect.clear();
  tolerance = 0.0;

//  if ( CheckAngleCriterionBetweenSideNormalAndEdge() )
//    dserror("This shouldn't happen at this point! -- hiermeier");

  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckAngleCriterionBetweenSideNormalAndEdge()
{
  // calculate the normal at the side element center
  const LINALG::Matrix< dimside, 1 > rst_side_center(
      DRT::UTILS::getLocalCenterPosition<dimside>( sidetype ) );

  LINALG::Matrix< probdim, numNodesSide > side_deriv1;
  LINALG::Matrix< probdim, probdim > xjm;
  LINALG::Matrix< probdim, 1 > normal_center;
  GEO::CUT::EvalDerivsInParameterSpace< probdim, sidetype >(
      xyze_surfaceElement_,rst_side_center,side_deriv1, xjm, NULL,
      &normal_center, NULL, true );

  // calcualte the direction vector of the edge element
  LINALG::Matrix< probdim, 1 > dedge( & xyze_lineElement_( 0, 0 ), false );
  const LINALG::Matrix< probdim, 1 > e_endpoint( & xyze_lineElement_( 0, 1 ), true );
  dedge.Update( -1.0, e_endpoint, 1.0 );
  const double e_nrm2 = dedge.Norm2();
  if ( e_nrm2 == 0.0 )
    dserror( "The 1-st edge length is zero!" );
  else
    dedge.Scale( 1.0 / e_nrm2 );

  // calculate the inner product and check the angle between the normal and
  // the edge

  const double inner_product = dedge.Dot( normal_center );
  /* If the angle between the normal and the edge is smaller than 89°
   * ( cos( 89° ) = 0.017452406... ), the vectors are definitely not parallel,
   * otherwise it's possible and we do another test. Note, that both vectors have
   * unit length! */
  {
    double a = std::acos( std::abs( inner_product ) ) * 180 / ( 2.0 * std::acos( 0.0 ) );
    std::cout << "angle between side normal and edge = " << a << "°" << std::endl;
  }
  return ( std::abs( inner_product ) <=  0.01745 );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckCollinearity(
         std::vector< LINALG::Matrix<dimside, 1> > & side_rs_corner_intersect,
         std::vector< LINALG::Matrix<dimedge, 1> > & edge_r_corner_intersect,
         double & tolerance )
{
  if ( numNodesEdge != 2 or numNodesSide != 2 )
    dserror( "Two line2 elements are expected, but instead a %s (edge) and %s (side) "
        "element were given.", DRT::DistypeToString( edgetype ).c_str(),
        DRT::DistypeToString( sidetype ).c_str() );

  side_rs_corner_intersect.clear();
  edge_r_corner_intersect.clear();


  // quick check first -- check angle between the two edges, if the check fails,
  // the two edges are definitely not parallel.
  if ( not CheckAngleCriterionBetweenTwoEdges() )
    return false;

  /* If the lines seem to be parallel, we do the expensive calculation and
   * check if the end points of the edge are on the side edge ( and are collinear ). */
  bool is_collinear = true;

  std::vector< LINALG::Matrix<probdim, 1> > side_xyz_corner_intersect;
  side_xyz_corner_intersect.reserve( 2 );

  side_rs_corner_intersect.reserve( 2 );
  edge_r_corner_intersect.reserve( 2 );

  LINALG::TMatrix<LINALG::Matrix<2,1>, 2, 1> e_corner_distance;
  for ( unsigned i = 0; i < 2; ++i )
  {
    const LINALG::Matrix<probdim, 1> e_cornerpoint( & xyze_lineElement_( 0, i ), true );
    Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create(
        xyze_surfaceElement_, e_cornerpoint, sidetype );

    bool withinlimits = pos->Compute( true );

    pos->Distance( e_corner_distance( i ) );
    const double pos_tolerance = pos->NewtonTolerance();

    // keep the largest up-coming tolerance
    if ( i == 0 or tolerance < pos_tolerance )
      tolerance = pos_tolerance;

    if ( e_corner_distance( i ).Norm2() < pos_tolerance )
    {
      // the edge end point seems to lie on the side edge, and the end-point is
      // also within the side limits
      if ( withinlimits )
      {
        side_xyz_corner_intersect.push_back( e_cornerpoint );
        side_rs_corner_intersect.push_back( LINALG::Matrix<dimside,1>( true ) );
        pos->LocalCoordinates( *( side_rs_corner_intersect.end()-1 ) );
      }
    }
    // if one distance is larger than the given tolerance, the edge cannot be
    // collinear
    else
      is_collinear = false;
  }
  // set number of cut points
  num_cut_points_ = side_rs_corner_intersect.size();

  // The two edges are collinear and we found no end points of the edge on
  // the side.
  // --> No more intersections possible for two edges besides the already
  // existing end-points of the side edge.
  if ( is_collinear and side_rs_corner_intersect.size() == 0 )
  {
    return true;
  }

  // End points of the edge on the side edge. In this case
  // we are done, as well, and it doesn't matter if the two objects
  // are parallel or not.
  if ( num_cut_points_ > 0 )
  {
    for ( unsigned i=0; i<num_cut_points_; ++i )
    {
      xsi_side_.Update( side_rs_corner_intersect[0] );
      if ( not FindLocalCoordinateOfEdgeEndPoint( xsi_edge_(0),
          side_xyz_corner_intersect[i], tolerance ) )
        dserror("We couldn't find the correct edge end-point!");
      if ( num_cut_points_ > 1 )
        edge_r_corner_intersect.push_back( LINALG::Matrix<dimedge,1>( xsi_edge_.A(), false ) );

      // safety
      if ( not ( LineWithinLimits() and SurfaceWithinLimits() ) )
        dserror("Something went wrong!");
    }

    return true;
  }

  if ( not is_collinear )
  {
    double prod = 0.0;
    // check the signs of the distance values
    for ( unsigned i=0; i< ( probdim-dimside ); ++i)
    {
      prod =  e_corner_distance( 0 )( i ) * e_corner_distance( 1 )( i );
      // if the off-set direction is exactly in one of the normal directions, we
      // use the second criterion to detect numerical artifacts ( values close
      // to zero but different signs )
      if ( ( prod < 0.0 ) and
           std::abs( e_corner_distance( 0 )( i ) ) > tolerance and
           std::abs( e_corner_distance( 1 )( i ) ) > tolerance )
      {
        // there is a intersection possible, although the lines are almost
        // parallel
        return false;
      }
    }
    // all distance values have the same sign ( or are close to zero ) and the
    // edges are not collinear
    // --> no intersection possible
    return true;
  }

  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::CheckAngleCriterionBetweenTwoEdges()
{
  LINALG::Matrix< probdim, 1 > deedge( & xyze_lineElement_( 0, 0 ), false );
  LINALG::Matrix< probdim, 1 > dsedge( & xyze_surfaceElement_( 0, 0 ), false );

  const LINALG::Matrix< probdim, 1 > e_endpoint( & xyze_lineElement_( 0, 1 ), true );
  deedge.Update( -1.0, e_endpoint, 1.0 );
  const double e_nrm2 = deedge.Norm2();
  if ( e_nrm2 == 0.0 )
    dserror( "The 1-st edge length is zero!" );

  const LINALG::Matrix< probdim, 1 > s_endpoint( & xyze_surfaceElement_( 0, 1 ), true );
  dsedge.Update( -1.0, s_endpoint, 1.0 );
  const double s_nrm2 = dsedge.Norm2();
  if ( s_nrm2 == 0.0 )
    dserror( "The 2-nd edge length is zero!" );

  const double inner_product = deedge.Dot( dsedge );
  /* If the angle between the two lines is larger than 1° ( cos( 1° ) = 0.99847695... ),
   * the two edges are definitely not parallel, otherwise it's possible and
   * we take a closer look. */
  if ( debug )
  {
    std::cout << e_nrm2 << ", " << s_nrm2 << ", " << inner_product << std::endl;
    double a = std::abs( inner_product ) / ( e_nrm2 * s_nrm2 );
    a = std::acos( ( a > 1.0 ? 1.0 : a ) )
               * 180 / ( 2.0 * std::acos( 0.0 ) );
    std::cout << "angle between edge and edge = " << a << "°" << std::endl;
  }
  return ( std::abs( inner_product ) >= e_nrm2 * s_nrm2 * 0.9998 );
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::FindLocalCoordinateOfEdgeEndPoint(
         double & pos,
         const LINALG::Matrix<probdim, 1> & xyz,
         const double & tolerance ) const
{
  const double r_endpoints[2] = { -1.0, 1.0 };
  LINALG::Matrix<probdim,1> dist;
  LINALG::Matrix<numNodesEdge,1> lineFunct;

  for ( unsigned i=0; i<2; ++i )
  {
    pos = r_endpoints[ i ];

    DRT::UTILS::shape_function_1D( lineFunct, pos, edgetype );
    dist.MultiplyNN( xyze_lineElement_, lineFunct );

    dist.Update( 1.0, xyz, -1.0 );

    if ( dist.Norm2() < tolerance )
      return true;
  }
  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge,unsigned dimside,
          unsigned numNodesEdge,unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::Intersect( PointSet & cuts )
{
  CheckInit();
  // ------------------------------------------------------------------------
  // (1) bounding box check
  // ------------------------------------------------------------------------
  bool debugbb = CheckBoundingBoxOverlap();
  if ( debugbb )
  {
    if (debug)
    {
      std::cout << "No Cut Point, as not even the Bounding Boxes are "
          "Overlapping! ( Codeline:" << __LINE__ << " )" << std::endl;
    }

    /* at the moment we want to test as many intersections as
     * possible, that's why bounding boxes are not used yet! */
    if ( UseBoundingBox() )
    {
      return false;   /* not even the bounding boxes are overlapping, no need to
                       * search for an intersection point! */
    }
  }

  // ------------------------------------------------------------------------
  // (2) line || side
  // ------------------------------------------------------------------------
  bool success = false;
  Point * nodalpoint = NULL;

  /* all information gathered by Compute Distance is gathered here!!!
   *   0... begin point
   *   1... end point of line!!! */
  bool zeroarea = false;
  //local coordinates are inside element
  LINALG::TMatrix<bool,2,1> lineendpoint_within_surfacelimits( true );
  //point is really inside element (normal distance = 0)
  LINALG::TMatrix<bool,2,1> lineendpoint_in_surface( true );
  LINALG::Matrix<2,1> lineendpoint_dist;
  LINALG::Matrix<2,1> lineendpoint_tol;
  LINALG::TMatrix<bool,2,1> lineendpoint_conv( true );
  LINALG::Matrix<dimside+dimedge,2> lineendpoint_xsi( true );

  /* store the intersection xsi before it is overwritten by
   * computelineendpoint_distance for later */
  LINALG::Matrix<dimside+dimedge,1> i_xsi;

  //loop over line endpoints
  for (unsigned lp = 0; lp < 2 ; ++lp)
  {
    GEO::CUT::Point * actpoint = 0;
    if (lp == 0)
      actpoint = GetEdge().BeginNode()->point();
    else if(lp == 1)
      actpoint = GetEdge().EndNode()->point();
    else
      dserror("You expect more than two Endpoints for a line?");

    switch (sidetype)
    {
      // --------------------------------------------------------------------
      /* handle parallel case for TRI3 elements in the 3-dimensional space
       * handle parallel case for 1-D elements in the 2-dimensional space */
      // --------------------------------------------------------------------
      case DRT::Element::tri3:
      case DRT::Element::line2:
      {
        lineendpoint_conv(lp) = ComputeDistance( actpoint, lineendpoint_dist(lp),
            lineendpoint_tol(lp), zeroarea , true ) ;
        lineendpoint_within_surfacelimits(lp) = SurfaceWithinLimits();
        // copy xsi_ in the lp column of lineendpoint_xsi
        std::copy(xsi_.A(),xsi_.A()+(dimside+dimedge),&lineendpoint_xsi(0,lp));

        break;
      } // end: case (DRT::Element::tri3) and 1-D elements
      // --------------------------------------------------------------------
      // handle parallel case for QUAD4 elements
      // --------------------------------------------------------------------
      case DRT::Element::quad4:
      {
        bool done = false;
        /* problem here is that ComputeDistance might change the normal direction if
         * the projected point lies outside the element, that why we prefer to compute
         * the distance onto tri3! */
        LINALG::Matrix<2,1>       tri_dist;
        LINALG::Matrix<2,1>       tri_tol;
        LINALG::TMatrix<bool,2,1> tri_conv;
        unsigned tri = 0;
        for ( tri = 0; tri < 2; ++tri )
        {
          tri_conv( tri ) = ComputeDistance( actpoint, tri_dist( tri ), tri_tol( tri ),
              zeroarea , true , tri );

          if (debug)
          {
            //just warping check!!!
            if (tri == 1)
            {
              if ( ! (std::abs( tri_dist(0) - tri_dist(1) ) < tri_tol(0) or
                      std::abs( tri_dist(0) - tri_dist(1) ) < tri_tol(1) ) )
                dserror("Expect Distance from both Triangles to be the same --> "
                    "is your QUAD4 warped??");
            }
          }

          /* In this case the point is not on the surface and we stop calculation of
           * distance right now! (could also be the second triangle if first one was
           * distorted --> very big tolerance) */
          if ( std::abs( tri_dist( tri ) ) >= tri_tol( tri ) )
          {
            /* REMARK: lineendpoint_within_surfacelimits[lp] is not correctly evaluated
             * in this case, as it is not needed!!! */
            done = true;
            break;
          }
          else
          {
            lineendpoint_within_surfacelimits(lp) = Tri3WithinLimits(tri_tol(tri));
            if (lineendpoint_within_surfacelimits(lp))
            {
              if (tri_conv(tri) == false)
                dserror("You found a case where the triangulation does not converge "
                    "and it is inside the element, basically removing this dserror "
                    "should be ok, but please contact me so that I can have a look into "
                    "it!!! (Christoph Ager)");
              break;
            }
          }
        }
        //calculate distance to triangles was enough
        if (done or tri == 2)
        {
          if (tri == 2)
          {
            //in this case both tri3s where outside!
            lineendpoint_within_surfacelimits(lp) = false;
            //set back to the last tri
            tri = 1;
          }
          lineendpoint_conv(lp) = tri_conv(tri);
          lineendpoint_dist(lp) = tri_dist(tri);
          lineendpoint_tol(lp) = tri_tol(tri);
          if (debug)
            std::copy(xsi_.A(),xsi_.A()+(dimside+dimedge),&lineendpoint_xsi(0,lp));
        }
        /* calculate distance to QUAD4, as the projected points lies inside and we
         * want to calculate local coords (should converge in 1 Newton step (theoretically
         * also tri3 didn't converge, but then there is anyway no hope) */
        else
        {
         // std::cout << "calculate distance from quad4 ..." << std::endl;
          lineendpoint_conv(lp) = ComputeDistance( actpoint, lineendpoint_dist(lp),
              lineendpoint_tol(lp), zeroarea , true ) ;
          lineendpoint_within_surfacelimits(lp) = SurfaceWithinLimits();
          std::copy(xsi_.A(),xsi_.A()+(dimside+dimedge),&lineendpoint_xsi(0,lp));
        }
        if (zeroarea)
          std::runtime_error("zeroarea handling for quad4 not implemented yet (not "
              "expected to happen here, as it should just occur in the triangulation "
              "procedure!)");

        break;
      } // end: case (DRT::Element::quad4)
      default:
      {
        dserror("Intersection Error: Other surfaces than TRI3 and QUAD4 are not supported "
            "to avoid huge problems!!!");
        break;
      }
    }
  }

  /* both points outside the plane on the same side --> no intersection
   * Note: distances of both end-points have the same sign. */
  if (std::abs(lineendpoint_dist(0)) >= lineendpoint_tol(0) and
      std::abs(lineendpoint_dist(1)) >= lineendpoint_tol(1) and
      lineendpoint_dist(0)*lineendpoint_dist(1)>0 and
      not zeroarea)
  {
    if (debug)
    {
      std::cout << "No Cut Point, because we detected following normal distance "
          "at the line end points: ... lineendpoint_dist[0]: " << lineendpoint_dist(0) <<
          " lineendpoint_dist[1]: " << lineendpoint_dist(1) << "(Codeline:" << __LINE__ <<
          ")" << std::endl;
    }
    return false; //there is no intersection between surface and line!!!
  }
  //for all paralell cases this should converge!!!
  // ------------------------------------------------------------------------
  // starting- and end-point of the edge are inside/on the side
  // ------------------------------------------------------------------------
  if (lineendpoint_conv(0) and lineendpoint_conv(1))
  {
    if ( std::abs(lineendpoint_dist(0)) < lineendpoint_tol(0) and
        lineendpoint_within_surfacelimits(0) )
      lineendpoint_in_surface(0) = true;

    if ( std::abs(lineendpoint_dist(1)) < lineendpoint_tol(1) and
        lineendpoint_within_surfacelimits(1) )
      lineendpoint_in_surface(1) = true;

    // both nodes of the line are inside/on the side
    if ( lineendpoint_in_surface(0) and
         lineendpoint_in_surface(1) )
    {
      InsertCut( GetEdge().BeginNode(), cuts );
      InsertCut( GetEdge().EndNode()  , cuts );
      TestSideEdges( GetEdge().BeginNode()->point(), &lineendpoint_xsi(0,0) );
      TestSideEdges( GetEdge().EndNode()  ->point(), &lineendpoint_xsi(0,1) );
      if (debug)
      {
        std::cout << "Cut points found, Begin & End Node of the Line are inside "
            "the surface: ... lineendpoint_dist[0]: " << lineendpoint_dist(0) <<
            " lineendpoint_dist[1]: " << lineendpoint_dist(1) << " (Codeline:" <<
            __LINE__ << ")" << std::endl;
      }
      if (debugbb)
        dserror("Bounding Boxes said that there is no intersection point!");
      return true;
    }
    // ----------------------------------------------------------------------
    // only the starting point of the edge is inside/on the side
    // ----------------------------------------------------------------------
    else if (lineendpoint_in_surface(0))
    {
      InsertCut( GetEdge().BeginNode(), cuts );
      nodalpoint = GetEdge().BeginNode()->point();
      TestSideEdges( nodalpoint, &lineendpoint_xsi(0,0));
      if (debug)
      {
        std::cout << "Cut points found, Begin Node of the Line are inside the surface: ... "
            "lineendpoint_dist[0]: " << lineendpoint_dist(0) << " lineendpoint_dist[1]: "
            << lineendpoint_dist(1) << " (Codeline:" << __LINE__ << ")" << std::endl;
      }
      success = true;
    }
    // ----------------------------------------------------------------------
    // only the end-point of the edge is inside/on the side
    // ----------------------------------------------------------------------
    else if (lineendpoint_in_surface(1))
    {
      InsertCut( GetEdge().EndNode(), cuts );
      nodalpoint = GetEdge().EndNode()->point();
      TestSideEdges( nodalpoint, &lineendpoint_xsi(0,1));
      if (debug)
      {
        std::cout << "Cut points found, End Node of the Line are inside the surface: ... "
            "lineendpoint_dist[0]: " << lineendpoint_dist(0) << " lineendpoint_dist[1]: "
            << lineendpoint_dist(1) << " (Codeline:" << __LINE__ << ")" << std::endl;
      }
      success = true;
    }
  }

  // ------------------------------------------------------------------------
  /* (3) see if edges of the given side need a special treatment
   *     --> intersection of side edges with given edge */
  // ------------------------------------------------------------------------
  // Search all side edges for a cut. We could cross the side.
  const std::vector<Edge*> & side_edges = GetSide().Edges();
  for ( std::vector<Edge*>::const_iterator i=side_edges.begin();
      i!=side_edges.end(); ++i )
  {
    Edge * e = *i;

    if ( nodalpoint!=NULL and nodalpoint->IsCut( e ) )
    {
      /* No need to do anything, since the side has a closed cycle of
       * edges. Any matching nodes will be found. */
      success = true;
      continue;
    }

    PointSet cut_points;
    GetEdge().GetCutPoints( e, cut_points );
    if ( cut_points.size() > 0 )
    {
      /* Nothing to be done. There are already cut points between these
       * edges. We cannot find new ones. */
      if (cut_points.size() == 1)
        cuts.insert(*(cut_points.begin()));
      else if (cut_points.size() == 2)
      {
        PointSet::iterator cp = cut_points.begin();
        cuts.insert(*cp);
        cp++;
        cuts.insert(*cp);
      }
      else
      {
        if (debug)
        {
          for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
            (*i)->Print(std::cout);
          GetEdge().BeginNode()->point()->Print(std::cout);
          GetEdge().EndNode()->point()->Print(std::cout);
          e->BeginNode()->point()->Print(std::cout);
          e->EndNode()->point()->Print(std::cout);
        }
        throw std::runtime_error("Two Edges have more than two cutpoint, sounds strange!");
      }
      if (debug)
      {
        std::cout << "Cut points found, by intersection of edges!"
            << " (Codeline:" << __LINE__ << ")" << std::endl;
      }
      success = true;
    } // if ( cut_points.size() > 0 )
    // given edge has currently no cut points
    else
    {
      double tolerance = 0.0;
      // check if the given edge intersects with one of the side edges
      if ( ComputeCut( e, GetEdgePtr(), cuts, tolerance ) )
      {
        if (debug)
        {
          std::cout << "Cut points found, by intersection of edges (edges of side "
              "and given edge)!" << " (Codeline:" << __LINE__ << ")" << std::endl;
        }
        success = true;
      }
    }
  } // end: loop over side edges

  //parallel cases handled!!!
  if (success)
  {
    if (debugbb)
      dserror("Bounding Boxes said that there is no intersection point!");
    return true;
  }
  else if ( lineendpoint_conv(0) and lineendpoint_conv(1) )
  {
    // parallel to the plane but No Cut Points found!!!
    if (std::abs(lineendpoint_dist(0)) < lineendpoint_tol(0) and
        std::abs(lineendpoint_dist(1)) < lineendpoint_tol(1))
    {
      if (debug)
      {
        std::cout << "No Cut Point, because we detected following normal distance at"
            " the line end points: ... lineendpoint_dist[0]: " << lineendpoint_dist(0)
            << " lineendpoint_dist[1]: " << lineendpoint_dist(1)
            << " and there was no edge intersection! (Codeline:" << __LINE__ << ")"
            << std::endl;
      }
      return false;
    }
  }
  else if (zeroarea)
  {
    if (debug)
    {
      std::cout << "No Cut Point, because surface has no area and edges don't have "
          "any intersection! (Codeline:" << __LINE__ << ")" << std::endl;
    }
    return false;
  }

  // ------------------------------------------------------------------------
  // (4) try to calculate the intersection point directly with Newton
  // ------------------------------------------------------------------------
  double itol = 0.0;
  bool conv = ComputeEdgeSideIntersection(itol,false);
  if (conv) //if newton converges we trust the result!
  {
    if ( SurfaceWithinLimits() and LineWithinLimits() )
    {
      FinalPoint();
      //case point lies on begin node of line!!!
      if ( std::abs( xsi_( dimside, 0 )+1.0 ) < REFERENCETOL )
      {
        InsertCut( GetEdge().BeginNode(), cuts );
        TestSideEdges( GetEdge().BeginNode()->point(), xsi_ );
        if (debug)
          /* Basically this case should already be caught in section (2)
           * for paralell intersections! */
          std::cout << "Cut points found by intersection, "
              "Begin Node is in surface! (Codeline:" << __LINE__ << ")"
              << std::endl;
        if (debugbb)
          dserror("Bounding Boxes said that ther is no intersection point!");
        return true;
      }
      //case point lies on end node of line!!!
      else if ( std::abs( xsi_( dimside, 0 )-1.0 ) < REFERENCETOL )
      {
        InsertCut( GetEdge().EndNode(), cuts );
        TestSideEdges( GetEdge().EndNode()->point(), xsi_ );
        if (debug)
        {
          /* Basically this case should already be caught in section
           * (2) for paralell intersections! */
          std::cout << "Cut points found by intersection, End Node is in surface! "
              "(Codeline:" << __LINE__ << ")" << std::endl;
        }
        return true;
      }
      else
      {
        Node * n = GetSide().OnNode( x_ ); //point lies on a node of the side???
        if ( n!=NULL ) //n is a new point
        {
          InsertCut( n, cuts );
          TestSideEdges( n->point(), xsi_ );
          if (debug)
          {
            /* Basically this case should already be caught in section (2)
             * for parallel intersections! */
            std::cout << "Cut points found by intersection, with local coords: xsi: "
                << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0)
                << " (Codeline:" << __LINE__ << ")" << std::endl;
          }
          if (debugbb)
            dserror("Bounding Boxes said that there is no intersection point!");
          return true;
        }
        else
        {
          Point * p = Point::NewPoint( GetMesh(), x_.A(), xsi_( dimside, 0 ),
              GetEdgePtr(), GetSidePtr(), itol );
          TestSideEdges( p, xsi_ );
          cuts.insert( p );
          if (debug)
          {
            std::cout << "Cut points found by intersection, with local coords: xsi: "
                << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0)
                << " (Codeline:" << __LINE__ << ")" << std::endl;
          }
          if (debugbb)
            dserror("Bounding Boxes said that there is no intersection point!");
          return true;
        }
      }
    }
    else
    {
      if (debug)
      {
        std::cout << "No Cut Point found by intersection, with local coords: xsi: "
            << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0)
            << " and a tolerance of " << itol << " (Codeline:" << __LINE__ << ")"
            << std::endl;
      }
      //there is no intersection between surface and line!!!
      return false;
    }
  }
  // ------------------------------------------------------------------------
  /* (5) treat special cases ...
   *     (no parallel case and no standard case) */
  // ------------------------------------------------------------------------
  else
  {
    i_xsi = xsi_;

    /* (5.1) point is outside the element and is not part of the interpolation
     * space (just for QUAD4 try triangulation) */
    if (sidetype == DRT::Element::quad4)
    {
      LINALG::Matrix<2,1> tri_tol;
      LINALG::Matrix<2,1> tri_conv;
      for (unsigned tri = 0; tri < 2; ++tri)
      {
        tri_conv( tri ) = ComputeEdgeTri3Intersection( tri_tol( tri ), tri );
        /* we expect the QUAD4 to converge in case that the projected point is
         * inside the QUAD4 or TRI3! */
        if ( tri_conv( tri ) and Tri3WithinLimits( tri_tol( tri ) ) )
        {
          std::cout << "Local Coordinates from ComputeDistance "
              "(lineendpoint_dist[0] = " << lineendpoint_dist(0)
              << "[tol=" << lineendpoint_tol(0) << "]) on surface: "
              << "xsi: " << lineendpoint_xsi(0,0) << " / eta: "
              << lineendpoint_xsi(1,0) << " / alpha: "
              << lineendpoint_xsi(2,0) << std::endl;

          std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[1] = "
              << lineendpoint_dist(1) << "[tol=" << lineendpoint_tol(1)
              << "]) on surface: " << "xsi: " << lineendpoint_xsi(0,1)
              << " / eta: "<< lineendpoint_xsi(1,1) << " / alpha: "
              << lineendpoint_xsi(2,1) << std::endl;

          std::cout << "Local Coordinates from ComputeIntersection (tol="
              << itol << ") on Quad4: " << "xsi: " << i_xsi(0,0)
              << " / eta: "<< i_xsi(1,0) << " / alpha: "<< i_xsi(2,0) << std::endl;
          std::cout << "Local Coordinates from ComputeIntersection (tol="
              << tri_tol(tri) << ") on Tri3 ( " << tri << " ) : " << "xsi: " << xsi_(0,0)
              << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << std::endl;
          dserror("ComputeEdgeSideIntersection for Quad4 didn't converge, "
              "but ComputeEdgeTri3Intersection for triangulation (id=%d) is inside the Element!",
              tri);
        }
      }
      // both converged and are outside!!! ( tested directly )
      if (tri_conv(0) and tri_conv(1))
      {
        if (debug)
        {
          std::cout << "No Cut Point found by intersection with triangulated quad4! "
              "(Codeline:" << __LINE__ << ")" << std::endl;
        }
        // point outside the interpolation space of QUAD4
        return false;
      }
      //else
        /* Looks like this is the result of an distorted element, here triangulation is not
         * better that the distance directly from the QUAD4 */
    }

    // (5.2) distorted elements
    /* now let's do what is possible for the rest, basically we can just find
     * out if there is no intersection, everything else would be pure speculation! */

    {
      //scale this tolerance into reference system!!!

      //calculate the tolerance in the reference coordinate system!!!
      double scaling;
      double ritol;
      {
        /* not definitely sure if this is a choise, but at least the same
         * as in the cut kernel! */
        scaling = xyze_surfaceElement_.NormInf();
        double linescale = xyze_lineElement_.NormInf();
        if (linescale > scaling)
          scaling = linescale;
      }
      ritol = itol/scaling;

      xsi_ = i_xsi; //set xsi_ back to the intersection xsi!!!
      if ( SurfaceWithinLimits(ritol) and LineWithinLimits(ritol) )
      {
        /* we think there shouldn't be an intersection point ...
         * to prove this we refine our bounding boxes!!! */
        if (sidetype == DRT::Element::tri3)
        {
          if (!RefinedBBOverlapCheck())
          {
            if (debug)
            {
              std::cout << "No Cut Point found by Refined Bounding Boxes! "
                  "(Codeline:" << __LINE__ << ")" << std::endl;
            }
            return false;
          }
          else
          {
            //Output of the available information!!!
            std::cout << "local coordinates: <plane: " << xsi_(0) << ", "
                << xsi_(1) << ">, <line: " << xsi_(2) << "> | relative tolerance: "
                << ritol << std::endl;

            std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[0] = "
                << lineendpoint_dist(0) << "[tol=" << lineendpoint_tol(0) << "], conv = "
                << lineendpoint_conv(0) << ") on surface: " << "xsi: "
                << lineendpoint_xsi(0,0) << " / eta: "<< lineendpoint_xsi(1,0)
                << " / alpha: "<< lineendpoint_xsi(2,0) << std::endl;

            std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[1] = "
                << lineendpoint_dist(1) << "[tol=" << lineendpoint_tol(1)
                << "], conv = " << lineendpoint_conv(1) << ") on surface: " << "xsi: "
                << lineendpoint_xsi(0,1) << " / eta: "<< lineendpoint_xsi(1,1)
                << " / alpha: "<< lineendpoint_xsi(2,1) << std::endl;
            {
              //just to create gmsh output of the failed intersection!!!
              KERNEL::ComputeIntersection<probdim,edgetype,sidetype> ci( xsi_ );
              ci( xyze_surfaceElement_, xyze_lineElement_ );

              std::string filename(OUTPUT::GenerateGmshOutputFilename(".intersection_CUTFAIL.pos"));
              std::ofstream file(filename.c_str());
              ci.WritetoGmsh(file);
              file.close();
            }
            /* (6) throw error! - It is not allowed to remove this error
             * (don't even do that in your local version, as the cut then just
             * generates something undefined)!!! */
            throw std::runtime_error("CRITICAL ERROR - GEO::CUT::Intersection::Intersect: "
                "Your Intersection seems to be a special case, which is not treated right "
                "yet!!! --> to fix this the intersection code in BACI has to be improved!!!");
          }
        }
      }
      else
      {
        if (debug)
        {
          std::cout << "No Cut Point found by intersection, with local coords: xsi: "
          << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0)
          << " and a tolerance of " << itol << " (Codeline:" << __LINE__ << ")" << std::endl;
        }
        return false;
      }
    }
  }

  /* (6) throw error! - It is not allowed to remove this error
   * (don't even do that in your local version, as the cut then just
   * generates something undefined)!!! */
  throw std::runtime_error("CRITICAL ERROR - GEO::CUT::Intersection::Intersect: "
      "Your Intersection seems to be a special case, which is not treated right yet!!! "
      "--> to fix this the intersection code in BACI has to be improved!!!");
  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim,
          DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          bool debug,
          unsigned dimedge, unsigned dimside,
          unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim,edgetype,sidetype,debug,dimedge,dimside,
     numNodesEdge,numNodesSide>::RefinedBBOverlapCheck(int maxstep)
{
  if (sidetype != DRT::Element::tri3)
    dserror("RefinedBBOverlapCheck is made for distored tri3s!");

  std::vector< LINALG::Matrix<3, 1> > surfpoints;
  std::vector< LINALG::Matrix<3, 1> > linepoints;
  LINALG::Matrix<3, 1> actpoint;
  for (unsigned nid = 0; nid < numNodesSide; nid ++)
  {
    (GetSide().Nodes()[nid])->Coordinates(actpoint.A());
    surfpoints.push_back(actpoint);
  }
  //find shortest edge on the surface!
  int smallestedge;
  {
    LINALG::Matrix<3, 1> d1, d2, d3;
    d1.Update(1.0, surfpoints[1],-1.0,surfpoints[2]);
    d2.Update(1.0, surfpoints[2],-1.0,surfpoints[0]);
    d3.Update(1.0, surfpoints[0],-1.0,surfpoints[1]);
    std::vector<double> lengths;
    lengths.push_back(d1.Norm2());
    lengths.push_back(d2.Norm2());
    lengths.push_back(d3.Norm2());
    if (lengths[0] < lengths[1] && lengths[0] < lengths[2]) //0 is smallest
      smallestedge = 0;
    else if(lengths[1] < lengths[2] && lengths[1] < lengths[0]) //1 is smallest
      smallestedge = 1;
    else //2 is smallest plus all the rest
      smallestedge = 2;
  }

  LINALG::Matrix<3, 1> v1(true);
  LINALG::Matrix<3, 1> v2(true);
  v1.Update(1.0, surfpoints[(smallestedge+1)%3],-1.0,surfpoints[smallestedge]);
  v2.Update(1.0, surfpoints[(smallestedge+2)%3],-1.0,surfpoints[smallestedge]);

  for (unsigned nid = 0; nid < numNodesEdge; nid ++)
  {
    (GetEdge().Nodes()[nid])->Coordinates(actpoint.A());
    linepoints.push_back(actpoint);
  }
  LINALG::Matrix<3, 1> v3(true);
  v3.Update(1.0, linepoints[1],-1.0,linepoints[0]);

  bool overlap = false; //we have an overlap of the refined bounding boxes?

  int pow_fac = 4;
  int max_steps = pow(pow_fac,maxstep-1);
  double dmax_steps = 1.0/max_steps;
  int act_steps = 0;
  uint act_boxidx = 0;
  std::vector< std::vector< int > > overlappingidx;
  std::vector< std::vector< int > > newoverlappingidx;
  std::vector< int > tmpoverlappingidx;
  tmpoverlappingidx.push_back(0); //min idx with overlap for tri3
  tmpoverlappingidx.push_back(max_steps); //max idx with overlap tri3
  tmpoverlappingidx.push_back(0); //min idx with overlap for line2
  tmpoverlappingidx.push_back(max_steps); //max idx with overlap for line2
  overlappingidx.push_back(tmpoverlappingidx);

  LINALG::Matrix<3, 1> p1,p2,p3,p4,lp1,lp2;

  for (int refinestep = 1; refinestep <= maxstep; refinestep++)
  {
    newoverlappingidx.clear();
    act_boxidx = 0;
    act_steps = pow(pow_fac,(maxstep - refinestep));

    overlap = false;
    for (int surfstep = 0; surfstep < max_steps; surfstep += act_steps)
    {
      tmpoverlappingidx.clear();
      while ((act_boxidx + 1) < overlappingidx.size() and (overlappingidx[act_boxidx + 1])[0] < surfstep) act_boxidx++;
      if (!(overlappingidx[act_boxidx][0] <= surfstep and overlappingidx[act_boxidx][1] >= surfstep + act_steps )) continue;
      double alpha = surfstep*dmax_steps;
      double alphap = (surfstep+act_steps)*dmax_steps;
      p1.Update(1.0, surfpoints[smallestedge], alpha,v1);
      p2.Update(1.0, surfpoints[smallestedge], alphap,v1);
      p3.Update(1.0, surfpoints[smallestedge], alpha,v2);
      p4.Update(1.0, surfpoints[smallestedge], alphap,v2);
      Teuchos::RCP<BoundingBox> sbb = Teuchos::rcp( BoundingBox::Create() );
      sbb->AddPoint(p1);
      sbb->AddPoint(p2);
      sbb->AddPoint(p3);
      sbb->AddPoint(p4);
      for (int linestep = 0; linestep < max_steps; linestep += act_steps)
      {
        if (!(overlappingidx[act_boxidx][2] <= linestep and overlappingidx[act_boxidx][3] >= linestep + act_steps )) continue;
        double lalpha = linestep*dmax_steps;
        double lalphap = (linestep+act_steps)*dmax_steps;
        lp1.Update(1.0, linepoints[0], lalpha,v3);
        lp2.Update(1.0, linepoints[0], lalphap,v3);
        Teuchos::RCP<BoundingBox> ebb = Teuchos::rcp( BoundingBox::Create() );
        ebb->AddPoint(lp1);
        ebb->AddPoint(lp2);

        if (sbb->Within(POSITIONTOL/BOXOVERLAP,*ebb))
        {
          overlap = true;
          if (!tmpoverlappingidx.size())
          {
            tmpoverlappingidx.push_back(surfstep);//min idx with overlap for tri3
            tmpoverlappingidx.push_back(surfstep+act_steps);//max idx with overlap for tri3
            tmpoverlappingidx.push_back(linestep);//min idx with overlap for line2
            tmpoverlappingidx.push_back(linestep+act_steps);//max idx with overlap for line2
          }
          else
          {
            if (tmpoverlappingidx[2] > linestep) tmpoverlappingidx[2] = linestep;
            if (tmpoverlappingidx[3] < linestep + act_steps) tmpoverlappingidx[3] = linestep + act_steps;
          }
        }
      }
      if (tmpoverlappingidx.size()) newoverlappingidx.push_back(tmpoverlappingidx);
    }
    overlappingidx = newoverlappingidx;

    std::cout << "RefinedBBOverlapCheck: Refinement Level " << refinestep << " there is " << overlap << " overlap!" << std::endl;
    if (!overlap)
      break;
  }
  return overlap;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::IntersectionBase>
GEO::CUT::IntersectionFactory::CreateIntersection(
    DRT::Element::DiscretizationType edge_type,
    DRT::Element::DiscretizationType side_type ) const
{
  const int probdim = DRT::Problem::Instance()->NDim();
  switch ( edge_type )
  {
    case DRT::Element::line2:
      return Teuchos::rcp( CreateIntersection<DRT::Element::line2>(
          side_type, probdim ) );
    default:
      dserror("Unsupported edgeType! If meaningful, add your edgeType here. \n"
          "Given edgeType = %s", DRT::DistypeToString( edge_type ).c_str() );
      break;
  }
  exit( EXIT_FAILURE );
}



template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::line2>;
template class GEO::CUT::Intersection<3,DRT::Element::line2,DRT::Element::line2>;
//template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad4>;
template class GEO::CUT::Intersection<3,DRT::Element::line2,DRT::Element::quad4>;
//template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad8>;
template class GEO::CUT::Intersection<3,DRT::Element::line2,DRT::Element::quad8>;
//template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad9>;
template class GEO::CUT::Intersection<3,DRT::Element::line2,DRT::Element::quad9>;
//template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::tri3>;
template class GEO::CUT::Intersection<3,DRT::Element::line2,DRT::Element::tri3>;
