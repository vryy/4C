
#include "cut_position2d.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_tetgen.H"

#if 0
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
#endif


GEO::CUT::Point::PointPosition GEO::CUT::IntegrationCell::VolumePosition( const std::set<Facet*> & facets )
{
  Point::PointPosition pos = Point::undecided;
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    pos = std::min( pos, static_cast<Point::PointPosition>( f->PositionSideId() ) );
  }
  return pos;
}

GEO::CUT::Hex8IntegrationCell * GEO::CUT::Hex8IntegrationCell::CreateCell( Mesh & mesh,
                                                                           VolumeCell * cell,
                                                                           const std::set<Facet*> & facets )
{
  if ( facets.size()==6 )
  {
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Equals( DRT::Element::quad4 ) )
      {
        return NULL;
      }
    }

    // Must not be concave?!

    // So we have a hex8 here. Create it.
    //
    // Find two sides that do not share a common point, find positive
    // directions for both, and find the edges between those sides.

    std::set<Facet*>::const_iterator i=facets.begin();
    Facet * bot = *i;
    Facet * top = NULL;
    for ( ++i; i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Touches( bot ) )
      {
        top = f;
        break;
      }
    }
    if ( top==NULL )
    {
      throw std::runtime_error( "illegal hex8 cell" );
    }

    const std::vector<Point*> & bot_points = bot->Points();
    const std::vector<Point*> & top_points = top->Points();

    std::vector<Point*> points( 8, static_cast<Point*>( NULL ) );

    std::copy( bot_points.begin(), bot_points.end(), &points[0] );

    std::vector<Point*>::iterator end = points.begin()+4;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f!=bot and f!=top )
      {
        const std::vector<Point*> & side_points = f->Points();
        for ( std::vector<Point*>::const_iterator i=side_points.begin();
              i!=side_points.end();
              ++i )
        {
          Point * p = *i;
          std::vector<Point*>::iterator pointpos1 = std::find( points.begin(), end, p );
          if ( pointpos1!=end )
          {
            ++i;
            if ( i==side_points.end() )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            std::vector<Point*>::iterator pointpos2 = std::find( points.begin(), end, *i );
            if ( pointpos2==end )
            {
              i = side_points.end();
              pointpos2 = std::find( points.begin(), end, side_points[3] );
              if ( pointpos2==end )
              {
                throw std::runtime_error( "illegal hex8 cell" );
              }
              std::swap( pointpos1, pointpos2 );
            }

            unsigned pos = i - side_points.begin() - 1;
            Point * top_point2 = side_points[( pos + 2 ) % 4];
            Point * top_point1 = side_points[( pos + 3 ) % 4];

            if ( std::find( top_points.begin(), top_points.end(), top_point1 )==top_points.end() or
                 std::find( top_points.begin(), top_points.end(), top_point2 )==top_points.end() )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            pos = ( pointpos1 - points.begin() ) + 4;
            if ( points[pos]==NULL )
            {
              points[pos] = top_point1;
            }
            else if ( points[pos] != top_point1 )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            pos = ( pointpos2 - points.begin() ) + 4;
            if ( points[pos]==NULL )
            {
              points[pos] = top_point2;
            }
            else if ( points[pos] != top_point2 )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            break;
          }
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3,4> bot_xyze;
    LINALG::Matrix<3,4> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top->CornerCoordinates( top_xyze.A() );

    int distance_counter = 0;
    for ( int i=0; i<4; ++i )
    {
      Position2d<DRT::Element::quad4> bot_distance( bot_xyze, LINALG::Matrix<3,1>( top_xyze.A()+3*i, true ) );

      bot_distance.Compute();

      if ( bot_distance.Distance() > 0 )
      {
        distance_counter += 1;
      }
      if ( bot_distance.Distance() < 0 )
      {
        distance_counter -= 1;
      }
    }

    if ( distance_counter==-4 or
         distance_counter==4 )
    {
      for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        f->NewQuad4Cell( mesh );
      }
    }

    if ( distance_counter==-4 )
    {
      std::vector<Point*> rpoints( 8 );
      rpoints[0] = points[3];
      rpoints[1] = points[2];
      rpoints[2] = points[1];
      rpoints[3] = points[0];
      rpoints[4] = points[7];
      rpoints[5] = points[6];
      rpoints[6] = points[5];
      rpoints[7] = points[4];
      return mesh.NewHex8Cell( VolumePosition( facets ), rpoints, cell );
    }
    else if ( distance_counter==4 )
    {
      return mesh.NewHex8Cell( VolumePosition( facets ), points, cell );
    }
    else
    {
      return NULL;
    }
  }
  return NULL;
}

GEO::CUT::Tet4IntegrationCell * GEO::CUT::Tet4IntegrationCell::CreateCell( Mesh & mesh,
                                                                           VolumeCell * cell,
                                                                           const std::set<Facet*> & facets )
{
  if ( facets.size()==4 )
  {
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Equals( DRT::Element::tri3 ) )
      {
        return NULL;
      }
    }

    // So we have a tet4 here. Create it.

    std::set<Facet*>::const_iterator i=facets.begin();
    Facet * bot = *i;

    const std::vector<Point*> & bot_points = bot->Points();
    Point* top_point = NULL;

    for ( ++i; i!=facets.end(); ++i )
    {
      Facet * f = *i;

      const std::vector<Point*> & side_points = f->Points();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        if ( std::find( bot_points.begin(), bot_points.end(), p )==bot_points.end() )
        {
          if ( top_point==NULL )
          {
            top_point = p;
          }
          else if ( top_point!=p )
          {
            throw std::runtime_error( "illegal tet4 cell" );
          }
        }
      }
    }

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      f->NewTri3Cell( mesh );
    }

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,1> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top_point->Coordinates( top_xyze.A() );

    Position2d<DRT::Element::tri3> bot_distance( bot_xyze, top_xyze );
    bot_distance.Compute();

    if ( bot_distance.Distance() > 0 )
    {
      std::vector<Point*> points;
      points.reserve( 4 );
      std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
      points.push_back( top_point );
      return mesh.NewTet4Cell( VolumePosition( facets ), points, cell );
    }
    else if ( bot_distance.Distance() < 0 )
    {
      std::vector<Point*> points;
      points.reserve( 4 );
      std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
      points.push_back( top_point );
      return mesh.NewTet4Cell( VolumePosition( facets ), points, cell );
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return NULL;
}

void GEO::CUT::Tet4IntegrationCell::CreateCells( Mesh & mesh,
                                                 Element * element,
                                                 VolumeCell * cell,
                                                 const std::set<Facet*> & facets,
                                                 std::set<IntegrationCell*> & integrationcells )
{
#ifdef QHULL
  const int dim = 3;
  tetgenio in;

  std::set<Point*, PointPidLess> points;
  for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet & f = **i;
    f.GetPoints( mesh, points );
  }

  // allocate pointlist
  in.numberofpoints = points.size();
  in.pointlist = new double[in.numberofpoints * dim];

  int pos = 0;
  for ( std::set<Point*, PointPidLess>::iterator i=points.begin();
        i!=points.end();
        ++i )
  {
    Point & p = **i;
    p.Coordinates( & in.pointlist[pos*dim] );
//     for ( int j=0; j<dim; ++j )
//     {
//       in.pointlist[pos*dim+j] *= TETGENPOINTSCALE;
//     }
    pos += 1;
  }

  in.pointmarkerlist = new int[in.numberofpoints];
  std::fill( in.pointmarkerlist, in.pointmarkerlist+in.numberofpoints, 0 );

  in.numberoffacets = 0;
  for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet & facet = **i;
    in.numberoffacets += facet.NumTetgenFacets( mesh );
  }
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  std::fill( in.facetmarkerlist, in.facetmarkerlist+in.numberoffacets, 0 );

  std::vector<Point*> pointlist;
  pointlist.reserve( points.size() );
  std::copy( points.begin(), points.end(), std::back_inserter( pointlist ) );

  pos = 0;
  for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet & facet = **i;
    int numtets = facet.NumTetgenFacets( mesh );
    for ( int j=0; j<numtets; ++j )
    {
      tetgenio::facet & f = in.facetlist[pos];

      facet.GenerateTetgen( mesh, element, f, j, in.facetmarkerlist[pos], pointlist );

      pos += 1;
    }
  }

//   if ( generator==NULL )
//   {
//     // debug
//     const char * name = "tet";
//     in.save_nodes( const_cast<char*>( name ) );
//     in.save_poly( const_cast<char*>( name ) );
//   }

  char switches[] = "pQ";    // pQ o2 Y R nn
  tetgenio out;

  try
  {
    tetrahedralize( switches, &in, &out );
  }
  catch ( int & err )
  {
    const char * name = "tet";
    in.save_nodes( const_cast<char*>( name ) );
    in.save_poly( const_cast<char*>( name ) );

    throw;
  }

//   if ( generator!=NULL )
//   {
//     generator->Generate( parent, out );
//   }
//   else
//   {
//     // debug
//     // Output mesh
//     const char * name = "tetout";
//     out.save_nodes( const_cast<char*>( name ) );
//     out.save_elements( const_cast<char*>( name ) );
//     out.save_faces( const_cast<char*>( name ) );
//   }

  std::map<int, int> nodeidmap;

  std::map<int, std::vector<Epetra_SerialDenseMatrix> > sides_xyz;

  // need to copy surface markers via side id
  for ( int i=0; i<out.numberoftrifaces; ++i )
  {
    if ( out.trifacemarkerlist[i] > -1 )
    {
      std::vector<Epetra_SerialDenseMatrix> & side_xyz = sides_xyz[out.trifacemarkerlist[i]];
      side_xyz.push_back( Epetra_SerialDenseMatrix( 3, 3 ) );
      Epetra_SerialDenseMatrix & xyz = side_xyz.back();
      for ( int j=0; j<3; ++j )
      {
        int pointidx = out.trifacelist[i*3+j] * 3;
        std::copy( &out.pointlist[pointidx], &out.pointlist[pointidx+3], &xyz( 0, j ) );
      }
//       xyz.Scale( 1./TETGENPOINTSCALE );
    }
  }

  std::map<int, Facet*> sides;

  for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    int sid = f->SideId();
    if ( sid > -1 )
    {
      sides[sid] = f;
    }
  }

  for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::iterator i=sides_xyz.begin();
        i!=sides_xyz.end();
        ++i )
  {
    std::map<int, Facet*>::iterator j = sides.find( i->first );
    if ( j!=sides.end() )
    {
      Facet * f = j->second;
      std::vector<Epetra_SerialDenseMatrix> & xyz = i->second;
      f->NewTri3Cells( mesh, xyz );
    }
    else
    {
      throw std::runtime_error( "side id not found" );
    }
  }

  const int numTetNodes = 4;
  Epetra_SerialDenseMatrix xyz( 3, 4 );

  Point::PointPosition position = VolumePosition( facets );

  for ( int i=0; i<out.numberoftetrahedra; ++i )
  {
    for ( int j=0; j<numTetNodes; ++j )
    {
      int pointidx = out.tetrahedronlist[i*out.numberofcorners+j] * 3;
      std::copy( &out.pointlist[pointidx], &out.pointlist[pointidx+3], &xyz( 0, j ) );
    }
//     xyz.Scale( 1./TETGENPOINTSCALE );
    integrationcells.insert( mesh.NewTet4Cell( position, xyz, cell ) );
  }

#endif
}

GEO::CUT::Wedge6IntegrationCell * GEO::CUT::Wedge6IntegrationCell::CreateCell( Mesh & mesh,
                                                                               VolumeCell * cell,
                                                                               const std::set<Facet*> & facets )
{
  if ( facets.size()==5 )
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Equals( DRT::Element::tri3 ) )
      {
        tris.push_back( f );
      }
      else if ( f->Equals( DRT::Element::quad4 ) )
      {
        quads.push_back( f );
      }
      else
      {
        return NULL;
      }
    }

    if ( tris.size()!=2 or quads.size()!=3 )
    {
      return NULL;
    }

    Facet * bot = tris[0];
    Facet * top = tris[1];

    const std::vector<Point*> & bot_points = bot->Points();
    const std::vector<Point*> & top_points = top->Points();

    std::vector<Point*> points( 6, static_cast<Point*>( NULL ) );

    std::copy( bot_points.begin(), bot_points.end(), &points[0] );

    std::vector<Point*>::iterator end = points.begin()+3;

    for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
    {
      Facet * f = *i;
      const std::vector<Point*> & side_points = f->Points();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        std::vector<Point*>::iterator pointpos1 = std::find( points.begin(), end, p );
        if ( pointpos1!=end )
        {
          ++i;
          if ( i==side_points.end() )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          std::vector<Point*>::iterator pointpos2 = std::find( points.begin(), end, *i );
          if ( pointpos2==end )
          {
            i = side_points.end();
            pointpos2 = std::find( points.begin(), end, side_points[3] );
            if ( pointpos2==end )
            {
              throw std::runtime_error( "illegal wedge6 cell" );
            }
            std::swap( pointpos1, pointpos2 );
          }

          unsigned pos = i - side_points.begin() - 1;
          Point * top_point2 = side_points[( pos + 2 ) % 4];
          Point * top_point1 = side_points[( pos + 3 ) % 4];

          if ( std::find( top_points.begin(), top_points.end(), top_point1 )==top_points.end() or
               std::find( top_points.begin(), top_points.end(), top_point2 )==top_points.end() )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          pos = ( pointpos1 - points.begin() ) + 3;
          if ( points[pos]==NULL )
          {
            points[pos] = top_point1;
          }
          else if ( points[pos] != top_point1 )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          pos = ( pointpos2 - points.begin() ) + 3;
          if ( points[pos]==NULL )
          {
            points[pos] = top_point2;
          }
          else if ( points[pos] != top_point2 )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          break;
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,3> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top->CornerCoordinates( top_xyze.A() );

    int distance_counter = 0;
    for ( int i=0; i<3; ++i )
    {
      Position2d<DRT::Element::tri3> bot_distance( bot_xyze, LINALG::Matrix<3,1>( top_xyze.A()+3*i, true ) );

      bot_distance.Compute();

      if ( bot_distance.Distance() > 0 )
      {
        distance_counter += 1;
      }
      if ( bot_distance.Distance() < 0 )
      {
        distance_counter -= 1;
      }
    }

    if ( distance_counter==-3 or
         distance_counter==3 )
    {
      for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
      {
        Facet * f = *i;
        f->NewTri3Cell( mesh );
      }
      for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
      {
        Facet * f = *i;
        f->NewQuad4Cell( mesh );
      }
    }

    if ( distance_counter==-3 )
    {
      std::vector<Point*> rpoints( 6 );
      rpoints[0] = points[2];
      rpoints[1] = points[1];
      rpoints[2] = points[0];
      rpoints[3] = points[5];
      rpoints[4] = points[4];
      rpoints[5] = points[3];
      return mesh.NewWedge6Cell( VolumePosition( facets ), rpoints, cell );
    }
    else if ( distance_counter==3 )
    {
      return mesh.NewWedge6Cell( VolumePosition( facets ), points, cell );
    }
    else
    {
      return NULL;
    }
  }
  return NULL;
}

GEO::CUT::Pyramid5IntegrationCell * GEO::CUT::Pyramid5IntegrationCell::CreateCell( Mesh & mesh,
                                                                                   VolumeCell * cell,
                                                                                   const std::set<Facet*> & facets )
{
  if ( facets.size()==5 )
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Equals( DRT::Element::tri3 ) )
      {
        tris.push_back( f );
      }
      else if ( f->Equals( DRT::Element::quad4 ) )
      {
        quads.push_back( f );
      }
      else
      {
        return NULL;
      }
    }

    if ( tris.size()!=4 or quads.size()!=1 )
    {
      return NULL;
    }

    Facet * bot = quads[0];

    const std::vector<Point*> & bot_points = bot->Points();
    Point* top_point = NULL;

    for ( std::vector<Facet*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Facet * f = *i;

      const std::vector<Point*> & side_points = f->Points();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        if ( std::find( bot_points.begin(), bot_points.end(), p )==bot_points.end() )
        {
          if ( top_point==NULL )
          {
            top_point = p;
          }
          else if ( top_point!=p )
          {
            throw std::runtime_error( "illegal pyramid5 cell" );
          }
        }
      }
    }

    for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Facet * f = *i;
      f->NewTri3Cell( mesh );
    }
    for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
    {
      Facet * f = *i;
      f->NewQuad4Cell( mesh );
    }

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,1> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top_point->Coordinates( top_xyze.A() );

    Position2d<DRT::Element::tri3> bot_distance( bot_xyze, top_xyze );
    bot_distance.Compute();

    if ( bot_distance.Distance() > 0 )
    {
      std::vector<Point*> points;
      points.reserve( 5 );
      std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
      points.push_back( top_point );
      return mesh.NewPyramid5Cell( VolumePosition( facets ), points, cell );
    }
    else if ( bot_distance.Distance() < 0 )
    {
      std::vector<Point*> points;
      points.reserve( 5 );
      std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
      points.push_back( top_point );
      return mesh.NewPyramid5Cell( VolumePosition( facets ), points, cell );
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return NULL;
}

void GEO::CUT::Hex8IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SH(";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Tet4IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SS(";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Wedge6IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SI(";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Pyramid5IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SP(";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

// double GEO::CUT::Hex8IntegrationCell::Volume() const
// {
//   throw std::runtime_error( "" );
// }

// double GEO::CUT::Tet4IntegrationCell::Volume() const
// {
//   throw std::runtime_error( "" );
// }

// double GEO::CUT::Wedge6IntegrationCell::Volume() const
// {
//   throw std::runtime_error( "" );
// }

// double GEO::CUT::Pyramid5IntegrationCell::Volume() const
// {
//   throw std::runtime_error( "" );
// }
