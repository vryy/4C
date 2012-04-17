#include "cut_kernel.H"
#include "cut_point.H"
#include "cut_position2d.H"

unsigned GEO::CUT::KERNEL::FindNextCornerPoint( const std::vector<Point*> & points,
                                                LINALG::Matrix<3,1> & x1,
                                                LINALG::Matrix<3,1> & x2,
                                                LINALG::Matrix<3,1> & x3,
                                                LINALG::Matrix<3,1> & b1,
                                                LINALG::Matrix<3,1> & b2,
                                                LINALG::Matrix<3,1> & b3,
                                                unsigned i )
{
  unsigned pointsize = points.size();
  unsigned j = ( i+1 ) % pointsize;
  if ( pointsize < 3 )
  {
    return j;
  }

  points[i]->Coordinates( x1.A() );
  points[j]->Coordinates( x2.A() );

  b1.Update( 1, x2, -1, x1, 0 );

  double norm = b1.Norm2();
  if ( norm < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  b1.Scale( 1./norm );

  if ( b1.Norm2() < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  i = j;
  for ( unsigned k=2; k<pointsize; ++k )
  {
    i = ( i+1 ) % pointsize;
    Point * p = points[i];
    p->Coordinates( x3.A() );

    b2.Update( 1, x3, -1, x1, 0 );

    norm = b2.Norm2();
    if ( norm < std::numeric_limits<double>::min() )
      throw std::runtime_error( "same point in facet not supported" );

    b2.Scale( 1./norm );

    // cross product to get the normal at the point
    b3( 0 ) = b1( 1 )*b2( 2 ) - b1( 2 )*b2( 1 );
    b3( 1 ) = b1( 2 )*b2( 0 ) - b1( 0 )*b2( 2 );
    b3( 2 ) = b1( 0 )*b2( 1 ) - b1( 1 )*b2( 0 );

    if ( b3.Norm2() > PLANARTOL )
    {
      // Found. Return last node on this line.
      return ( i+pointsize-1 ) % pointsize;
    }
  }

  // All on one line. Return first and last point.
  if ( j==0 )
  {
    return 0;
  }
  else
  {
    return pointsize-1;
  }
}

void GEO::CUT::KERNEL::FindCornerPoints( const std::vector<Point*> & points,
                                         std::vector<Point*> & corner_points )
{
  LINALG::Matrix<3,1> x1;
  LINALG::Matrix<3,1> x2;
  LINALG::Matrix<3,1> x3;
  LINALG::Matrix<3,1> b1;
  LINALG::Matrix<3,1> b2;
  LINALG::Matrix<3,1> b3;

  for ( unsigned i = FindNextCornerPoint( points, x1, x2, x3, b1, b2, b3, 0 );
        true;
        i = FindNextCornerPoint( points, x1, x2, x3, b1, b2, b3, i ) )
  {
    Point * p = points[i];
    if ( corner_points.size()>0 and corner_points.front()==p )
      break;
    corner_points.push_back( p );
  }
}

bool GEO::CUT::KERNEL::IsValidQuad4( const std::vector<Point*> & points )
{
  if ( points.size()==4 )
  {
    LINALG::Matrix<3,3> xyze;
    LINALG::Matrix<3,1> xyz;
    for ( int i=0; i<4; ++i )
    {
      points[( i+0 )%4]->Coordinates( &xyze( 0, 0 ) );
      points[( i+1 )%4]->Coordinates( &xyze( 0, 1 ) );
      points[( i+2 )%4]->Coordinates( &xyze( 0, 2 ) );
      points[( i+3 )%4]->Coordinates( &xyz( 0, 0 ) );

      Position2d<DRT::Element::tri3> pos( xyze, xyz );
      if ( pos.Compute() )
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

DRT::Element::DiscretizationType GEO::CUT::KERNEL::CalculateShape( const std::vector<Point*> & points,
                                                                   std::vector<Point*> & line_points )
{
  FindCornerPoints( points, line_points );

  if ( IsValidTri3( line_points ) )
  {
    return DRT::Element::tri3;
  }
  else if ( IsValidQuad4( line_points ) )
  {
    return DRT::Element::quad4;
  }

  return DRT::Element::dis_none;
}

/*-----------------------------------------------------------------------------------------------------*
  Check whether three points are lying on the same line by checking whether the cross product is zero
                                                                                          Sudhakar 04/12
*------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsOnLine( Point* & pt1, Point* & pt2, Point* & pt3 )
{
  LINALG::Matrix<3,1> x1,x2,x3;
  LINALG::Matrix<3,1> pt1pt2,pt1pt3,cross;
  pt1->Coordinates(x1.A());
  pt2->Coordinates(x2.A());
  pt3->Coordinates(x3.A());

  pt1pt2.Update(1,x2,-1,x1,0);
  pt1pt3.Update(1,x3,-1,x1,0);

  cross(0,0) = pt1pt2(1,0)*pt1pt3(2,0)-pt1pt2(2,0)*pt1pt3(1,0);
  cross(1,0) = pt1pt2(0,0)*pt1pt3(2,0)-pt1pt2(2,0)*pt1pt3(0,0);
  cross(2,0) = pt1pt2(1,0)*pt1pt3(0,0)-pt1pt2(0,0)*pt1pt3(1,0);

  if(cross.NormInf()<1e-8)  //if the cross product is zero - on the same line
    return true;
  return false;
}

/*---------------------------------------------------------------------------------------------------------*
      Check whether the list of points given forms a convex polygon                         Sudhakar 04/12
      If any 3 points fall along the line, this will delete the middle point and return new pt
      Intially the polygon is projected into the given plane
*----------------------------------------------------------------------------------------------------------*/
std::vector<int> GEO::CUT::KERNEL::CheckConvexity( std::vector<Point*>& ptlist, std::string& geomType )
{
  if( ptlist.size()<5 )
    dserror( "The number of points < 5. Is it called for appropriate facet?" );

  std::string projPlane;

  for( unsigned i=0;i<ptlist.size();i++ )
  {
    Point* pt2 = ptlist[i];
    Point* pt3 = ptlist[(i+1)%ptlist.size()];
    unsigned ind = i-1;
    if( i==0 )
      ind = ptlist.size()-1;
    Point* pt1 = ptlist[ind];

    bool isline = IsOnLine( pt1, pt2, pt3 );
    if( isline==false )
    {
      std::vector<double> eqn;
      eqn = EqnPlane( pt1, pt2, pt3 );
      if( fabs(eqn[0])>fabs(eqn[1]) && fabs(eqn[0])>fabs(eqn[2]) )
        projPlane = "x";
      else if( fabs(eqn[1])>fabs(eqn[2]) && fabs(eqn[1])>fabs(eqn[0]) )
        projPlane = "y";
      else
        projPlane = "z";
      break;
    }
    else
      dserror( "Inline checking for facets not done before calling this" );
  }

  int ind1=0,ind2=0;
  if( projPlane=="x" )
  {
    ind1 = 1;
    ind2 = 2;
  }
  else if( projPlane=="y" )
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if( projPlane=="z" )
  {
    ind1 = 0;
    ind2 = 1;
  }
  else
    dserror( "unspecified projection type" );

  double x1[3],x2[3],x3[3],xtemp[3];
  std::vector<std::string> dirMove;
  for( unsigned i=0;i<ptlist.size();i++ )
  {
    Point* pt2 = ptlist[i];
    Point* pt3 = ptlist[(i+1)%ptlist.size()];
    unsigned ind = i-1;
    if(i==0)
      ind = ptlist.size()-1;
    Point* pt1 = ptlist[ind];

    pt1->Coordinates(x1);
    pt2->Coordinates(x2);
    pt3->Coordinates(x3);

    for( unsigned j=0;j<3;j++ )
      xtemp[j] = x2[j]-x1[j];

    double res = x3[ind1]*xtemp[ind2]-x3[ind2]*xtemp[ind1]+
                 xtemp[ind1]*x1[ind2]-xtemp[ind2]*x1[ind1];

    if(res<0.0)
      dirMove.push_back("left");
    else
      dirMove.push_back("right");
  }

  std::vector<int> leftind,rightind;
  for( unsigned i=0;i<dirMove.size();i++ )
  {
    if(dirMove[i]=="left")
      leftind.push_back(i);
    else
      rightind.push_back(i);
  }

  if( leftind.size()==0 )
  {
    geomType = "convex";
    return leftind;
  }
  else if( rightind.size()==0 )
  {
    geomType = "convex";
    return rightind;
  }
  else if( leftind.size()==1 )
  {
    geomType = "1ptConcave";
    return leftind;
  }
  else if( rightind.size()==1 )
  {
    geomType = "1ptConcave";
    return rightind;
  }
  else
  {
    geomType = "concave";
    if( leftind.size() < rightind.size() )
      return leftind;
    else
      return rightind;
  }
  return leftind;
}

/*-----------------------------------------------------------------------------------------------------*
            Find the equation of plane that contains these non-collinear points
                                                                                          Sudhakar 04/12
*------------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::KERNEL::EqnPlane( Point* & pt1, Point* & pt2, Point* & pt3 )
{
  std::vector<double> eqn_plane(4);
  double x1[3],x2[3],x3[3];

  pt1->Coordinates(x1);
  pt2->Coordinates(x2);
  pt3->Coordinates(x3);

  eqn_plane[0] = x1[1]*(x2[2]-x3[2])+x2[1]*(x3[2]-x1[2])+x3[1]*(x1[2]-x2[2]);
  eqn_plane[1] = x1[2]*(x2[0]-x3[0])+x2[2]*(x3[0]-x1[0])+x3[2]*(x1[0]-x2[0]);
  eqn_plane[2] = x1[0]*(x2[1]-x3[1])+x2[0]*(x3[1]-x1[1])+x3[0]*(x1[1]-x2[1]);
  eqn_plane[3] = x1[0]*(x2[1]*x3[2]-x3[1]*x2[2])+x2[0]*(x3[1]*x1[2]-x1[1]*x3[2])+x3[0]*(x1[1]*x2[2]-x2[1]*x1[2]);

  return eqn_plane;
}
