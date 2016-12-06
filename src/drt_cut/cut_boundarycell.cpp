/*---------------------------------------------------------------------*/
/*!
\file cut_boundarycell.cpp

\brief cut boundary cell

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_boundarycell.H"
#include "cut_volumecell.H"


#if 0
void GEO::CUT::Tri3BoundaryCell::CollectCoordinates( const std::vector<Point*> & side,
                                                     std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  if ( side.size()!=3 )
  {
    throw std::runtime_error( "wrong number of points" );
  }

  plain_facet_set facets;
  FindCommonFacets( side[0], side[1], side[2], facets );

  Facet * f;
  if ( facets.size()==1 )
  {
    f = *facets.begin();
  }
  else
  {
    Facet * found = NULL;
    for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->IsTriangle( side ) )
      {
        if ( found==NULL )
        {
          found = f;
        }
        else
        {
          throw std::runtime_error( "not unique" );
        }
      }
    }

    if ( found==NULL )
    {
      throw std::runtime_error( "not a valid cut side" );
    }

    f = found;
  }

  if ( not f->OnCutSide() )
  {
    return;
  }

  std::vector<Epetra_SerialDenseMatrix> & side_xyz = sides_xyz[f];
  side_xyz.push_back( Epetra_SerialDenseMatrix( 3, 3 ) );
  Epetra_SerialDenseMatrix & xyz = side_xyz.back();

  for ( int j=0; j<3; ++j )
  {
    Point * p = side[j];
    p->Coordinates( &xyz( 0, j ) );
  }
}

void GEO::CUT::Quad4BoundaryCell::CollectCoordinates( const std::vector<Point*> & side,
                                                      std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  if ( side.size()!=4 )
  {
    throw std::runtime_error( "wrong number of points" );
  }

  plain_facet_set facets;
  FindCommonFacets( side[0], side[1], side[2], side[3], facets );

  if ( facets.size()!=1 )
  {
    throw std::runtime_error( "not a valid cut side" );
  }

  Facet * f = *facets.begin();
  if ( not f->OnCutSide() )
  {
    return;
  }

  std::vector<Epetra_SerialDenseMatrix> & side_xyz = sides_xyz[f];
  side_xyz.push_back( Epetra_SerialDenseMatrix( 3, 4 ) );
  Epetra_SerialDenseMatrix & xyz = side_xyz.back();

  for ( int j=0; j<4; ++j )
  {
    Point * p = side[j];
    p->Coordinates( &xyz( 0, j ) );
  }
}
#endif

void GEO::CUT::Tri3BoundaryCell::DumpGmsh( std::ofstream & file, int * value )
{
  file.precision(16);

  file << "ST(";
  for ( int i=0; i<3; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<3; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << facet_->SideId();
  }
  file << "};\n";
}

void GEO::CUT::Quad4BoundaryCell::DumpGmsh( std::ofstream & file, int * value )
{
  file.precision(16);

  file << "SQ(";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << facet_->SideId();
  }
  file << "};\n";
}

void GEO::CUT::ArbitraryBoundaryCell::DumpGmsh( std::ofstream & file, int * value )
{
  //TO DO: implement gmsh output for arbitrarily shaped bcell
//  dserror("not implemented");
}

void GEO::CUT::Tri3BoundaryCell::DumpGmshNormal( std::ofstream & file )
{

  file.precision(16);

  file << "VP(";
  LINALG::Matrix<3,1> midpoint_triag(true);
  for ( int i=0; i<3; ++i )
  {
    LINALG::Matrix<3,1> cur(true);
    cur(0,0) = xyz_(0,i);
    cur(1,0) = xyz_(1,i);
    cur(2,0) = xyz_(2,i);
    midpoint_triag.Update(1.0,cur,1.0);
  }

  midpoint_triag.Scale(1.0/3.0);
  file << midpoint_triag(0,0) << ","
      << midpoint_triag(1,0) << ","
      << midpoint_triag(2,0);
  file << "){";

  //Choose midpoint of triangle as normal for now. Not best choice possibly.
  LINALG::Matrix<2,1> eta;
  eta(0,0) = 0.5;
  eta(1,0) = 0.5;
  LINALG::Matrix<3,1> normal;
  this->Normal(eta, normal);
  file << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );

  file << "};\n";
}

void GEO::CUT::Quad4BoundaryCell::DumpGmshNormal( std::ofstream & file )
{
  file.precision(16);

    file << "VP(";
    LINALG::Matrix<3,1> midpoint_quad(true);
    for ( int i=0; i<4; ++i )
    {
      LINALG::Matrix<3,1> cur(true);
      cur(0,0) = xyz_(0,i);
      cur(1,0) = xyz_(1,i);
      cur(2,0) = xyz_(2,i);
      midpoint_quad.Update(1.0,cur,1.0);
    }

    midpoint_quad.Scale(1.0/4.0);
    file << midpoint_quad(0,0) << ","
        << midpoint_quad(1,0) << ","
        << midpoint_quad(2,0);
    file << "){";

    //Choose midpoint of triangle as normal for now. Not best choice possibly.
    LINALG::Matrix<2,1> eta;
    eta(0,0) = 0.5;
    eta(1,0) = 0.5;
    LINALG::Matrix<3,1> normal;
    this->Normal(eta, normal);
    file << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );

    file << "};\n";
}

void GEO::CUT::ArbitraryBoundaryCell::DumpGmshNormal( std::ofstream & file )
{
  //TO DO: implement gmsh output for arbitrarily shaped bcell
//  dserror("not implemented");
}


void GEO::CUT::Tri3BoundaryCell::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const
{
#if 1
  // get derivatives at pos
  LINALG::Matrix<3,3> side_xyze( xyz_.A(), true );

  LINALG::Matrix<2,3> deriv;
  LINALG::Matrix<2,3> A;

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::tri3);
  A.MultiplyNT( deriv, side_xyze );

  // cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );

  double norm = normal.Norm2();
  normal.Scale( 1./norm );
#else
  LINALG::Matrix<3,1> x1( &xyz_( 0, 0 ) );
  LINALG::Matrix<3,1> x2( &xyz_( 0, 1 ) );
  LINALG::Matrix<3,1> x3( &xyz_( 0, 2 ) );

  LINALG::Matrix<3,1> r1( true);
  LINALG::Matrix<3,1> r2( true );

  r1.Update( 1., x1, -1., x3, 0.);
  r2.Update( 1., x2, -1., x1, 0.);

//  r1.Update( +1., x3, -1., x2, 0.);
//  r2.Update( +1., x1, -1., x3, 0.);

//  r1.Update( +1., x2, -1., x1, 0.);
//  r2.Update( +1., x3, -1., x2, 0.);


  r2.Scale(1.0/r2.Norm2());
  r1.Scale(1.0/r1.Norm2());

  // cross product to get the normal at the point
  normal( 0 ) = r1( 1 )*r2( 2 ) - r1( 2 )*r2( 1 );
  normal( 1 ) = r1( 2 )*r2( 0 ) - r1( 0 )*r2( 2 );
  normal( 2 ) = r1( 0 )*r2( 1 ) - r1( 1 )*r2( 0 );

  double norm = normal.Norm2();
  normal.Scale( 1./norm );
#endif
}

void GEO::CUT::Quad4BoundaryCell::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const
{
  // get derivatives at pos
  LINALG::Matrix<3,4> side_xyze( xyz_.A(), true );
  //Position2d<DRT::Element::quad4> position( side_xyze, xsi );
  //position.Normal( xsi, normal );

  LINALG::Matrix<2,4> deriv;
  LINALG::Matrix<2,3> A;

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::quad4);
  A.MultiplyNT( deriv, side_xyze );

  // cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );

  double norm = normal.Norm2();
  normal.Scale( 1./norm );
}

void GEO::CUT::ArbitraryBoundaryCell::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const
{


  dserror("Call GetNormalVector() to get normal for arbitrary boundarycells");
  exit(1);
  /*// cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );*/

//   double norm = normal.Norm2();
//   normal.Scale( 1./norm );
}

DRT::UTILS::GaussIntegration GEO::CUT::Tri3BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi( DRT::Element::tri3, CubatureDegree( DRT::Element::tri3 ) );
  return gi;
}

DRT::UTILS::GaussIntegration GEO::CUT::Quad4BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi( DRT::Element::quad4, CubatureDegree( DRT::Element::quad4 ) );
  return gi;
}

DRT::UTILS::GaussIntegration GEO::CUT::ArbitraryBoundaryCell::gaussRule()
{
  return gaussRule_;
}

void GEO::CUT::ArbitraryBoundaryCell::ElementCenter(LINALG::Matrix<3,1> & midpoint)
{
  dserror("Element Center for ArbitraryBoundaryCells not implemented!");
}


/*------------------------------------------------------------------------*
  | Center of Tri3 boundarycell                shahmiri 07/12
  *-----------------------------------------------------------------------*/
void GEO::CUT::Tri3BoundaryCell::ElementCenter( LINALG::Matrix<3,1> &  midpoint )
{
  LINALG::Matrix<3,1> center;
  center(0,0) = 0.25;
  center(1,0) = 0.25;
  center(2,0) = 0.0;
  MyElementCenter<DRT::Element::tri3>(center,midpoint);
}

/*------------------------------------------------------------------------*
  | Center of Quad4 boundarycell                shahmiri 07/12
  *-----------------------------------------------------------------------*/
void GEO::CUT::Quad4BoundaryCell::ElementCenter( LINALG::Matrix<3,1> & midpoint )
{
  LINALG::Matrix<3,1> center;
  center(0,0) = 0.0;
  center(1,0) = 0.0;
  center(2,0) = 0.0;
  MyElementCenter<DRT::Element::quad4>(center,midpoint);
}

LINALG::Matrix<3,1> GEO::CUT::Tri3BoundaryCell::GetNormalVector()
{
  dserror("Call Transform function to get normal for Tri3 boundarycell");
  LINALG::Matrix<3,1> normal;
  return normal;
}

LINALG::Matrix<3,1> GEO::CUT::Quad4BoundaryCell::GetNormalVector()
{
  dserror("Call Transform function to get normal for Quad4 boundarycell");
  LINALG::Matrix<3,1> normal;
  return normal;
}

LINALG::Matrix<3,1> GEO::CUT::ArbitraryBoundaryCell::GetNormalVector()
{
  return normal_;
}

/*------------------------------------------------------------------------*
 * Print end points of boundarycell (just for debugging)          sudhakar 10/14
 *------------------------------------------------------------------------*/
void GEO::CUT::BoundaryCell::Print()
{
  std::cout<<"The coordinates of boundary cell are \n";
  for( unsigned i=0; i<points_.size(); i++ )
  {
    std::cout<< "pID: " << points_()[i]->Id() << ", coord: " << xyz_(0,i)<<"\t"<<xyz_(1,i)<<"\t"<<xyz_(2,i)<<"\n";
  }
}

/*------------------------------------------------------------------------*
 * To make it easier with DD data structures (should be removed at some point)
 *                               winter 07/2015
 *------------------------------------------------------------------------*/
const std::vector<std::vector<double> > GEO::CUT::BoundaryCell::CoordinatesV()
{
  std::vector<std::vector<double> > corners;
  for(std::vector<Point*> ::const_iterator j=points_().begin();
      j!=points_().end(); ++j )
  {
    std::vector<double> cornerLocal(3);
    cornerLocal.assign((*j)->X(), (*j)->X() + 3);
    corners.push_back(cornerLocal);
  }
  return corners;
}

/*------------------------------------------------------------------------*
 *           A first step to validate if a boundary cell is valid.
 *                               winter 11/2015
 *------------------------------------------------------------------------*/
bool GEO::CUT::Tri3BoundaryCell::IsValidBoundaryCell()
{
  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

  const std::vector<Point*> points = this->Points();

  // create planes consisting of 3 nodes each
  LINALG::Matrix<numnodes,1> p0( points[0]->X() );
  LINALG::Matrix<numnodes,1> p1( points[1]->X() );
  LINALG::Matrix<numnodes,1> p2( points[2]->X() );

  LINALG::Matrix<numnodes,1> v01;
  LINALG::Matrix<numnodes,1> v02;
  LINALG::Matrix<numnodes,1> v12;

  v01.Update( 1, p1, -1, p0, 0 );
  v02.Update( 1, p2, -1, p0, 0 );
  v12.Update( 1, p1, -1, p2, 0 );

  //Get distance to origin
  LINALG::Matrix<numnodes,1> temp(true);
  temp(0,0)=p0.Norm2();         //Distance of points to origin
  temp(1,0)=p1.Norm2();
  temp(2,0)=p2.Norm2();

  // This is to scale the tolerance, it determines our maximum precision (i.e. machine precision)
//  double max_dist_to_orgin = temp.NormInf();

  //Distance between points in triangle
  temp(0,0)=v01.Norm2();
  temp(1,0)=v02.Norm2();
  temp(2,0)=v12.Norm2();

  double min_dist_in_tri = temp.MinValue();
  //We want to test with this one I think... But might lead to problems.
//  double tolerance = LINSOLVETOL*max_dist_to_orgin;

  temp(0,0)=points[0]->Tolerance();
  temp(1,0)=points[1]->Tolerance();
  temp(2,0)=points[2]->Tolerance();
  double max_tol_points = temp.MaxValue();


//  std::cout << "min_dist_in_tri: " << min_dist_in_tri << std::endl;
//  std::cout << "tolerance: " << tolerance << std::endl;
//  std::cout << "max_tol_points: " << max_tol_points << std::endl;

  if(min_dist_in_tri < max_tol_points)
  {
    return false;
  }

  return true;

}
