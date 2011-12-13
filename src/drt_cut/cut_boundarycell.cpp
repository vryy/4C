
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_facet.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


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

//   double norm = normal.Norm2();
//   normal.Scale( 1./norm );
#else
  LINALG::Matrix<3,1> x1( &xyz_( 0, 0 ) );
  LINALG::Matrix<3,1> x2( &xyz_( 0, 1 ) );
  LINALG::Matrix<3,1> x3( &xyz_( 0, 2 ) );

  x2.Update( -1., x1, 1. );
  x1.Update( -1., x3, 1. );

  // cross product to get the normal at the point
  normal( 0 ) = x1( 1 )*x2( 2 ) - x1( 2 )*x2( 1 );
  normal( 1 ) = x1( 2 )*x2( 0 ) - x1( 0 )*x2( 2 );
  normal( 2 ) = x1( 0 )*x2( 1 ) - x1( 1 )*x2( 0 );

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

//   double norm = normal.Norm2();
//   normal.Scale( 1./norm );
}

void GEO::CUT::ArbitraryBoundaryCell::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const
{


  dserror("Normal unavailable for this type of boundarycell");
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
