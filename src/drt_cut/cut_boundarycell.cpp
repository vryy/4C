
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_facet.H"


void GEO::CUT::BoundaryCell::CreateCells( Mesh & mesh,
                                          VolumeCell * cell,
                                          const std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  for ( std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> >::const_iterator i=sides_xyz.begin();
        i!=sides_xyz.end();
        ++i )
  {
    Facet * f = i->first;
    const std::vector<Epetra_SerialDenseMatrix> & xyz = i->second;
    switch ( xyz.size() )
    {
    case 0:
      throw std::runtime_error( "no surfaces" );
    case 1:
    {
      const Epetra_SerialDenseMatrix & x = xyz[0];
      if ( x.N()==3 )
      {
        cell->NewTri3Cell( mesh, f, x );
      }
      else if ( x.N()==4 )
      {
        cell->NewQuad4Cell( mesh, f, x );
      }
      else
      {
        throw std::runtime_error( "illegal number of point in boundary cell" );
      }
    }
    default:
      cell->NewTri3Cells( mesh, f, xyz );
    }
  }
}

void GEO::CUT::Tri3BoundaryCell::CollectCoordinates( const std::vector<Point*> & side,
                                                     std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  if ( side.size()!=3 )
  {
    throw std::runtime_error( "wrong number of points" );
  }

  std::set<Facet*> facets = side[0]->Facets();
  side[1]->Intersection( facets );
  side[2]->Intersection( facets );

  if ( facets.size()!=1 )
  {
    throw std::runtime_error( "not a valid cut side" );
  }

  Facet * f = *facets.begin();
  if ( f->SideId() < 0 )
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

  std::set<Facet*> facets = side[0]->Facets();
  side[1]->Intersection( facets );
  side[2]->Intersection( facets );
  side[3]->Intersection( facets );

  if ( facets.size()!=1 )
  {
    throw std::runtime_error( "not a valid cut side" );
  }

  Facet * f = *facets.begin();
  if ( f->SideId() < 0 )
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

double GEO::CUT::Tri3BoundaryCell::Area() const
{
  throw std::runtime_error( "" );
}

double GEO::CUT::Quad4BoundaryCell::Area() const
{
  throw std::runtime_error( "" );
}
