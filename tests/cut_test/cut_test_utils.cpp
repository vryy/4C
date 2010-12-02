
#include "cut_test_utils.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_meshintersection.H"

int numnode;
int numele;

GEO::CUT::Element* create_hex8( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 8 );
  for ( int i=0; i<8; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateHex8( numele++, nids );
}

GEO::CUT::Element* create_hex20( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateHex20( numele++, nids );
}

GEO::CUT::Element* create_hex27( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 27 );
  for ( int i=0; i<27; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateHex27( numele++, nids );
}

GEO::CUT::Element* create_tet4( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 4 );
  for ( int i=0; i<4; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateTet4( numele++, nids );
}

GEO::CUT::Element* create_wedge6( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 6 );
  for ( int i=0; i<6; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateWedge6( numele++, nids );
}

GEO::CUT::Element* create_pyramid5( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 5 );
  for ( int i=0; i<5; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreatePyramid5( numele++, nids );
}

GEO::CUT::Side* create_quad4( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 4 );
  for ( int i=0; i<4; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateQuad4( numele++, nids );
}

GEO::CUT::Side* create_quad8( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 8 );
  for ( int i=0; i<8; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateQuad8( numele++, nids );
}

GEO::CUT::Side* create_quad9( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze )
{
  std::vector<int> nids;
  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    mesh.GetNode( numnode, &xyze( 0, i ) );
    nids.push_back( numnode++ );
  }

  return mesh.CreateQuad9( numele++, nids );
}

GEO::CUT::Element* create_hex8( GEO::CUT::Mesh & mesh, double dx, double dy, double dz )
{
  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 0;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0;
  xyze( 1, 3 ) = 1;
  xyze( 2, 3 ) = 0;

  xyze( 0, 4 ) = 0;
  xyze( 1, 4 ) = 0;
  xyze( 2, 4 ) = 1;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = 0;
  xyze( 2, 5 ) = 1;

  xyze( 0, 6 ) = 1;
  xyze( 1, 6 ) = 1;
  xyze( 2, 6 ) = 1;

  xyze( 0, 7 ) = 0;
  xyze( 1, 7 ) = 1;
  xyze( 2, 7 ) = 1;

  for ( int i=0; i<8; ++i )
  {
    xyze( 0, i ) += dx;
    xyze( 1, i ) += dy;
    xyze( 2, i ) += dz;
  }

  return create_hex8( mesh, xyze );
}

void create_hex8_mesh( GEO::CUT::Mesh & mesh, int rows, int cols, int depth )
{
  for (int i=0; i<rows+1; ++i)
  {
    for (int j=0; j<cols+1; ++j)
    {
      for ( int k=0; k<depth+1; ++k )
      {
        int id = i+j*(rows+1)+k*(rows+1)*(cols+1);
        double coord[3];

        coord[0] = 1./rows*i;
        coord[1] = 1./cols*j;
        coord[2] = 1./depth*k;

        mesh.GetNode( numnode+id, coord );
      }
    }
  }

  for (int i=0; i<rows; ++i)
  {
    for (int j=0; j<cols; ++j)
    {
      for ( int k=0; k<depth; ++k )
      {
        std::vector<int> nids;
        nids.reserve( 8 );
        nids.push_back( numnode+i+    j*(rows+1)  +k*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+    j*(rows+1)+1+k*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+(j+1)*(rows+1)+1+k*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+(j+1)*(rows+1)  +k*(rows+1)*(cols+1) );

        nids.push_back( numnode+i+    j*(rows+1)  +( k+1 )*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+    j*(rows+1)+1+( k+1 )*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+(j+1)*(rows+1)+1+( k+1 )*(rows+1)*(cols+1) );
        nids.push_back( numnode+i+(j+1)*(rows+1)  +( k+1 )*(rows+1)*(cols+1) );

        mesh.CreateHex8( numele++, nids );
      }
    }
  }

  numnode += (rows+1)*(cols+1)*(depth+1);
}

void create_quad4_mesh( GEO::CUT::Mesh & mesh, int rows, int cols, std::vector<GEO::CUT::Side*> & sides )
{
  double sqrt2 = 1./sqrt( 2. );

  for (int i=0; i<rows+1; ++i)
  {
    for (int j=0; j<cols+1; ++j)
    {
      int id = i+j*(rows+1);
      double coord[3];

      double x = ( 2./rows*i - 1 );
      double y = ( 2./cols*j - 1 );

      coord[0] = x*sqrt2 - y*sqrt2 + 0.5;
      coord[1] = x*sqrt2 + y*sqrt2 + 0.5;
      coord[2] = 0.5;

      mesh.GetNode( numnode+id, coord );
    }
  }

  for (int i=0; i<rows; ++i)
  {
    for (int j=0; j<cols; ++j)
    {
      std::vector<int> nids;
      nids.reserve( 4 );
      nids.push_back( numnode+i+    j*(rows+1)   );
      nids.push_back( numnode+i+    j*(rows+1)+1 );
      nids.push_back( numnode+i+(j+1)*(rows+1)+1 );
      nids.push_back( numnode+i+(j+1)*(rows+1)   );

      sides.push_back( mesh.CreateQuad4( numele++, nids ) );
    }
  }

  numnode += (rows+1)*(cols+1);
}

void create_quad4_cylinder_mesh( GEO::CUT::MeshIntersection & intersection, double x, double y, int rows, int cols )
{
  double r = 1.;

  int rownodes = rows+1;
  int colnodes = cols+1;

  for (int i=0; i<rownodes; ++i)
  {
    for (int j=0; j<colnodes; ++j)
    {
      int id = i+j*rownodes;
      double coord[3];

      double alpha = static_cast<double>( i ) / rows;

      coord[0] = x + r*cos( 2*M_PI*alpha );
      coord[1] = y + r*sin( 2*M_PI*alpha );
      coord[2] = static_cast<double>( j ) / cols;

      intersection.CutMesh().GetNode( numnode+id, coord );
    }
  }

  for (int i=0; i<rows+1; ++i)
  {
    for (int j=0; j<cols; ++j)
    {
      std::vector<int> nids;
      nids.reserve( 4 );
      nids.push_back( numnode+( ( i   )%( rows ) )+    j*rownodes );
      nids.push_back( numnode+( ( i+1 )%( rows ) )+    j*rownodes );
      nids.push_back( numnode+( ( i+1 )%( rows ) )+(j+1)*rownodes );
      nids.push_back( numnode+( ( i   )%( rows ) )+(j+1)*rownodes );

      intersection.AddCutSide( numele++, nids, DRT::Element::quad4 );
    }
  }

  numnode += rownodes*colnodes;
}
