
#include "cut_test_utils.H"
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_element.H"

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
