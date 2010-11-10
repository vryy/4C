
#include <iostream>

#include "cut_test_utils.H"

#include "../drt_cut/cut_meshintersection.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_hex27_quad9_simple()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad9_xyze( 3, 9 );
  Epetra_SerialDenseMatrix hex27_xyze( 3, 27 );

  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<2; ++j )
    {
      quad9_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][j];
    }
    quad9_xyze( 2, i ) = 0;
  }

  intersection.AddCutSide( 1, nids, quad9_xyze, DRT::Element::quad9 );

  nids.clear();
  nids.reserve( 27 );
  for ( int i=0; i<27; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<3; ++j )
    {
      hex27_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[i][j];
    }
  }

  // shrink hex27 element

  for ( int i=0; i<27; ++i )
  {
    for ( int j=0; j<3; ++j )
    {
      hex27_xyze( j, i ) *= 0.5;
    }
  }

  intersection.AddElement( 1, nids, hex27_xyze, DRT::Element::hex27 );

  intersection.Cut( NULL );
}

void test_hex20_quad9_simple()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad9_xyze( 3, 9 );
  Epetra_SerialDenseMatrix hex20_xyze( 3, 20 );

  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<2; ++j )
    {
      quad9_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][j];
    }
    quad9_xyze( 2, i ) = 0;
  }

  intersection.AddCutSide( 1, nids, quad9_xyze, DRT::Element::quad9 );

  nids.clear();
  nids.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<3; ++j )
    {
      hex20_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[i][j];
    }
  }

  // shrink hex20 element

  for ( int i=0; i<20; ++i )
  {
    for ( int j=0; j<3; ++j )
    {
      hex20_xyze( j, i ) *= 0.5;
    }
  }

  intersection.AddElement( 1, nids, hex20_xyze, DRT::Element::hex20 );

  intersection.Cut( NULL );
}

void test_hex20_quad9_moved()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad9_xyze( 3, 9 );
  Epetra_SerialDenseMatrix hex20_xyze( 3, 20 );

  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<2; ++j )
    {
      quad9_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][j];
    }
    quad9_xyze( 2, i ) = 0;
  }

  // move quad9 element

  for ( int i=0; i<9; ++i )
  {
    quad9_xyze( 2, i ) = 0.1 + 0.5*quad9_xyze( 0, i );
    quad9_xyze( 0, i ) += 0.1;
    quad9_xyze( 1, i ) += 0.1;
  }

  intersection.AddCutSide( 1, nids, quad9_xyze, DRT::Element::quad9 );

  nids.clear();
  nids.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<3; ++j )
    {
      hex20_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[i][j];
    }
  }

  // shrink hex20 element

  for ( int i=0; i<20; ++i )
  {
    for ( int j=0; j<3; ++j )
    {
      hex20_xyze( j, i ) *= 0.5;
    }
  }

  intersection.AddElement( 1, nids, hex20_xyze, DRT::Element::hex20 );

  intersection.Cut( NULL );
}

void test_tet10_quad9_simple()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad9_xyze( 3, 9 );
  Epetra_SerialDenseMatrix tet10_xyze( 3, 10 );

  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<2; ++j )
    {
      quad9_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][j];
    }
    quad9_xyze( 2, i ) = 0.2;
  }

  intersection.AddCutSide( 1, nids, quad9_xyze, DRT::Element::quad9 );

  nids.clear();
  nids.reserve( 10 );
  for ( int i=0; i<10; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<3; ++j )
    {
      tet10_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_tet10_nodes_reference[i][j];
    }
  }

  // shrink tet10 element

  for ( int i=0; i<10; ++i )
  {
    for ( int j=0; j<3; ++j )
    {
      tet10_xyze( j, i ) *= 0.5;
    }
  }

  intersection.AddElement( 1, nids, tet10_xyze, DRT::Element::tet10 );

  intersection.Cut( NULL );
}

void test_tet10_quad9_moved()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad9_xyze( 3, 9 );
  Epetra_SerialDenseMatrix tet10_xyze( 3, 10 );

  nids.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<2; ++j )
    {
      quad9_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[i][j];
    }
    quad9_xyze( 2, i ) = 0;
  }

  // move quad9 element

  for ( int i=0; i<9; ++i )
  {
    quad9_xyze( 2, i ) = 0.1 + 0.5*quad9_xyze( 0, i );
    quad9_xyze( 0, i ) += 0.1;
    quad9_xyze( 1, i ) += 0.1;
  }

  intersection.AddCutSide( 1, nids, quad9_xyze, DRT::Element::quad9 );

  nids.clear();
  nids.reserve( 10 );
  for ( int i=0; i<10; ++i )
  {
    nids.push_back( i );
    for ( int j=0; j<3; ++j )
    {
      tet10_xyze( j, i ) = DRT::UTILS::eleNodeNumbering_tet10_nodes_reference[i][j];
    }
  }

  // shrink tet10 element

  for ( int i=0; i<10; ++i )
  {
    for ( int j=0; j<3; ++j )
    {
      tet10_xyze( j, i ) *= 0.5;
    }
  }

  intersection.AddElement( 1, nids, tet10_xyze, DRT::Element::tet10 );

  intersection.Cut( NULL );
}
