
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_element.H"
#include "../drt_cut/cut_levelsetintersection.H"
#include "cut_test_utils.H"

void test_ls_hex8_simple()
{
  GEO::CUT::LevelSetIntersection lsi;

  // simple hex8 element
  std::vector<int> nids( 8 );
  std::vector<double> lsvs( 8 );
  Epetra_SerialDenseMatrix xyze( 3, 8 );

  for ( int i=0; i<8; ++i )
  {
    nids[i] = i;
  }

  std::fill( &lsvs[0], &lsvs[4], -1. );
  std::fill( &lsvs[4], &lsvs[8],  1. );

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

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_simple2()
{
  GEO::CUT::LevelSetIntersection lsi;

  // simple hex8 element
  std::vector<int> nids( 8 );
  std::vector<double> lsvs( 8, -1 );
  Epetra_SerialDenseMatrix xyze( 3, 8 );

  for ( int i=0; i<8; ++i )
  {
    nids[i] = i;
  }

  lsvs[3] = 1;
  lsvs[7] = 1;

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

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_simple3()
{
  GEO::CUT::LevelSetIntersection lsi;

  // simple hex8 element
  std::vector<int> nids( 8 );
  std::vector<double> lsvs( 8, -1 );
  Epetra_SerialDenseMatrix xyze( 3, 8 );

  for ( int i=0; i<8; ++i )
  {
    nids[i] = i;
  }

  lsvs[3] = 1;
  lsvs[4] = 1;

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

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_simple4()
{
  GEO::CUT::LevelSetIntersection lsi;

  // simple hex8 element
  std::vector<int> nids( 8 );
  std::vector<double> lsvs( 8, -1 );
  Epetra_SerialDenseMatrix xyze( 3, 8 );

  for ( int i=0; i<8; ++i )
  {
    nids[i] = i;
  }

  lsvs[1] = 1;
  lsvs[3] = 1;
  lsvs[4] = 1;
  lsvs[6] = 1;

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

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}
