
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_levelsetintersection.H"
#include "cut_test_utils.H"

void test_ls_hex8_florian1()
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

  lsvs[0] = 0.331317;
  lsvs[1] = 0.331355;
  lsvs[2] = -0.0718546;
  lsvs[3] = -0.0719582;
  lsvs[4] = -1.06152;
  lsvs[5] = -1.06155;
  lsvs[6] = 0.178732;
  lsvs[7] = 0.17888;

  xyze( 0, 0 ) = 92.3077;
  xyze( 1, 0 ) = 122.481;
  xyze( 2, 0 ) = 0.75   ;

  xyze( 0, 1 ) = 92.3077;
  xyze( 1, 1 ) = 122.481;
  xyze( 2, 1 ) = -0.75  ;

  xyze( 0, 2 ) = 92.3077;
  xyze( 1, 2 ) = 124.031;
  xyze( 2, 2 ) = -0.75  ;

  xyze( 0, 3 ) = 92.3077;
  xyze( 1, 3 ) = 124.031;
  xyze( 2, 3 ) = 0.75   ;

  xyze( 0, 4 ) = 93.8462;
  xyze( 1, 4 ) = 122.481;
  xyze( 2, 4 ) = 0.75   ;

  xyze( 0, 5 ) = 93.8462;
  xyze( 1, 5 ) = 122.481;
  xyze( 2, 5 ) = -0.75  ;

  xyze( 0, 6 ) = 93.8462;
  xyze( 1, 6 ) = 124.031;
  xyze( 2, 6 ) = -0.75  ;

  xyze( 0, 7 ) = 93.8462;
  xyze( 1, 7 ) = 124.031;
  xyze( 2, 7 ) = 0.75   ;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_florian2()
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

  lsvs[0] = -1.06142;
  lsvs[1] = -1.06145;
  lsvs[2] = 0.178845;
  lsvs[3] = 0.178992;
  lsvs[4] = 0.331525;
  lsvs[5] = 0.331563;
  lsvs[6] = -0.071624;
  lsvs[7] = -0.0717273;

  xyze( 0, 0 ) = 6.15385;
  xyze( 1, 0 ) = 122.481;
  xyze( 2, 0 ) = 0.75   ;

  xyze( 0, 1 ) = 6.15385;
  xyze( 1, 1 ) = 122.481;
  xyze( 2, 1 ) = -0.75  ;

  xyze( 0, 2 ) = 6.15385;
  xyze( 1, 2 ) = 124.031;
  xyze( 2, 2 ) = -0.75  ;

  xyze( 0, 3 ) = 6.15385;
  xyze( 1, 3 ) = 124.031;
  xyze( 2, 3 ) = 0.75   ;

  xyze( 0, 4 ) = 7.69231;
  xyze( 1, 4 ) = 122.481;
  xyze( 2, 4 ) = 0.75   ;

  xyze( 0, 5 ) = 7.69231;
  xyze( 1, 5 ) = 122.481;
  xyze( 2, 5 ) = -0.75  ;

  xyze( 0, 6 ) = 7.69231;
  xyze( 1, 6 ) = 124.031;
  xyze( 2, 6 ) = -0.75  ;

  xyze( 0, 7 ) = 7.69231;
  xyze( 1, 7 ) = 124.031;
  xyze( 2, 7 ) = 0.75   ;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_florian3()
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

  lsvs[0] = -0.0524131;
  lsvs[1] = -0.0566379;
  lsvs[2] = 0.0747609;
  lsvs[3] = 0.0696167;
  lsvs[4] = 0.0151905;
  lsvs[5] = 0.00655406;
  lsvs[6] = -0.0626401;
  lsvs[7] = -0.11719;

  xyze( 0, 0 ) = 49.2308;
  xyze( 1, 0 ) = 100.775;
  xyze( 2, 0 ) = 0.75   ;

  xyze( 0, 1 ) = 49.2308;
  xyze( 1, 1 ) = 100.775;
  xyze( 2, 1 ) = -0.75  ;

  xyze( 0, 2 ) = 49.2308;
  xyze( 1, 2 ) = 102.326;
  xyze( 2, 2 ) = -0.75  ;

  xyze( 0, 3 ) = 49.2308;
  xyze( 1, 3 ) = 102.326;
  xyze( 2, 3 ) = 0.75   ;

  xyze( 0, 4 ) = 50.7692;
  xyze( 1, 4 ) = 100.775;
  xyze( 2, 4 ) = 0.75   ;

  xyze( 0, 5 ) = 50.7692;
  xyze( 1, 5 ) = 100.775;
  xyze( 2, 5 ) = -0.75  ;

  xyze( 0, 6 ) = 50.7692;
  xyze( 1, 6 ) = 102.326;
  xyze( 2, 6 ) = -0.75  ;

  xyze( 0, 7 ) = 50.7692;
  xyze( 1, 7 ) = 102.326;
  xyze( 2, 7 ) = 0.75   ;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_florian4()
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

  lsvs[0] =   -0.0102494;
  lsvs[1] =   0.0126085;
  lsvs[2] =   -0.0312269;
  lsvs[3] =   0.0440133;
  lsvs[4] =   0.0628933;
  lsvs[5] =   0.0952024;
  lsvs[6] =   -0.366343;
  lsvs[7] =   -0.322604;

  xyze( 0, 0 ) = 47.0588;
  xyze( 1, 0 ) = 109.091;
  xyze( 2, 0 ) = 2.94118;

  xyze( 0, 1 ) = 47.0588 ;
  xyze( 1, 1 ) = 109.091 ;
  xyze( 2, 1 ) = -2.94118;

  xyze( 0, 2 ) = 47.0588 ;
  xyze( 1, 2 ) = 115.152 ;
  xyze( 2, 2 ) = -2.94118;

  xyze( 0, 3 ) = 47.0588;
  xyze( 1, 3 ) = 115.152;
  xyze( 2, 3 ) = 2.94118;

  xyze( 0, 4 ) = 52.9412;
  xyze( 1, 4 ) = 109.091;
  xyze( 2, 4 ) = 2.94118;

  xyze( 0, 5 ) = 52.9412 ;
  xyze( 1, 5 ) = 109.091 ;
  xyze( 2, 5 ) = -2.94118;

  xyze( 0, 6 ) = 52.9412 ;
  xyze( 1, 6 ) = 115.152 ;
  xyze( 2, 6 ) = -2.94118;

  xyze( 0, 7 ) = 52.9412;
  xyze( 1, 7 ) = 115.152;
  xyze( 2, 7 ) = 2.94118;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_florian5()
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

  lsvs[0] = 0.0986988;
  lsvs[1] = 0.00224424;
  lsvs[2] = -0.00388199;
  lsvs[3] = 0.17228;
  lsvs[4] = -0.116327;
  lsvs[5] = 0.0639577;
  lsvs[6] = 0.042944;
  lsvs[7] = -0.307789;

  xyze( 0, 0 ) = 49.2308 ;
  xyze( 1, 0 ) = 83.7209 ;
  xyze( 2, 0 ) = 0.769231;

  xyze( 0, 1 ) = 49.2308  ;
  xyze( 1, 1 ) = 83.7209  ;
  xyze( 2, 1 ) = -0.769231;

  xyze( 0, 2 ) = 49.2308  ;
  xyze( 1, 2 ) = 85.2713  ;
  xyze( 2, 2 ) = -0.769231;

  xyze( 0, 3 ) = 49.2308 ;
  xyze( 1, 3 ) = 85.2713 ;
  xyze( 2, 3 ) = 0.769231;

  xyze( 0, 4 ) = 50.7692 ;
  xyze( 1, 4 ) = 83.7209 ;
  xyze( 2, 4 ) = 0.769231;

  xyze( 0, 5 ) = 50.7692  ;
  xyze( 1, 5 ) = 83.7209  ;
  xyze( 2, 5 ) = -0.769231;

  xyze( 0, 6 ) = 50.7692  ;
  xyze( 1, 6 ) = 85.2713  ;
  xyze( 2, 6 ) = -0.769231;

  xyze( 0, 7 ) = 50.7692 ;
  xyze( 1, 7 ) = 85.2713 ;
  xyze( 2, 7 ) = 0.769231;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

void test_ls_hex8_florian6()
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

  lsvs[0] = -0.0780914;
  lsvs[1] = 0;
  lsvs[2] = -0.0798239;
  lsvs[3] = 0.214841;
  lsvs[4] = -0.390728;
  lsvs[5] = -0.381449;
  lsvs[6] = 0.29728;
  lsvs[7] = 0.333705;

  xyze( 0, 0 ) = 49.2308 ;
  xyze( 1, 0 ) = 83.7209 ;
  xyze( 2, 0 ) = 0.769231;

  xyze( 0, 1 ) = 49.2308  ;
  xyze( 1, 1 ) = 83.7209  ;
  xyze( 2, 1 ) = -0.769231;

  xyze( 0, 2 ) = 49.2308  ;
  xyze( 1, 2 ) = 85.2713  ;
  xyze( 2, 2 ) = -0.769231;

  xyze( 0, 3 ) = 49.2308 ;
  xyze( 1, 3 ) = 85.2713 ;
  xyze( 2, 3 ) = 0.769231;

  xyze( 0, 4 ) = 50.7692 ;
  xyze( 1, 4 ) = 83.7209 ;
  xyze( 2, 4 ) = 0.769231;

  xyze( 0, 5 ) = 50.7692  ;
  xyze( 1, 5 ) = 83.7209  ;
  xyze( 2, 5 ) = -0.769231;

  xyze( 0, 6 ) = 50.7692  ;
  xyze( 1, 6 ) = 85.2713  ;
  xyze( 2, 6 ) = -0.769231;

  xyze( 0, 7 ) = 50.7692 ;
  xyze( 1, 7 ) = 85.2713 ;
  xyze( 2, 7 ) = 0.769231;

  lsi.AddElement( 1, nids, xyze, &lsvs[0], DRT::Element::hex8 );
  lsi.Cut();
}

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

  // this is the impossible (undefined) case

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
