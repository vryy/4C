
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_element.H"
#include "cut_test_utils.H"

GEO::CUT::Element* create_tet4( GEO::CUT::Mesh & mesh )
{
  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 2;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 2;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 1;

  xyze( 0, 2 ) = 2;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 0.5;

  return create_tet4( mesh, xyze );
}

GEO::CUT::Side* create_quad4( GEO::CUT::Mesh & mesh, double x, double dx, double dz )
{
  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) =  x - dx;
  xyze( 1, 0 ) = -0.5;
  xyze( 2, 0 ) = -0.5 - dz;

  xyze( 0, 1 ) =  x + dx;
  xyze( 1, 1 ) = -0.5;
  xyze( 2, 1 ) =  1.5 + dz;

  xyze( 0, 2 ) =  x + dx;
  xyze( 1, 2 ) =  1.5;
  xyze( 2, 2 ) =  1.5 + dz;

  xyze( 0, 3 ) =  x - dx;
  xyze( 1, 3 ) =  1.5;
  xyze( 2, 3 ) = -0.5 - dz;

  return create_quad4( mesh, xyze );
}

void test_hex8_simple()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * e = create_hex8( mesh );
  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 0.1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_tet4_simple()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1;

  GEO::CUT::Element * e = create_tet4( mesh, xyze );
  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 0.1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_pyramid5_simple()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 5 );

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

  xyze( 0, 4 ) = 0.5;
  xyze( 1, 4 ) = 1;
  xyze( 2, 4 ) = 1;

  GEO::CUT::Element * e = create_pyramid5( mesh, xyze );
  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 0.1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_wedge6_simple()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 6 );

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
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1;

  xyze( 0, 4 ) = 1;
  xyze( 1, 4 ) = 0;
  xyze( 2, 4 ) = 1;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = 1;
  xyze( 2, 5 ) = 1;

  GEO::CUT::Element * e = create_wedge6( mesh, xyze );
  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 0.1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}


void test_hex8_fullside()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * e = create_hex8( mesh );
  GEO::CUT::Side * s = create_quad4( mesh, 1, 0, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_diagonal()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * e = create_hex8( mesh );
  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_tet4()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );
  GEO::CUT::Element * tet4 = create_tet4( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  // add second cut to be able to find nodal positions

  xyze( 0, 0 ) = -0.1;
  xyze( 1, 0 ) =  0.1;
  xyze( 2, 0 ) = -0.1;

  xyze( 0, 1 ) =  1.1;
  xyze( 1, 1 ) =  0.1;
  xyze( 2, 1 ) = -0.1;

  xyze( 0, 2 ) =  1.1;
  xyze( 1, 2 ) = -0.1;
  xyze( 2, 2 ) =  0.1;

  xyze( 0, 3 ) = -0.1;
  xyze( 1, 3 ) = -0.1;
  xyze( 2, 3 ) =  0.1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=tet4->Sides().begin();
        i!=tet4->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();

  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8_1 = create_hex8( mesh );
  GEO::CUT::Element * hex8_2 = create_hex8( mesh, 0.5, 0.5, 0.5 );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_touch()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8_1 = create_hex8( mesh );
  GEO::CUT::Element * hex8_2 = create_hex8( mesh, 1, 0, 0 );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_touch2()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8_1 = create_hex8( mesh );
  GEO::CUT::Element * hex8_2 = create_hex8( mesh, 1, 0.5, 0.5 );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_schraeg()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = 1;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 1;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_tet4_touch()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 2;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 2;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 1;

  xyze( 0, 2 ) = 2;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 0.5;

  GEO::CUT::Element * tet4 = create_tet4( mesh, xyze );

  // add second cut to be able to find nodal positions

  xyze( 0, 0 ) = -0.1;
  xyze( 1, 0 ) =  0.1;
  xyze( 2, 0 ) = -0.1;

  xyze( 0, 1 ) =  1.1;
  xyze( 1, 1 ) =  0.1;
  xyze( 2, 1 ) = -0.1;

  xyze( 0, 2 ) =  1.1;
  xyze( 1, 2 ) = -0.1;
  xyze( 2, 2 ) =  0.1;

  xyze( 0, 3 ) = -0.1;
  xyze( 1, 3 ) = -0.1;
  xyze( 2, 3 ) =  0.1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=tet4->Sides().begin();
        i!=tet4->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_tet4_touch2()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.5;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = -0.5;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = -0.5;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1.5;

  GEO::CUT::Element * tet4 = create_tet4( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=tet4->Sides().begin();
        i!=tet4->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_mesh()
{
  GEO::CUT::Mesh mesh;

  create_hex8_mesh( mesh, 10, 10, 10 );

  GEO::CUT::Side * s = create_quad4( mesh, 0.5, 0.5, 0 );

  mesh.Status();

  std::set<GEO::CUT::Element*> done;
  std::set<GEO::CUT::Element*> elements_done;
  mesh.Cut( *dynamic_cast<GEO::CUT::LinearSide*>( s ), done, elements_done );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_double()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * e = create_hex8( mesh );
  GEO::CUT::Side * s1 = create_quad4( mesh, 0.4, 0.1, 0 );
  GEO::CUT::Side * s2 = create_quad4( mesh, 0.6, 0.1, 0 );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s1 ) );
  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s2 ) );

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_multiple()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * e = create_hex8( mesh );

  for ( int i=0; i<11; ++i )
  {
    double x = 0.1*i;
    GEO::CUT::Side * s = create_quad4( mesh, x, 0.1, 0 );

    mesh.Status();

    e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );

    mesh.Status();
  }

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_bad1()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 0.7291666666666666;
  xyze( 1, 0 ) = 0.5208333333332368;
  xyze( 2, 0 ) = 0.02500000000896939;

  xyze( 0, 1 ) = 0.7291666666666667;
  xyze( 1, 1 ) = 0.5208333333333334;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 0.75;
  xyze( 1, 2 ) = 0.5208333333333334;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0.7499999999999999;
  xyze( 1, 3 ) = 0.5208333333332476;
  xyze( 2, 3 ) = 0.02500000000797485;

  xyze( 0, 4 ) = 0.7291666666666667;
  xyze( 1, 4 ) = 0.5;
  xyze( 2, 4 ) = 0.025;

  xyze( 0, 5 ) = 0.7291666666666667;
  xyze( 1, 5 ) = 0.5;
  xyze( 2, 5 ) = 0;

  xyze( 0, 6 ) = 0.75;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = 0;

  xyze( 0, 7 ) = 0.75;
  xyze( 1, 7 ) = 0.5;
  xyze( 2, 7 ) = 0.025;

  GEO::CUT::Element * e = create_hex8( mesh, xyze );

  xyze( 0, 0 ) = 0.75;
  xyze( 1, 0 ) = 0.5010108360343256;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 0.7435592801990288;
  xyze( 1, 1 ) = 0.5208333333333334;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 0.7435592801990578;
  xyze( 1, 2 ) = 0.5208333333332442;
  xyze( 2, 2 ) = 0.02500000000828232;

  xyze( 0, 3 ) = 0.75;
  xyze( 1, 3 ) = 0.5010108360343257;
  xyze( 2, 3 ) = 0.02500000000038694;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();
  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_bad2()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 1.05556;
  xyze( 1, 0 ) = 0.444444;
  xyze( 2, 0 ) = -4.82103e-20;

  xyze( 0, 1 ) = 1.05556;
  xyze( 1, 1 ) = 0.444444;
  xyze( 2, 1 ) = -0.05;

  xyze( 0, 2 ) = 1.05556;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = -0.05;

  xyze( 0, 3 ) = 1.05556;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 0;

  xyze( 0, 4 ) = 1.11111;
  xyze( 1, 4 ) = 0.444444;
  xyze( 2, 4 ) = 1.41172e-22;

  xyze( 0, 5 ) = 1.11111;
  xyze( 1, 5 ) = 0.444444;
  xyze( 2, 5 ) = -0.05;

  xyze( 0, 6 ) = 1.11111;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = -0.05;

  xyze( 0, 7 ) = 1.11111;
  xyze( 1, 7 ) = 0.5;
  xyze( 2, 7 ) = 0;

  GEO::CUT::Element * e = create_hex8( mesh, xyze );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.0505;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1.05714;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = -9.3343e-19;

  xyze( 0, 3 ) = 1.05714;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = -0.0505;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();
  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_bad3()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 1.05556;
  xyze( 1, 0 ) = 0.444444;
  xyze( 2, 0 ) = 0.05;

  xyze( 0, 1 ) = 1.05556;
  xyze( 1, 1 ) = 0.444444;
  xyze( 2, 1 ) = -4.82103e-20;

  xyze( 0, 2 ) = 1.05556;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 1.05556;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 0.05;

  xyze( 0, 4 ) = 1.11111;
  xyze( 1, 4 ) = 0.444444;
  xyze( 2, 4 ) = 0.05;

  xyze( 0, 5 ) = 1.11111;
  xyze( 1, 5 ) = 0.444444;
  xyze( 2, 5 ) = 1.41172e-22;

  xyze( 0, 6 ) = 1.11111;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = 0;

  xyze( 0, 7 ) = 1.11111;
  xyze( 1, 7 ) = 0.5;
  xyze( 2, 7 ) = 0.05;

  GEO::CUT::Element * e = create_hex8( mesh, xyze );

  xyze( 0, 0 ) = 1.05714;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -9.3343e-19;

  xyze( 0, 1 ) = 1.05714;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = 0.0505;

  xyze( 0, 2 ) = 1.11429;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 0.0505;

  xyze( 0, 3 ) = 1.11429;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1.60089e-18;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();
  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

/*
 * (0.944444,0,0.05), (0.944444,0,6.55604e-19),
 * (0.944444,0.0555556,1.29045e-18), (0.944444,0.0555556,0.05), (1,0,0.05),
 * (1,0,0), (1,0.0555556,3.40507e-19), (1,0.0555556,0.05)
 *
 * {(1,0,0.0505), (1,0.0555556,0.0505), (1,0.0555556,3.40507e-19), (1,0,0)}
 */

void test_hex8_bad4()
{
  GEO::CUT::Mesh mesh;

  double hex8_xyz[24] = {
    0.944444,0,0.05,
    0.944444,0,6.55604e-19,
    0.944444,0.0555556,1.29045e-18,
    0.944444,0.0555556,0.05,
    1,0,0.05,
    1,0,0,
    1,0.0555556,3.40507e-19,
    1,0.0555556,0.05
  };

  double quad4_xyz[12] = {
    1,0,0.0505,
    1,0.0555556,0.0505,
    1,0.0555556,3.40507e-19,
    1,0,0
  };

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  std::copy( hex8_xyz, hex8_xyz+24, &xyze( 0, 0 ) );

  GEO::CUT::Element * e = create_hex8( mesh, xyze );

  std::copy( quad4_xyz, quad4_xyz+12, &xyze( 0, 0 ) );

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  e->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();
  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}

void test_hex8_wedge6()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 6 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 2;
  xyze( 2, 0 ) = -0.5;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = -0.5;

  xyze( 0, 2 ) = 3;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = -0.5;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 2;
  xyze( 2, 3 ) = 1.5;

  xyze( 0, 4 ) = 0.5;
  xyze( 1, 4 ) = 0.5;
  xyze( 2, 4 ) = 1.5;

  xyze( 0, 5 ) = 3;
  xyze( 1, 5 ) = 0.5;
  xyze( 2, 5 ) = 1.5;

  GEO::CUT::Element * wedge6 = create_wedge6( mesh, xyze );

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=wedge6->Sides().begin();
        i!=wedge6->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 1.5;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1.5;
  xyze( 2, 2 ) = 1.5;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch2()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 1.5;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1.5;
  xyze( 2, 2 ) = 1.5;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch3()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.5;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 1.5;
  xyze( 2, 1 ) = -0.5;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1.5;
  xyze( 2, 2 ) = 1.5;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_cut()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = 1;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 0.5;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 1;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_gedreht()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

#if 0
  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.2;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = 1.2;
  xyze( 2, 1 ) = 0.5;

  xyze( 0, 2 ) = 0.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 1.2;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = -0.2;
  xyze( 2, 3 ) = 0.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );
#endif

  std::vector<GEO::CUT::Side*> sides;
  create_quad4_mesh( mesh, 2, 2, sides );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    GEO::CUT::Side* quad4 = *i;
    hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );
  }

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8_durchstoss()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8_1 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = -0.5;
  xyze( 1, 0 ) =  0.2;
  xyze( 2, 0 ) =  0.2;

  xyze( 0, 1 ) = -0.5;
  xyze( 1, 1 ) =  0.8;
  xyze( 2, 1 ) =  0.2;

  xyze( 0, 2 ) = -0.5;
  xyze( 1, 2 ) =  0.8;
  xyze( 2, 2 ) =  0.8;

  xyze( 0, 3 ) = -0.5;
  xyze( 1, 3 ) =  0.2;
  xyze( 2, 3 ) =  0.8;

  xyze( 0, 4 ) =  1.5;
  xyze( 1, 4 ) =  0.2;
  xyze( 2, 4 ) =  0.2;

  xyze( 0, 5 ) =  1.5;
  xyze( 1, 5 ) =  0.8;
  xyze( 2, 5 ) =  0.2;

  xyze( 0, 6 ) =  1.5;
  xyze( 1, 6 ) =  0.8;
  xyze( 2, 6 ) =  0.8;

  xyze( 0, 7 ) =  1.5;
  xyze( 1, 7 ) =  0.2;
  xyze( 2, 7 ) =  0.8;

  GEO::CUT::Element * hex8_2 = create_hex8( mesh, xyze );

  // add second cut to be able to find nodal positions

  xyze( 0, 0 ) = -0.1;
  xyze( 1, 0 ) =  0.1;
  xyze( 2, 0 ) = -0.1;

  xyze( 0, 1 ) =  1.1;
  xyze( 1, 1 ) =  0.1;
  xyze( 2, 1 ) = -0.1;

  xyze( 0, 2 ) =  1.1;
  xyze( 1, 2 ) = -0.1;
  xyze( 2, 2 ) =  0.1;

  xyze( 0, 3 ) = -0.1;
  xyze( 1, 3 ) = -0.1;
  xyze( 2, 3 ) =  0.1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8_onside()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8_1 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0.2;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0.1;
  xyze( 2, 1 ) = 0.2;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 0.2;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.9;
  xyze( 2, 3 ) = 0.2;

  xyze( 0, 4 ) = 0.5;
  xyze( 1, 4 ) = 0.5;
  xyze( 2, 4 ) = 0.8;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = 0.1;
  xyze( 2, 5 ) = 0.8;

  xyze( 0, 6 ) = 1.5;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = 0.8;

  xyze( 0, 7 ) = 1;
  xyze( 1, 7 ) = 0.9;
  xyze( 2, 7 ) = 0.8;

  GEO::CUT::Element * hex8_2 = create_hex8( mesh, xyze );

  // add second cut to be able to find nodal positions

  xyze( 0, 0 ) = -0.1;
  xyze( 1, 0 ) =  0.1;
  xyze( 2, 0 ) = -0.1;

  xyze( 0, 1 ) =  1.1;
  xyze( 1, 1 ) =  0.1;
  xyze( 2, 1 ) = -0.1;

  xyze( 0, 2 ) =  1.1;
  xyze( 1, 2 ) = -0.1;
  xyze( 2, 2 ) =  0.1;

  xyze( 0, 3 ) = -0.1;
  xyze( 1, 3 ) = -0.1;
  xyze( 2, 3 ) =  0.1;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8_internal()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = -1;
  xyze( 1, 0 ) = -1;
  xyze( 2, 0 ) = -1;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = -1;
  xyze( 2, 1 ) = -1;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = -1;

  xyze( 0, 3 ) = -1;
  xyze( 1, 3 ) = 1;
  xyze( 2, 3 ) = -1;

  xyze( 0, 4 ) = -1;
  xyze( 1, 4 ) = -1;
  xyze( 2, 4 ) = 1;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = -1;
  xyze( 2, 5 ) = 1;

  xyze( 0, 6 ) = 1;
  xyze( 1, 6 ) = 1;
  xyze( 2, 6 ) = 1;

  xyze( 0, 7 ) = -1;
  xyze( 1, 7 ) = 1;
  xyze( 2, 7 ) = 1;

  GEO::CUT::Element * hex8_1 = create_hex8( mesh, xyze );

  xyze( 0, 0 ) = -1.5;
  xyze( 1, 0 ) = -0.5;
  xyze( 2, 0 ) = 0.707107;

  xyze( 0, 1 ) = -0.5;
  xyze( 1, 1 ) = -1.5;
  xyze( 2, 1 ) = -0.707107;

  xyze( 0, 2 ) = -0.207107;
  xyze( 1, 2 ) = 0.207107;
  xyze( 2, 2 ) = -1.70711;

  xyze( 0, 3 ) = -1.20711;
  xyze( 1, 3 ) = 1.20711;
  xyze( 2, 3 ) = -0.292893;

  xyze( 0, 4 ) = 0.207107;
  xyze( 1, 4 ) = -0.207107;
  xyze( 2, 4 ) = 1.70711;

  xyze( 0, 5 ) = 1.20711;
  xyze( 1, 5 ) = -1.20711;
  xyze( 2, 5 ) = 0.292893;

  xyze( 0, 6 ) = 1.5;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = -0.707107;

  xyze( 0, 7 ) = 0.5;
  xyze( 1, 7 ) = 1.5;
  xyze( 2, 7 ) = 0.707107;

  GEO::CUT::Element * hex8_2 = create_hex8( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );
  hex8_1->FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8_sideintersection()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = -1;
  xyze( 1, 0 ) = -1;
  xyze( 2, 0 ) = -1;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = -1;
  xyze( 2, 1 ) = -1;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = -1;

  xyze( 0, 3 ) = -1;
  xyze( 1, 3 ) = 1;
  xyze( 2, 3 ) = -1;

  xyze( 0, 4 ) = -1;
  xyze( 1, 4 ) = -1;
  xyze( 2, 4 ) = 1;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = -1;
  xyze( 2, 5 ) = 1;

  xyze( 0, 6 ) = 1;
  xyze( 1, 6 ) = 1;
  xyze( 2, 6 ) = 1;

  xyze( 0, 7 ) = -1;
  xyze( 1, 7 ) = 1;
  xyze( 2, 7 ) = 1;

  GEO::CUT::Element * hex8_1 = create_hex8( mesh, xyze );

  xyze( 0, 0 ) =  0.5;
  xyze( 1, 0 ) = -0.5;
  xyze( 2, 0 ) = -0.5;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = -0.5;
  xyze( 2, 1 ) = -0.5;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = -0.5;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = -0.5;

  xyze( 0, 4 ) = 0.5;
  xyze( 1, 4 ) = -0.5;
  xyze( 2, 4 ) = 0.5;

  xyze( 0, 5 ) = 1.5;
  xyze( 1, 5 ) = -0.5;
  xyze( 2, 5 ) = 0.5;

  xyze( 0, 6 ) = 1.5;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = 0.5;

  xyze( 0, 7 ) = 0.5;
  xyze( 1, 7 ) = 0.5;
  xyze( 2, 7 ) = 0.5;


  GEO::CUT::Element * hex8_2 = create_hex8( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );

  mesh.FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_hex8_inside()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = -1;
  xyze( 1, 0 ) = -1;
  xyze( 2, 0 ) = -1;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = -1;
  xyze( 2, 1 ) = -1;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = -1;

  xyze( 0, 3 ) = -1;
  xyze( 1, 3 ) = 1;
  xyze( 2, 3 ) = -1;

  xyze( 0, 4 ) = -1;
  xyze( 1, 4 ) = -1;
  xyze( 2, 4 ) = 1;

  xyze( 0, 5 ) = 1;
  xyze( 1, 5 ) = -1;
  xyze( 2, 5 ) = 1;

  xyze( 0, 6 ) = 1;
  xyze( 1, 6 ) = 1;
  xyze( 2, 6 ) = 1;

  xyze( 0, 7 ) = -1;
  xyze( 1, 7 ) = 1;
  xyze( 2, 7 ) = 1;

  GEO::CUT::Element * hex8_1 = create_hex8( mesh, xyze );

  xyze( 0, 0 ) = -0.5;
  xyze( 1, 0 ) = -0.5;
  xyze( 2, 0 ) = -0.5;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = -0.5;
  xyze( 2, 1 ) = -0.5;

  xyze( 0, 2 ) = 0.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = -0.5;

  xyze( 0, 3 ) = -0.5;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = -0.5;

  xyze( 0, 4 ) = -0.5;
  xyze( 1, 4 ) = -0.5;
  xyze( 2, 4 ) = 0.5;

  xyze( 0, 5 ) = 0.5;
  xyze( 1, 5 ) = -0.5;
  xyze( 2, 5 ) = 0.5;

  xyze( 0, 6 ) = 0.5;
  xyze( 1, 6 ) = 0.5;
  xyze( 2, 6 ) = 0.5;

  xyze( 0, 7 ) = -0.5;
  xyze( 1, 7 ) = 0.5;
  xyze( 2, 7 ) = 0.5;


  GEO::CUT::Element * hex8_2 = create_hex8( mesh, xyze );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::const_iterator i=hex8_2->Sides().begin();
        i!=hex8_2->Sides().end();
        ++i )
  {
    GEO::CUT::Side * s = *i;
    s->SetId( 1 );
    hex8_1->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( s ) );
  }

  mesh.Status();

  hex8_1->MakeFacets( mesh );

  mesh.FindNodePositions();

  hex8_1->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_schnitt()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.2;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = -0.2;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 1.2;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 0.5;
  xyze( 2, 3 ) = 1.2;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch4()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.2;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 0;
  xyze( 2, 2 ) = 1.2;

  xyze( 0, 3 ) = 0.2;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1.2;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch5()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.2;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 0;
  xyze( 2, 2 ) = 1.2;

  xyze( 0, 3 ) = 1.2;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 1.2;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch6()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 1;
  xyze( 2, 1 ) = 0.5;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 1;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0;
  xyze( 2, 3 ) = 0.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  // add second cut to be able to find nodal positions

  xyze( 0, 0 ) =  0.1;
  xyze( 1, 0 ) = -0.1;
  xyze( 2, 0 ) = -0.1;

  xyze( 0, 1 ) =  0.1;
  xyze( 1, 1 ) =  1.1;
  xyze( 2, 1 ) = -0.1;

  xyze( 0, 2 ) = -0.1;
  xyze( 1, 2 ) =  1.1;
  xyze( 2, 2 ) =  0.1;

  xyze( 0, 3 ) = -0.1;
  xyze( 1, 3 ) = -0.1;
  xyze( 2, 3 ) =  0.1;

  GEO::CUT::Side* quad4_2 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );
  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4_2 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

void test_hex8_quad4_touch7()
{
  GEO::CUT::Mesh mesh;
  GEO::CUT::Element * hex8 = create_hex8( mesh );

  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 1;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0.8;
  xyze( 2, 1 ) = 0.5;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 1;

  xyze( 0, 3 ) = 1;
  xyze( 1, 3 ) = 0.2;
  xyze( 2, 3 ) = 0.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );

  mesh.Status();

  hex8->Cut( mesh, *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ) );

  mesh.Status();

  hex8->MakeFacets( mesh );
  hex8->FindNodePositions();
  hex8->GenerateTetgen( mesh, NULL );
}

#if 0
void test_quad4_quad4_simple()
{
  GEO::CUT::Mesh mesh;

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

  GEO::CUT::Side* quad4_1 = create_quad4( mesh, xyze );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 1.5;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 1.5;
  xyze( 2, 3 ) = 0;

  GEO::CUT::Side* quad4_2 = create_quad4( mesh, xyze );

  mesh.Status();

  quad4_1->Cut( mesh, *quad4_2, NULL );

  mesh.Status();
}
#endif

void test_hex8_quad4_mesh()
{
  GEO::CUT::Mesh mesh;

  create_hex8_mesh( mesh, 2, 2, 2 );

  std::vector<GEO::CUT::Side*> sides;
  create_quad4_mesh( mesh, 3, 3, sides );

  mesh.Status();

  for ( std::vector<GEO::CUT::Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    GEO::CUT::Side* quad4 = *i;
    std::set<GEO::CUT::Element*> done;
    std::set<GEO::CUT::Element*> elements_done;
    mesh.Cut( *dynamic_cast<GEO::CUT::LinearSide*>( quad4 ), done, elements_done );
  }

  mesh.Status();

  mesh.MakeFacets();
  mesh.FindNodePositions();
  mesh.GenerateTetgen();
}
