#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_position2d.H"
#include "cut_test_utils.H"
#include "../../src/drt_cut/cut_triangulateFacet.H"

void check4nodedInline( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check4nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check5nodedInline( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check5nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check5nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check5nodedAdjacentconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check6nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check6nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check7nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check7nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check8nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check8nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check8nodedAdjacentconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check9nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check9nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check10nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check10nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );

void check5nodedTwinConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check6nodedTwinConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check8nodedTriConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );
void check8nodedTriConcaveGenPlane( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );

void check13nodedConvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );

void check7nodedconti3concave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );

void check10nodedShiftEarClipToSplit( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s );

GEO::CUT::Side* Create_quad4( GEO::CUT::Mesh & mesh, double x, double dx, double dz, bool reverse=false )
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

  if ( reverse )
  {
    std::swap( xyze( 0, 1 ), xyze( 0, 3 ) );
    std::swap( xyze( 1, 1 ), xyze( 1, 3 ) );
    std::swap( xyze( 2, 1 ), xyze( 2, 3 ) );
  }

  return create_quad4( mesh, xyze );
}

void test_facet_split()
{
  std::cout<<"checking triangulation of facets...\n";
  GEO::CUT::Options options;
  GEO::CUT::Mesh mesh(options);
  GEO::CUT::Element * e = create_hex8( mesh );
  GEO::CUT::Side * s = Create_quad4( mesh, 0.5, 0.1, 0 );

  check4nodedInline( mesh, e, s );
  check4nodedconcave( mesh, e, s );
  check5nodedInline( mesh, e, s );
  check5nodedconvex( mesh, e, s );
  check5nodedconcave( mesh, e, s );
  check5nodedAdjacentconcave( mesh, e, s );
  check6nodedconvex( mesh, e, s );
  check6nodedconcave( mesh, e, s );
  check7nodedconvex( mesh, e, s );
  check7nodedconcave( mesh, e, s );
  check8nodedconvex( mesh, e, s );
  check8nodedconcave( mesh, e, s );
  check8nodedAdjacentconcave( mesh, e, s );
  check9nodedconvex( mesh, e, s );
  check9nodedconcave( mesh, e, s );
  check10nodedconvex( mesh, e, s );
  check10nodedconcave( mesh, e, s );

  check5nodedTwinConcave( mesh, e, s );
  check6nodedTwinConcave( mesh, e, s );
  check8nodedTriConcave( mesh, e, s );
  check8nodedTriConcaveGenPlane( mesh, e, s );

  check13nodedConvex( mesh, e, s );

  check7nodedconti3concave( mesh, e, s );

  check10nodedShiftEarClipToSplit( mesh, e, s );
}

void check4nodedInline( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check4nodedInline...\n";
  std::vector<GEO::CUT::Point*> ptlist(4);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.5;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p4 )
    dserror( "triangulation failed for check4nodedInline" );
}

void check4nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check4nodedConcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(4);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.5;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 2.0;x[1] = 0.0;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p2 || cell1[1]!=p3 || cell1[2]!=p4 )
    dserror( "triangulation failed for check4nodedInline" );
  if( cell2[0]!=p2 || cell2[1]!=p4 || cell2[2]!=p1 )
    dserror( "triangulation failed for check4nodedInline" );
}

void check5nodedInline( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check5nodedInline...\n";
  std::vector<GEO::CUT::Point*> ptlist(5);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.5;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.0;x[1] = 1.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p4 ||cell1[3]!=p5 )
    dserror( "triangulation failed for check5nodedInline" );
}

void check5nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check5nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(5);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.2;x[1] = 0.9;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 0.5;x[1] = 1.5;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3  )
    dserror( "triangulation failed for check5nodedconvex" );
  if( cell2[0]!=p3 || cell2[1]!=p4 || cell2[2]!=p5 || cell2[3]!=p1  )
    dserror( "triangulation failed for check5nodedconvex" );
}

void check5nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check5nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(5);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.2;x[1] = 0.9;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 0.5;x[1] = 0.7;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p4 || cell1[1]!=p5 || cell1[2]!=p1  )
    dserror( "triangulation failed for check5nodedconcave" );
  if( cell2[0]!=p1 || cell2[1]!=p2 || cell2[2]!=p3 || cell2[3]!=p4  )
    dserror( "triangulation failed for check5nodedconcave" );
}

void check5nodedAdjacentconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check5nodedAdjacentconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(5);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.4;x[1] = 0.2;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.6;x[1] = 0.2;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 2.0;x[1] = 0.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

/*  for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p5 || cell1[1]!=p1 || cell1[2]!=p2 )
    dserror( "triangulation failed for check5nodedAdjacentconcave" );
  if( cell2[0]!=p3 || cell2[1]!=p4 || cell2[2]!=p5 )
    dserror( "triangulation failed for check5nodedAdjacentconcave" );
  if( cell3[0]!=p3 || cell3[1]!=p5 || cell3[2]!=p2 )
    dserror( "triangulation failed for check5nodedAdjacentconcave" );
}

void check6nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check6nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(6);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.2;x[1] = 0.9;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.5;x[1] = 1.5;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3 || cell1[3]!=p4 )
    dserror( "triangulation failed for check6nodedconvex" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p6 || cell2[3]!=p1  )
    dserror( "triangulation failed for check6nodedconvex" );
}

void check6nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check6nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(6);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.2;x[1] = 0.9;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.5;x[1] = 0.7;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p5 || cell1[1]!=p6 || cell1[2]!=p1 || cell1[3]!=p2 )
    dserror( "triangulation failed for check6nodedconcave" );
  if( cell2[0]!=p2 || cell2[1]!=p3 || cell2[2]!=p4 || cell2[3]!=p5  )
    dserror( "triangulation failed for check6nodedconcave" );
}

void check7nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check7nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(7);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.4;x[1] = 0.6;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.1;x[1] = 0.9;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.3;x[1] = 1.5;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3 || cell1[3]!=p4 )
    dserror( "triangulation failed for check7nodedconvex" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p1  )
    dserror( "triangulation failed for check7nodedconvex" );
  if( cell3[0]!=p5 || cell3[1]!=p6 || cell3[2]!=p7 || cell3[3]!=p1 )
    dserror( "triangulation failed for check7nodedconvex" );
}

void check7nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check7nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(7);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.2;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.4;x[1] = 0.6;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.1;x[1] = 0.9;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.3;x[1] = 0.7;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.1;x[1] = 1.0;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p6 || cell1[1]!=p7 || cell1[2]!=p1 || cell1[3]!=p2 )
    dserror( "triangulation failed for check7nodedconcave" );
  if( cell2[0]!=p2 || cell2[1]!=p3 || cell2[2]!=p6 )
    dserror( "triangulation failed for check7nodedconcave" );
  if( cell3[0]!=p3 || cell3[1]!=p4 || cell3[2]!=p5 || cell3[3]!=p6 )
    dserror( "triangulation failed for check7nodedconcave" );
}

void check8nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check8nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(8);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.5;x[1] = -0.4;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = -0.05;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.1;x[1] = 0.4;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.8;x[1] = 0.8;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = -0.1;x[1] = 0.1;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3 || cell1[3]!=p4 )
    dserror( "triangulation failed for check8nodedconvex" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p6 || cell2[3]!=p1  )
    dserror( "triangulation failed for check8nodedconvex" );
  if( cell3[0]!=p6 || cell3[1]!=p7 || cell3[2]!=p8 || cell3[3]!=p1 )
    dserror( "triangulation failed for check8nodedconvex" );
}

void check8nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check8nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(8);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.5;x[1] = -0.4;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = -0.05;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.1;x[1] = 0.4;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.8;x[1] = 0.8;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.4;x[1] = 0.1;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p8 || cell1[1]!=p1 || cell1[2]!=p2 || cell1[3]!=p3 )
    dserror( "triangulation failed for check8nodedconcave" );
  if( cell2[0]!=p3 || cell2[1]!=p4 || cell2[2]!=p5 || cell2[3]!=p8 )
    dserror( "triangulation failed for check8nodedconcave" );
  if( cell3[0]!=p5 || cell3[1]!=p6 || cell3[2]!=p7 || cell3[3]!=p8 )
    dserror( "triangulation failed for check8nodedconcave" );
}

void check8nodedAdjacentconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check8nodedAdjacentconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(8);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 2.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 2.0;x[1] = 0.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 3.0;x[1] = 0.0;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 3.0;x[1] = 2.0;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.0;x[1] = 2.0;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  if( cell1[0]!=p8 || cell1[1]!=p1 || cell1[2]!=p2 )
    dserror( "triangulation failed for check8nodedAdjacentconcave" );
  if( cell2[0]!=p8 || cell2[1]!=p2 || cell2[2]!=p3 )
    dserror( "triangulation failed for check8nodedAdjacentconcave" );
  if( cell3[0]!=p4 || cell3[1]!=p5 || cell3[2]!=p6 || cell3[3]!=p7 )
    dserror( "triangulation failed for check8nodedAdjacentconcave" );
  if( cell4[0]!=p4 || cell4[1]!=p7 || cell4[2]!=p8 || cell4[3]!=p3 )
    dserror( "triangulation failed for check8nodedAdjacentconcave" );
}

void check9nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check9nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(9);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.5;x[1] = -0.4;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = -0.05;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.15;x[1] = 0.3;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 1.1;x[1] = 0.4;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.8;x[1] = 0.8;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  x[0] = -0.1;x[1] = 0.1;
  GEO::CUT::Point * p9 = mesh.NewPoint( x, NULL, s );
  ptlist[8] = p9;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3 || cell1[3]!=p4 )
    dserror( "triangulation failed for check9nodedconvex" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p6 || cell2[3]!=p1  )
    dserror( "triangulation failed for check9nodedconvex" );
  if( cell3[0]!=p6 || cell3[1]!=p7 || cell3[2]!=p8 || cell3[3]!=p1 )
    dserror( "triangulation failed for check9nodedconvex" );
  if( cell4[0]!=p8 || cell4[1]!=p9 || cell4[2]!=p1 )
    dserror( "triangulation failed for check9nodedconvex" );
}

void check9nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check9nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(9);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.5;x[1] = -0.4;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = -0.05;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.15;x[1] = 0.3;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 1.1;x[1] = 0.4;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.8;x[1] = 0.8;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  x[0] = 0.1;x[1] = 0.1;
  GEO::CUT::Point * p9 = mesh.NewPoint( x, NULL, s );
  ptlist[8] = p9;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  if( cell1[0]!=p9 || cell1[1]!=p1 || cell1[2]!=p2 || cell1[3]!=p3 )
    dserror( "triangulation failed for check9nodedconcave" );
  if( cell2[0]!=p3 || cell2[1]!=p4 || cell2[2]!=p5 || cell2[3]!=p9 )
    dserror( "triangulation failed for check9nodedconcave" );
  if( cell3[0]!=p5 || cell3[1]!=p6 || cell3[2]!=p7 || cell3[3]!=p9 )
    dserror( "triangulation failed for check9nodedconcave" );
  if( cell4[0]!=p7 || cell4[1]!=p8 || cell4[2]!=p9 )
    dserror( "triangulation failed for check9nodedconcave" );
}

void check10nodedconvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check10nodedconvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(10);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.3;x[1] = 0.2;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 1.0;x[1] = 0.4;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.2;x[1] = 0.7;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  x[0] = 0.0;x[1] = 0.4;
  GEO::CUT::Point * p9 = mesh.NewPoint( x, NULL, s );
  ptlist[8] = p9;

  x[0] = -0.1;x[1] = 0.2;
  GEO::CUT::Point * p10 = mesh.NewPoint( x, NULL, s );
  ptlist[9] = p10;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  if( cell1[0]!=p1 || cell1[1]!=p2 || cell1[2]!=p3 || cell1[3]!=p4 )
    dserror( "triangulation failed for check10nodedconvex" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p6 || cell2[3]!=p1  )
    dserror( "triangulation failed for check10nodedconvex" );
  if( cell3[0]!=p6 || cell3[1]!=p7 || cell3[2]!=p8 || cell3[3]!=p1 )
    dserror( "triangulation failed for check10nodedconvex" );
  if( cell4[0]!=p8 || cell4[1]!=p9 || cell4[2]!=p10 || cell4[3]!=p1 )
    dserror( "triangulation failed for check10nodedconvex" );
}

void check10nodedconcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check10nodedconcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(10);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = -0.3;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.4;x[1] = -0.3;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 1.3;x[1] = 0.2;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 1.0;x[1] = 0.4;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.4;x[1] = 0.7;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.2;x[1] = 0.7;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  x[0] = 0.0;x[1] = 0.4;
  GEO::CUT::Point * p9 = mesh.NewPoint( x, NULL, s );
  ptlist[8] = p9;

  x[0] = 0.1;x[1] = 0.3;
  GEO::CUT::Point * p10 = mesh.NewPoint( x, NULL, s );
  ptlist[9] = p10;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  if( cell1[0]!=p10 || cell1[1]!=p1 || cell1[2]!=p2 || cell1[3]!=p3 )
    dserror( "triangulation failed for check10nodedconcave" );
  if( cell2[0]!=p3 || cell2[1]!=p4 || cell2[2]!=p5 || cell2[3]!=p10 )
    dserror( "triangulation failed for check10nodedconcave" );
  if( cell3[0]!=p5 || cell3[1]!=p6 || cell3[2]!=p7 || cell3[3]!=p10 )
    dserror( "triangulation failed for check10nodedconcave" );
  if( cell4[0]!=p7 || cell4[1]!=p8 || cell4[2]!=p9 || cell4[3]!=p10 )
    dserror( "triangulation failed for check10nodedconcave" );
}

void check5nodedTwinConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check5nodedTwinConcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(5);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.2;x[1] = 0.1;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 0.3;x[1] = 0.2;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.2;x[1] = 0.7;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

/*  std::cout<<"split size = "<<split.size()<<"\n";
  for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/


  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p2 || cell1[1]!=p3 || cell1[2]!=p4 )
    dserror( "triangulation failed for check5nodedTwinConcave" );
  if( cell2[0]!=p4 || cell2[1]!=p5 || cell2[2]!=p1 )
    dserror( "triangulation failed for check5nodedTwinConcave" );
}

void check6nodedTwinConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check6nodedTwinConcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(6);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 2.0;x[1] = 0.0;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.3;x[1] = 0.7;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 2.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.0;x[1] = 1.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.9;x[1] = 0.7;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  if( cell1[0]!=p3 || cell1[1]!=p4 || cell1[2]!=p5 || cell1[3]!=p6 )
    dserror( "triangulation failed for check6nodedTwinConcave" );
  if( cell2[0]!=p6 || cell2[1]!=p1 || cell2[2]!=p2 || cell2[3]!=p3  )
    dserror( "triangulation failed for check6nodedTwinConcave" );
}

void check8nodedTriConcave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check8nodedTriConcave...\n";
  std::vector<GEO::CUT::Point*> ptlist(8);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.5;x[1] = 0.2;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.5;x[1] = 0.7;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.0;x[1] = 1.0;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 0.0;x[1] = 0.8;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 0.2;x[1] = 0.5;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p2 || cell1[1]!=p3 || cell1[2]!=p4 || cell1[3]!=p5 )
    dserror( "triangulation failed for check8nodedTriConcave" );
  if( cell2[0]!=p5 || cell2[1]!=p6 || cell2[2]!=p7 || cell2[3]!=p8  )
    dserror( "triangulation failed for check8nodedTriConcave" );
  if( cell3[0]!=p8 || cell3[1]!=p1 || cell3[2]!=p2 )
    dserror( "triangulation failed for check8nodedTriConcave" );
}

// the points are considered in the plane 2x+3y+4z=7
void check8nodedTriConcaveGenPlane( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check8nodedTriConcaveGenPlane...\n";
  std::vector<GEO::CUT::Point*> ptlist(8);
  double x[3];

  x[1] = 0.0;x[2] = 0.0;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[1] = 0.5;x[2] = 0.2;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[1] = 1.0;x[2] = 0.0;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[1] = 1.0;x[2] = 1.0;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[1] = 0.5;x[2] = 0.7;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[1] = 0.0;x[2] = 1.0;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[1] = 0.0;x[2] = 0.8;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[1] = 0.2;x[2] = 0.5;x[0]=0.5*(7.0-3.0*x[1]-4.0*x[2]);
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  if( cell1[0]!=p2 || cell1[1]!=p3 || cell1[2]!=p4 || cell1[3]!=p5 )
    dserror( "triangulation failed for check8nodedTriConcaveGenPlane" );
  if( cell2[0]!=p5 || cell2[1]!=p6 || cell2[2]!=p7 || cell2[3]!=p8  )
    dserror( "triangulation failed for check8nodedTriConcaveGenPlane" );
  if( cell3[0]!=p8 || cell3[1]!=p1 || cell3[2]!=p2 )
    dserror( "triangulation failed for check8nodedTriConcaveGenPlane" );
}

void check13nodedConvex( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check13nodedConvex...\n";
  std::vector<GEO::CUT::Point*> ptlist(13);
  double x[3];
  x[2] = 0.0;

  //std::cout<<"the points are  \n";
  for( int i=0;i<13;i++ )
  {
    double theta = 2.0*22.0/7.0/13*i;
    x[0] = cos(theta);
    x[1] = sin(theta);

    GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
    ptlist[i] = p1;
    //std::cout<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\n";
  }

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

/*  for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  if( split.size()!=6 )
    dserror( "triangulation failed for check13nodedConvex" );

  for( int i=0;i<6;i++ )
  {
    std::vector<GEO::CUT::Point*> cell = split[i];
    if( i==5 )
    {
      if( cell[0]!=ptlist[0] || cell[1]!=ptlist[i*2+1] || cell[2]!=ptlist[i*2+2] )
        dserror( "triangulation failed for check13nodedConvex" );
    }
    else
    {
      if( cell[0]!=ptlist[0] || cell[1]!=ptlist[i*2+1] || cell[2]!=ptlist[i*2+2] || cell[3]!=ptlist[i*2+3] )
        dserror( "triangulation failed for check13nodedConvex" );
    }
  }
}

void check7nodedconti3concave( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check7nodedconti3concave...\n";
  std::vector<GEO::CUT::Point*> ptlist(7);
  double x[3];
  x[2] = 0.0;

  x[0] = 1.0;x[1] = 0.728885;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 0.992803;x[1] = 0.726274;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 0.958516;x[1] = 0.710801;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 0.950656;x[1] = 0.706818;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 0.916667;x[1] = 0.671407;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 0.916667;x[1] = 0.75;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 1.0;x[1] = 0.75;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

/*  for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  std::vector<GEO::CUT::Point*> cell5 = split[4];
  if( cell1[0]!=p7 || cell1[1]!=p1 || cell1[2]!=p2 )
    dserror( "triangulation failed for check7nodedconti3concave" );
  if( cell2[0]!=p7 || cell2[1]!=p2 || cell2[2]!=p3 )
    dserror( "triangulation failed for check7nodedconti3concave" );
  if( cell3[0]!=p4 || cell3[1]!=p5 || cell3[2]!=p6 )
    dserror( "triangulation failed for check7nodedconti3concave" );
  if( cell4[0]!=p4 || cell4[1]!=p6 || cell4[2]!=p7 )
    dserror( "triangulation failed for check7nodedconti3concave" );
  if( cell5[0]!=p4 || cell5[1]!=p7 || cell5[2]!=p3 )
    dserror( "triangulation failed for check7nodedconti3concave" );
}

void check10nodedShiftEarClipToSplit( GEO::CUT::Mesh& mesh, GEO::CUT::Element * e, GEO::CUT::Side * s )
{
  std::cout<<"check10nodedShiftEarClipToSplit...\n";
  std::vector<GEO::CUT::Point*> ptlist(10);
  double x[3];
  x[2] = 0.0;

  x[0] = 0.0;x[1] = 0.0;
  GEO::CUT::Point * p1 = mesh.NewPoint( x, NULL, s );
  ptlist[0] = p1;

  x[0] = 1.0;x[1] = 0.0;
  GEO::CUT::Point * p2 = mesh.NewPoint( x, NULL, s );
  ptlist[1] = p2;

  x[0] = 1.0;x[1] = 1.0;
  GEO::CUT::Point * p3 = mesh.NewPoint( x, NULL, s );
  ptlist[2] = p3;

  x[0] = 2.0;x[1] = 1.0;
  GEO::CUT::Point * p4 = mesh.NewPoint( x, NULL, s );
  ptlist[3] = p4;

  x[0] = 2.0;x[1] = 0.0;
  GEO::CUT::Point * p5 = mesh.NewPoint( x, NULL, s );
  ptlist[4] = p5;

  x[0] = 3.0;x[1] = 0.0;
  GEO::CUT::Point * p6 = mesh.NewPoint( x, NULL, s );
  ptlist[5] = p6;

  x[0] = 2.5;x[1] = 1.0;
  GEO::CUT::Point * p7 = mesh.NewPoint( x, NULL, s );
  ptlist[6] = p7;

  x[0] = 3.0;x[1] = 2.0;
  GEO::CUT::Point * p8 = mesh.NewPoint( x, NULL, s );
  ptlist[7] = p8;

  x[0] = 1.5;x[1] = 1.5;
  GEO::CUT::Point * p9 = mesh.NewPoint( x, NULL, s );
  ptlist[8] = p9;

  x[0] = 0.0;x[1] = 2.0;
  GEO::CUT::Point * p10 = mesh.NewPoint( x, NULL, s );
  ptlist[9] = p10;

  GEO::CUT::Facet face1( mesh, ptlist, s, false );
  GEO::CUT::TriangulateFacet tf( &face1, mesh, ptlist );
  tf.SplitFacet();

  std::vector<std::vector<GEO::CUT::Point*> > split;
  split = tf.GetSplitCells();

  /*for( std::vector<std::vector<GEO::CUT::Point*> >::iterator i=split.begin();i!=split.end();i++ )
  {
    std::cout<<"cell\n";
    std::vector<GEO::CUT::Point*> cell = *i;;
    for( std::vector<GEO::CUT::Point*>::iterator j=cell.begin();j!=cell.end();j++ )
    {
      GEO::CUT::Point* pt = *j;
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/

  std::vector<GEO::CUT::Point*> cell1 = split[0];
  std::vector<GEO::CUT::Point*> cell2 = split[1];
  std::vector<GEO::CUT::Point*> cell3 = split[2];
  std::vector<GEO::CUT::Point*> cell4 = split[3];
  std::vector<GEO::CUT::Point*> cell5 = split[4];
  if( cell1[0]!=p10 || cell1[1]!=p1 || cell1[2]!=p2 )
    dserror( "triangulation failed for check10nodedShiftEarClipToSplit" );
  if( cell2[0]!=p10 || cell2[1]!=p2 || cell2[2]!=p3 )
    dserror( "triangulation failed for check10nodedShiftEarClipToSplit" );
  if( cell3[0]!=p4 || cell3[1]!=p5 || cell3[2]!=p6 || cell3[3]!=p7 )
    dserror( "triangulation failed for check10nodedShiftEarClipToSplit" );
  if( cell4[0]!=p7 || cell4[1]!=p8 || cell4[2]!=p9 )
    dserror( "triangulation failed for check10nodedShiftEarClipToSplit" );
  if( cell5[0]!=p9 || cell5[1]!=p10 || cell5[2]!=p3 || cell5[3]!=p4 )
    dserror( "triangulation failed for check10nodedShiftEarClipToSplit" );
}
