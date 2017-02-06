/*---------------------------------------------------------------------------*/
/*!
\file cut_test_intersection.cpp

\brief cut test cpp file

\level 1

\maintainer Benedikt Schott, Christoph Ager

*/
/*---------------------------------------------------------------------------*/

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_intersection.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

#define nxyz( x, y, z ) { if ( nodeids.count( #x#y#z )==0 ) { nodeids[#x#y#z] = nodeids.size(); }}
#define nxyz1( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 0 )=x; quad4_xyze( 1, 0 )=y; quad4_xyze( 2, 0 )=z; nids[0] = nodeids[#x#y#z];
#define nxyz2( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 1 )=x; quad4_xyze( 1, 1 )=y; quad4_xyze( 2, 1 )=z; nids[1] = nodeids[#x#y#z];
#define nxyz3( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 2 )=x; quad4_xyze( 1, 2 )=y; quad4_xyze( 2, 2 )=z; nids[2] = nodeids[#x#y#z];
#define nxyz4( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 3 )=x; quad4_xyze( 1, 3 )=y; quad4_xyze( 2, 3 )=z; nids[3] = nodeids[#x#y#z];


void test_quad4_line2( double x1, double y1,
                       double x2, double y2 )
{
#if 0
  std::cout << "(" << x1 << "," << y1 << ")"
            << "--"
            << "(" << x2 << "," << y2 << ")\n";
#endif

  GEO::CUT::Options options;
  GEO::CUT::Mesh mesh(options);

  std::vector<double> xyz( 3 );

  xyz[0] = x1;
  xyz[1] = y1;
  xyz[2] = 0;
  GEO::CUT::Node * n1 = mesh.GetNode( 1, &xyz[0] );

  xyz[0] = x2;
  xyz[1] = y2;
  GEO::CUT::Node * n2 = mesh.GetNode( 2, &xyz[0] );

  GEO::CUT::Edge * edge = mesh.GetEdge( n1, n2 );

  std::vector<GEO::CUT::Node*> nodes;
  xyz[0] = 0;
  xyz[1] = 0;
  nodes.push_back( mesh.GetNode( 3, &xyz[0] ) );
  xyz[0] = 1;
  xyz[1] = 0;
  nodes.push_back( mesh.GetNode( 4, &xyz[0] ) );
  xyz[0] = 1;
  xyz[1] = 1;
  nodes.push_back( mesh.GetNode( 5, &xyz[0] ) );
  xyz[0] = 0;
  xyz[1] = 1;
  nodes.push_back( mesh.GetNode( 6, &xyz[0] ) );

  GEO::CUT::Side * side = mesh.GetSide( 1, nodes, shards::getCellTopologyData< shards::Quadrilateral<4> >() );

  Teuchos::RCP<GEO::CUT::IntersectionBase> intersection =
      GEO::CUT::IntersectionBase::Create( DRT::Element::line2, DRT::Element::quad4 );
  intersection->Init( & mesh, edge, side, false, false, false );

#if 0
  GEO::CUT::PointSet cuts;
  bool does = intersection.Intersect( cuts );
  std::cout << "does intersect: " << does << "  " << cuts.size() << "\n";
  for ( GEO::CUT::PointSet::iterator i=cuts.begin(); i!=cuts.end(); ++i )
  {
    GEO::CUT::Point * p = *i;
    p->Plot( std::cout );
  }
#else
  GEO::CUT::PointSet cuts;
  intersection->Intersect( cuts );
#endif

  mesh.Status();

  if ( cuts.size()!=2 )
  {
    std::stringstream str;
    str << "two cuts expected but got " << cuts.size() << ": ";
    std::copy( cuts.begin(), cuts.end(), std::ostream_iterator<GEO::CUT::Point*>( str, " " ) );
    throw std::runtime_error( str.str() );
  }
}

void test_quad4_line2()
{
  test_quad4_line2(  0.5,  1,   1.5, 1   );
  test_quad4_line2(  0.5,  0.5, 1.5, 0.5 );
  test_quad4_line2( -0.5,  0.5, 1.5, 0.5 );
  test_quad4_line2( -0.5,  0  , 1.5, 0   );
  test_quad4_line2( -0.5, -0.5, 1.5, 1.5 );
}

void test_hex8_quad4_qhull1()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.821504414, 0.524610698, 0.589920461);
  nxyz2( 0.842310429, 0.4672589, 0.58473295);
  nxyz3( 0.85763216, 0.4672589, 0.646184981);
  nxyz4( 0.836826086, 0.524610698, 0.651372552);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.806182683, 0.524610698, 0.52846843);
  nxyz2( 0.826988697, 0.4672589, 0.523280859);
  nxyz3( 0.842310429, 0.4672589, 0.58473295);
  nxyz4( 0.821504414, 0.524610698, 0.589920461);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.75;
  hex8_xyze(1,0) = 0.416666657;
  hex8_xyze(2,0) = 0.666666687;
  hex8_xyze(0,1) = 0.75;
  hex8_xyze(1,1) = 0.416666657;
  hex8_xyze(2,1) = 0.583333313;
  hex8_xyze(0,2) = 0.75;
  hex8_xyze(1,2) = 0.5;
  hex8_xyze(2,2) = 0.583333313;
  hex8_xyze(0,3) = 0.75;
  hex8_xyze(1,3) = 0.5;
  hex8_xyze(2,3) = 0.666666687;
  hex8_xyze(0,4) = 0.833333313;
  hex8_xyze(1,4) = 0.416666657;
  hex8_xyze(2,4) = 0.666666687;
  hex8_xyze(0,5) = 0.833333313;
  hex8_xyze(1,5) = 0.416666657;
  hex8_xyze(2,5) = 0.583333313;
  hex8_xyze(0,6) = 0.833333313;
  hex8_xyze(1,6) = 0.5;
  hex8_xyze(2,6) = 0.583333313;
  hex8_xyze(0,7) = 0.833333313;
  hex8_xyze(1,7) = 0.5;
  hex8_xyze(2,7) = 0.666666687;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

//INFO: This Cut test will not work with DD in LOCAL, as the hex8 does not fullfill the restrictions of this method!
void test_hex8_quad4_alex1()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.69467040756725095818,0.4021667090047654769,-0.025250000899999999748);
  nxyz2(0.70968393530890683252,0.41573627765474635565,-0.025250000899999999748);
  nxyz3(0.70968393546204655564,0.41573627770040033624,0.025250000899999999748);
  nxyz4(0.69467040760366560725,0.40216670904568263545,0.025250000899999999748);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.70968393530890683252,0.41573627765474635565,-0.025250000899999999748);
  nxyz2(0.72516543896079255216,0.42926679214049395794,-0.025250000899999999748);
  nxyz3(0.72516543893872920101,0.42926679226635205966,0.025250000899999999748);
  nxyz4(0.70968393546204655564,0.41573627770040033624,0.025250000899999999748);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.72516543896079255216,0.42926679214049395794,-0.025250000899999999748);
  nxyz2(0.74086661529827402273,0.44309748664700010501,-0.025250000899999999748);
  nxyz3(0.74086661523377261851,0.4430974866683990987,0.025250000899999999748);
  nxyz4(0.72516543893872920101,0.42926679226635205966,0.025250000899999999748);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.70729452400000003554;
  hex8_xyze(1,0) = 0.39701718099999999678;
  hex8_xyze(2,0) = -0.025000000399999999789;
  hex8_xyze(0,1) = 0.7415404319999999716;
  hex8_xyze(1,1) = 0.39539372900000002709;
  hex8_xyze(2,1) = -0.025000000399999999789;
  hex8_xyze(0,2) = 0.73677337200000003747;
  hex8_xyze(1,2) = 0.42250511099999998876;
  hex8_xyze(2,2) = -0.025000000399999999789;
  hex8_xyze(0,3) = 0.70965522500000000061;
  hex8_xyze(1,3) = 0.43864250199999998969;
  hex8_xyze(2,3) = -0.025000000399999999789;
  hex8_xyze(0,4) = 0.70729452400000003554;
  hex8_xyze(1,4) = 0.39701718099999999678;
  hex8_xyze(2,4) = 0.025000000399999999789;
  hex8_xyze(0,5) = 0.7415404319999999716;
  hex8_xyze(1,5) = 0.39539372900000002709;
  hex8_xyze(2,5) = 0.025000000399999999789;
  hex8_xyze(0,6) = 0.73677337200000003747;
  hex8_xyze(1,6) = 0.42250511099999998876;
  hex8_xyze(2,6) = 0.025000000399999999789;
  hex8_xyze(0,7) = 0.70965522500000000061;
  hex8_xyze(1,7) = 0.43864250199999998969;
  hex8_xyze(2,7) = 0.025000000399999999789;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex2()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(1.159049510000000005,0.081087201799999994445,-0.025250000899999999748);
  nxyz2(1.1400269300000001049,0.08774351330000000615,-0.025250000899999999748);
  nxyz3(1.1400269300000001049,0.08774351330000000615,0.025250000899999999748);
  nxyz4(1.159049510000000005,0.081087201799999994445,0.025250000899999999748);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.1400269300000001049,0.08774351330000000615,-0.025250000899999999748);
  nxyz2(1.1200000000000001066,0.090000003600000003101,-0.025250000899999999748);
  nxyz3(1.1200000000000001066,0.090000003600000003101,0.025250000899999999748);
  nxyz4(1.1400269300000001049,0.08774351330000000615,0.025250000899999999748);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.1200000000000001066,0.10999999900000000108,-0.025250000899999999748);
  nxyz2(1.1400269300000001049,0.11225649000000000044,-0.025250000899999999748);
  nxyz3(1.1400269300000001049,0.11225649000000000044,0.025250000899999999748);
  nxyz4(1.1200000000000001066,0.10999999900000000108,0.025250000899999999748);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.1400269300000001049,0.11225649000000000044,-0.025250000899999999748);
  nxyz2(1.159049510000000005,0.11891280099999999853,-0.025250000899999999748);
  nxyz3(1.159049510000000005,0.11891280099999999853,0.025250000899999999748);
  nxyz4(1.1400269300000001049,0.11225649000000000044,0.025250000899999999748);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1.16999995999999995;
  hex8_xyze(1,0) = 0.11428571499999999606;
  hex8_xyze(2,0) = 0.025000000399999999789;
  hex8_xyze(0,1) = 1.1399999899999999631;
  hex8_xyze(1,1) = 0.11428571499999999606;
  hex8_xyze(2,1) = 0.025000000399999999789;
  hex8_xyze(0,2) = 1.1399999899999999631;
  hex8_xyze(1,2) = 0.11428571499999999606;
  hex8_xyze(2,2) = -0.025000000399999999789;
  hex8_xyze(0,3) = 1.16999995999999995;
  hex8_xyze(1,3) = 0.11428571499999999606;
  hex8_xyze(2,3) = -0.025000000399999999789;
  hex8_xyze(0,4) = 1.16999995999999995;
  hex8_xyze(1,4) = 0.085714288099999993986;
  hex8_xyze(2,4) = 0.025000000399999999789;
  hex8_xyze(0,5) = 1.1399999899999999631;
  hex8_xyze(1,5) = 0.085714288099999993986;
  hex8_xyze(2,5) = 0.025000000399999999789;
  hex8_xyze(0,6) = 1.1399999899999999631;
  hex8_xyze(1,6) = 0.085714288099999993986;
  hex8_xyze(2,6) = -0.025000000399999999789;
  hex8_xyze(0,7) = 1.16999995999999995;
  hex8_xyze(1,7) = 0.085714288099999993986;
  hex8_xyze(2,7) = -0.025000000399999999789;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex3()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.85005472577752938523,0.01745303543601192317,0.31893504400781552688);
  nxyz2(0.87011280475064345374,0.017412162235328237925,0.31899942251609070265);
  nxyz3(0.87017171851365637814,0.036824462160089978247,0.31899931623323690699);
  nxyz4(0.85011194245098542499,0.036917863796646292751,0.31892277220611009447);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.85007500599999996638,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.87015003000000001876,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.87011280475064345374,0.017412162235328237925,0.31899942251609070265);
  nxyz4(0.85005472577752938523,0.01745303543601192317,0.31893504400781552688);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89022499300000001998,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.91030001599999998962,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.91033437387081139747,0.017328954203733981654,0.31912150736715011456);
  nxyz4(0.89020365758509467646,0.017371082628566521244,0.31906397688495774512);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89020365758509467646,0.017371082628566521244,0.31906397688495774512);
  nxyz2(0.91033437387081139747,0.017328954203733981654,0.31912150736715011456);
  nxyz3(0.91036696497657321192,0.036638768964142386098,0.31915135370820452154);
  nxyz4(0.8902602101632242082,0.036730974597990441455,0.31907598394014535792);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.87011280475064345374,0.017412162235328237925,0.31899942251609070265);
  nxyz2(0.89020365758509467646,0.017371082628566521244,0.31906397688495774512);
  nxyz3(0.8902602101632242082,0.036730974597990441455,0.31907598394014535792);
  nxyz4(0.87017171851365637814,0.036824462160089978247,0.31899931623323690699);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.87015003000000001876,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.89022499300000001998,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.89020365758509467646,0.017371082628566521244,0.31906397688495774512);
  nxyz4(0.87011280475064345374,0.017412162235328237925,0.31899942251609070265);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.86956518900000001615;
  hex8_xyze(1,0) = 0;
  hex8_xyze(2,0) = 0.30000001199999998258;
  hex8_xyze(0,1) = 0.89130437399999995485;
  hex8_xyze(1,1) = 0;
  hex8_xyze(2,1) = 0.30000001199999998258;
  hex8_xyze(0,2) = 0.89130437399999995485;
  hex8_xyze(1,2) = 0.022727273400000001097;
  hex8_xyze(2,2) = 0.30000001199999998258;
  hex8_xyze(0,3) = 0.86956518900000001615;
  hex8_xyze(1,3) = 0.022727273400000001097;
  hex8_xyze(2,3) = 0.30000001199999998258;
  hex8_xyze(0,4) = 0.86956518900000001615;
  hex8_xyze(1,4) = 0;
  hex8_xyze(2,4) = 0.34999999399999998095;
  hex8_xyze(0,5) = 0.89130437399999995485;
  hex8_xyze(1,5) = 0;
  hex8_xyze(2,5) = 0.34999999399999998095;
  hex8_xyze(0,6) = 0.89130437399999995485;
  hex8_xyze(1,6) = 0.022727273400000001097;
  hex8_xyze(2,6) = 0.34999999399999998095;
  hex8_xyze(0,7) = 0.86956518900000001615;
  hex8_xyze(1,7) = 0.022727273400000001097;
  hex8_xyze(2,7) = 0.34999999399999998095;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex4()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.85007500599999996638,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.87015003000000001876,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.87013650387801810826,0.017411815601464999959,0.31899986664172991224);
  nxyz4(0.85006757708168168008,0.017426803601911681346,0.3189762725947100086);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.85006757708168168008,0.017426803601911681346,0.3189762725947100086);
  nxyz2(0.87013650387801810826,0.017411815601464999959,0.31899986664172991224);
  nxyz3(0.87015820344018968147,0.036823660065671247332,0.31899984417057236641);
  nxyz4(0.8500885937777643564,0.036857880365897197072,0.31897182314140559711);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89021738600324706869,0.017396800949507920664,0.31902348251218615172);
  nxyz2(0.91031287743832123294,0.017381346031252110917,0.31904449943988899552);
  nxyz3(0.91032508978127524291,0.036755638330250986479,0.31905536627092928592);
  nxyz4(0.8902383031165566063,0.036789440073439365342,0.31902787500164930812);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89022499300000001998,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.91030001599999998962,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.91031287743832123294,0.017381346031252110917,0.31904449943988899552);
  nxyz4(0.89021738600324706869,0.017396800949507920664,0.31902348251218615172);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.87013650387801810826,0.017411815601464999959,0.31899986664172991224);
  nxyz2(0.89021738600324706869,0.017396800949507920664,0.31902348251218615172);
  nxyz3(0.8902383031165566063,0.036789440073439365342,0.31902787500164930812);
  nxyz4(0.87015820344018968147,0.036823660065671247332,0.31899984417057236641);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.87015003000000001876,-0.0020000000900000001157,0.31900000600000000262);
  nxyz2(0.89022499300000001998,-0.0020000000900000001157,0.31900000600000000262);
  nxyz3(0.89021738600324706869,0.017396800949507920664,0.31902348251218615172);
  nxyz4(0.87013650387801810826,0.017411815601464999959,0.31899986664172991224);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.86956518900000001615;
  hex8_xyze(1,0) = 0;
  hex8_xyze(2,0) = 0.30000001199999998258;
  hex8_xyze(0,1) = 0.89130437399999995485;
  hex8_xyze(1,1) = 0;
  hex8_xyze(2,1) = 0.30000001199999998258;
  hex8_xyze(0,2) = 0.89130437399999995485;
  hex8_xyze(1,2) = 0.022727273400000001097;
  hex8_xyze(2,2) = 0.30000001199999998258;
  hex8_xyze(0,3) = 0.86956518900000001615;
  hex8_xyze(1,3) = 0.022727273400000001097;
  hex8_xyze(2,3) = 0.30000001199999998258;
  hex8_xyze(0,4) = 0.86956518900000001615;
  hex8_xyze(1,4) = 0;
  hex8_xyze(2,4) = 0.34999999399999998095;
  hex8_xyze(0,5) = 0.89130437399999995485;
  hex8_xyze(1,5) = 0;
  hex8_xyze(2,5) = 0.34999999399999998095;
  hex8_xyze(0,6) = 0.89130437399999995485;
  hex8_xyze(1,6) = 0.022727273400000001097;
  hex8_xyze(2,6) = 0.34999999399999998095;
  hex8_xyze(0,7) = 0.86956518900000001615;
  hex8_xyze(1,7) = 0.022727273400000001097;
  hex8_xyze(2,7) = 0.34999999399999998095;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
  intersection.PrintCellStats();
}

void test_hex8_quad4_alex5()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.84430031334490374118,0.25029534684190529736,0.080176760981080896773);
  nxyz2(0.8513578595329502896,0.33261253137199714436,0.080186838347331437782);
  nxyz3(0.85130492844640159866,0.33260358891054925268,-0.00010000000000000000479);
  nxyz4(0.84426747403350743681,0.2502909361877762584,-0.00010000000000000000479);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.85130492844640159866,0.33260358891054925268,-0.00010000000000000000479);
  nxyz2(0.8513578595329502896,0.33261253137199714436,0.080186838347331437782);
  nxyz3(0.93136053102012328342,0.32577117056017301788,0.080116848605698884334);
  nxyz4(0.93130992716486038496,0.32578977984609908125,-0.00010000000000000000479);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.8513578595329502896,0.33261253137199714436,0.080186838347331437782);
  nxyz2(0.85147921498170020538,0.33264274398766413565,0.16043979097814206614);
  nxyz3(0.93147585299191071329,0.32571895020837232648,0.1603653495145052621);
  nxyz4(0.93136053102012328342,0.32577117056017301788,0.080116848605698884334);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.84435827369202076653,0.25031356557852207922,0.16041800389897961643);
  nxyz2(0.85147921498170020538,0.33264274398766413565,0.16043979097814206614);
  nxyz3(0.8513578595329502896,0.33261253137199714436,0.080186838347331437782);
  nxyz4(0.84430031334490374118,0.25029534684190529736,0.080176760981080896773);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.83333333333333348136;
  hex8_xyze(1,0) = 0.32000000000000000666;
  hex8_xyze(2,0) = 0.055555555555555427572;
  hex8_xyze(0,1) = 0.86111111111111116045;
  hex8_xyze(1,1) = 0.32000000000000000666;
  hex8_xyze(2,1) = 0.055555555555555427572;
  hex8_xyze(0,2) = 0.86111111111111127148;
  hex8_xyze(1,2) = 0.34666666666666667851;
  hex8_xyze(2,2) = 0.055555555555555427572;
  hex8_xyze(0,3) = 0.83333333333333348136;
  hex8_xyze(1,3) = 0.34666666666666667851;
  hex8_xyze(2,3) = 0.055555555555555427572;
  hex8_xyze(0,4) = 0.83333333333333337034;
  hex8_xyze(1,4) = 0.32000000000000000666;
  hex8_xyze(2,4) = 0.11111111111111085514;
  hex8_xyze(0,5) = 0.86111111111111116045;
  hex8_xyze(1,5) = 0.32000000000000000666;
  hex8_xyze(2,5) = 0.11111111111111085514;
  hex8_xyze(0,6) = 0.86111111111111116045;
  hex8_xyze(1,6) = 0.34666666666666667851;
  hex8_xyze(2,6) = 0.11111111111111085514;
  hex8_xyze(0,7) = 0.83333333333333337034;
  hex8_xyze(1,7) = 0.34666666666666667851;
  hex8_xyze(2,7) = 0.11111111111111085514;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex6()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.83,0.0824,0.3209);
  nxyz2(0.9103,0.0824,0.3209);
  nxyz3(0.9103,0.1649,0.3209);
  nxyz4(0.83,0.1649,0.3209);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.0824,0.3209);
  nxyz2(0.9103,0.0824,0.24065);
  nxyz3(0.9103,0.1649,0.24065);
  nxyz4(0.9103,0.1649,0.3209);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.83,0.1649,0.3209);
  nxyz2(0.9103,0.1649,0.3209);
  nxyz3(0.9103,0.2474,0.3209);
  nxyz4(0.83,0.2474,0.3209);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.1649,0.3209);
  nxyz2(0.9103,0.1649,0.24065);
  nxyz3(0.9103,0.2474,0.24065);
  nxyz4(0.9103,0.2474,0.3209);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.888889;
  hex8_xyze(1,0) = 0.16;
  hex8_xyze(2,0) = 0.277778;
  hex8_xyze(0,1) = 0.916667;
  hex8_xyze(1,1) = 0.16;
  hex8_xyze(2,1) = 0.277778;
  hex8_xyze(0,2) = 0.916667;
  hex8_xyze(1,2) = 0.186667;
  hex8_xyze(2,2) = 0.277778;
  hex8_xyze(0,3) = 0.888889;
  hex8_xyze(1,3) = 0.186667;
  hex8_xyze(2,3) = 0.277778;
  hex8_xyze(0,4) = 0.888889;
  hex8_xyze(1,4) = 0.16;
  hex8_xyze(2,4) = 0.333333;
  hex8_xyze(0,5) = 0.916667;
  hex8_xyze(1,5) = 0.16;
  hex8_xyze(2,5) = 0.333333;
  hex8_xyze(0,6) = 0.916667;
  hex8_xyze(1,6) = 0.186667;
  hex8_xyze(2,6) = 0.333333;
  hex8_xyze(0,7) = 0.888889;
  hex8_xyze(1,7) = 0.186667;
  hex8_xyze(2,7) = 0.333333;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex7()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.832,0.716,0.32);
  nxyz2(0.912,0.719,0.321);
  nxyz3(0.91,0.8,0.321);
  nxyz4(0.83,0.8,0.321);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.912,0.719,0.321);
  nxyz2(0.913,0.719,0.241);
  nxyz3(0.91,0.8,0.241);
  nxyz4(0.91,0.8,0.321);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.917,0.638,0.321);
  nxyz2(0.918,0.638,0.241);
  nxyz3(0.913,0.719,0.241);
  nxyz4(0.912,0.719,0.321);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.837,0.632,0.32);
  nxyz2(0.917,0.638,0.321);
  nxyz3(0.912,0.719,0.321);
  nxyz4(0.832,0.716,0.32);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.889;
  hex8_xyze(1,0) = 0.693;
  hex8_xyze(2,0) = 0.278;
  hex8_xyze(0,1) = 0.917;
  hex8_xyze(1,1) = 0.693;
  hex8_xyze(2,1) = 0.278;
  hex8_xyze(0,2) = 0.917;
  hex8_xyze(1,2) = 0.72;
  hex8_xyze(2,2) = 0.278;
  hex8_xyze(0,3) = 0.889;
  hex8_xyze(1,3) = 0.72;
  hex8_xyze(2,3) = 0.278;
  hex8_xyze(0,4) = 0.889;
  hex8_xyze(1,4) = 0.693;
  hex8_xyze(2,4) = 0.333;
  hex8_xyze(0,5) = 0.917;
  hex8_xyze(1,5) = 0.693;
  hex8_xyze(2,5) = 0.333;
  hex8_xyze(0,6) = 0.917;
  hex8_xyze(1,6) = 0.72;
  hex8_xyze(2,6) = 0.333;
  hex8_xyze(0,7) = 0.889;
  hex8_xyze(1,7) = 0.72;
  hex8_xyze(2,7) = 0.333;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex8()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.83717493470789061671,0.16765129807116907479,0.32048082642379321916);
  nxyz2(0.91725380481827356594,0.16190652493710860904,0.32130926372145246495);
  nxyz3(0.92409481657760572659,0.24358544117829519782,0.32111191840552105736);
  nxyz4(0.84409422239187015258,0.25043070880620033059,0.32068621380180872826);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91725380481827356594,0.16190652493710860904,0.32130926372145246495);
  nxyz2(0.91773307135174408966,0.16188415202742093002,0.24078373402031172379);
  nxyz3(0.92434065980863056033,0.24360126684155039567,0.24070854672435987309);
  nxyz4(0.92409481657760572659,0.24358544117829519782,0.32111191840552105736);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91217166751135125846,0.080675690213319439792,0.32142094357479067401);
  nxyz2(0.91267024464855073251,0.080630229657985358349,0.24082622257457167447);
  nxyz3(0.91773307135174408966,0.16188415202742093002,0.24078373402031172379);
  nxyz4(0.91725380481827356594,0.16190652493710860904,0.32130926372145246495);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.83195664187558648894,0.084150641456315830591,0.3203718853205224093);
  nxyz2(0.91217166751135125846,0.080675690213319439792,0.32142094357479067401);
  nxyz3(0.91725380481827356594,0.16190652493710860904,0.32130926372145246495);
  nxyz4(0.83717493470789061671,0.16765129807116907479,0.32048082642379321916);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.94444444444444441977;
  hex8_xyze(1,0) = 0.16000000000000014211;
  hex8_xyze(2,0) = 0.27777777777777834523;
  hex8_xyze(0,1) = 0.94444444444444441977;
  hex8_xyze(1,1) = 0.18666666666666684171;
  hex8_xyze(2,1) = 0.27777777777777834523;
  hex8_xyze(0,2) = 0.91666666666666662966;
  hex8_xyze(1,2) = 0.18666666666666684171;
  hex8_xyze(2,2) = 0.27777777777777834523;
  hex8_xyze(0,3) = 0.91666666666666662966;
  hex8_xyze(1,3) = 0.16000000000000014211;
  hex8_xyze(2,3) = 0.27777777777777834523;
  hex8_xyze(0,4) = 0.94444444444444430875;
  hex8_xyze(1,4) = 0.16000000000000016986;
  hex8_xyze(2,4) = 0.33333333333333331483;
  hex8_xyze(0,5) = 0.94444444444444430875;
  hex8_xyze(1,5) = 0.18666666666666684171;
  hex8_xyze(2,5) = 0.33333333333333331483;
  hex8_xyze(0,6) = 0.91666666666666651864;
  hex8_xyze(1,6) = 0.18666666666666684171;
  hex8_xyze(2,6) = 0.33333333333333331483;
  hex8_xyze(0,7) = 0.91666666666666651864;
  hex8_xyze(1,7) = 0.16000000000000016986;
  hex8_xyze(2,7) = 0.33333333333333331483;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_tet4_quad4_alex9()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(1.1373961924174400551,0.43783278788480001209,0.015000075099999999467);
  nxyz2(1.1372558011327500438,0.43786377186006097961,1.3891340300000000002e-19);
  nxyz3(1.1633039901419599538,0.51827976384601803783,0);
  nxyz4(1.1633580248937700485,0.51825684289835904917,0.015000075099999999467);

  GEO::CUT::SideHandle * s1 = intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.1372558011327500438,0.43786377186006097961,1.3891340300000000002e-19);
  nxyz2(1.1373961924222799613,0.43783278790287999405,-0.015000075099999999467);
  nxyz3(1.1633580249045900601,0.5182568428569539476,-0.015000075099999999467);
  nxyz4(1.1633039901419599538,0.51827976384601803783,0);

  GEO::CUT::SideHandle * s2 = intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.0863570485981299818,0.54164967699933597167,5.1690185600000003046e-19);
  nxyz2(1.0862518404856800203,0.54179921727229496398,0.015000075099999999467);
  nxyz3(1.1633580248937700485,0.51825684289835904917,0.015000075099999999467);
  nxyz4(1.1633039901419599538,0.51827976384601803783,0);

  GEO::CUT::SideHandle * s3 = intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.0862518404756100754,0.54179921727284197086,-0.015000075099999999467);
  nxyz2(1.0863570485981299818,0.54164967699933597167,5.1690185600000003046e-19);
  nxyz3(1.1633039901419599538,0.51827976384601803783,0);
  nxyz4(1.1633580249045900601,0.5182568428569539476,-0.015000075099999999467);

  GEO::CUT::SideHandle * s4 = intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 1.1271155608629783718;
  tet4_xyze(1,0) = 0.52928799201809972885;
  tet4_xyze(2,0) = 0.0050000250499999997558;
  tet4_xyze(0,1) = 1.1633039901419599538;
  tet4_xyze(1,1) = 0.51827976384601803783;
  tet4_xyze(2,1) = 0;
  tet4_xyze(0,2) = 1.1530126441206562937;
  tet4_xyze(1,2) = 0.48640879576505308135;
  tet4_xyze(2,2) = 0.0050000250499999997558;
  tet4_xyze(0,3) = 1.127115560864781596;
  tet4_xyze(1,3) = 0.52928799201042275868;
  tet4_xyze(2,3) = -0.0050000250499999988885;

  nids.clear();
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  GEO::CUT::ElementHandle * e = intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  GEO::CUT::plain_side_set cut_sides;
  s1->CollectSides( cut_sides );
  s2->CollectSides( cut_sides );
  s3->CollectSides( cut_sides );
  s4->CollectSides( cut_sides );

  std::vector<std::vector<int> > tets( 1 );
  std::vector<int> & tet = tets.back();
  tet.reserve( 4 );
  std::vector<int> accept_tets( 1, true );

  const std::vector<GEO::CUT::Node*> & nodes = e->Nodes();
  std::vector<GEO::CUT::Point*> points;
  points.reserve( 4 );
  for ( std::vector<GEO::CUT::Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    GEO::CUT::Node * n = *i;
    tet.push_back( points.size() );
    points.push_back( n->point() );
  }

#if 0
  GEO::CUT::TetMeshIntersection tmi( tets, accept_tets, points, cut_sides );
  tmi.CutTest_Cut();
#else
  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
#endif
}

void test_tet4_quad4_alex10()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(1.1373961924174400551,0.43783278788480001209,0.015000075099999999467);
  nxyz2(1.1372558011327500438,0.43786377186006097961,1.3891340300000000002e-19);
  nxyz3(1.1633039901419599538,0.51827976384601803783,0);
  nxyz4(1.1633580248937700485,0.51825684289835904917,0.015000075099999999467);

  GEO::CUT::SideHandle * s1 = intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.1372558011327500438,0.43786377186006097961,1.3891340300000000002e-19);
  nxyz2(1.1373961924222799613,0.43783278790287999405,-0.015000075099999999467);
  nxyz3(1.1633580249045900601,0.5182568428569539476,-0.015000075099999999467);
  nxyz4(1.1633039901419599538,0.51827976384601803783,0);

  GEO::CUT::SideHandle * s2 = intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.0863570485981299818,0.54164967699933597167,5.1690185600000003046e-19);
  nxyz2(1.0862518404856800203,0.54179921727229496398,0.015000075099999999467);
  nxyz3(1.1633580248937700485,0.51825684289835904917,0.015000075099999999467);
  nxyz4(1.1633039901419599538,0.51827976384601803783,0);

  GEO::CUT::SideHandle * s3 = intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.0862518404756100754,0.54179921727284197086,-0.015000075099999999467);
  nxyz2(1.0863570485981299818,0.54164967699933597167,5.1690185600000003046e-19);
  nxyz3(1.1633039901419599538,0.51827976384601803783,0);
  nxyz4(1.1633580249045900601,0.5182568428569539476,-0.015000075099999999467);

  GEO::CUT::SideHandle * s4 = intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 1.1530126441231649537;
  tet4_xyze(1,0) = 0.4864087957581522681;
  tet4_xyze(2,0) = -0.0050000250499999988885;
  tet4_xyze(0,1) = 1.127115560864781596;
  tet4_xyze(1,1) = 0.52928799201042275868;
  tet4_xyze(2,1) = -0.0050000250499999988885;
  tet4_xyze(0,2) = 1.0909091200000000654;
  tet4_xyze(1,2) = 0.45454546800000000806;
  tet4_xyze(2,2) = -0.010000050099999999512;
  tet4_xyze(0,3) = 1.1530126441206562937;
  tet4_xyze(1,3) = 0.48640879576505308135;
  tet4_xyze(2,3) = 0.0050000250499999997558;

  nids.clear();
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  GEO::CUT::ElementHandle * e = intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  GEO::CUT::plain_side_set cut_sides;
  s1->CollectSides( cut_sides );
  s2->CollectSides( cut_sides );
  s3->CollectSides( cut_sides );
  s4->CollectSides( cut_sides );

  std::vector<std::vector<int> > tets( 1 );
  std::vector<int> & tet = tets.back();
  tet.reserve( 4 );
  std::vector<int> accept_tets( 1, true );

  const std::vector<GEO::CUT::Node*> & nodes = e->Nodes();
  std::vector<GEO::CUT::Point*> points;
  points.reserve( 4 );
  for ( std::vector<GEO::CUT::Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    GEO::CUT::Node * n = *i;
    tet.push_back( points.size() );
    points.push_back( n->point() );
  }

#if 0
  GEO::CUT::TetMeshIntersection tmi( tets, accept_tets, points, cut_sides );
  tmi.CutTest_Cut();
#else
  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
#endif
}

void test_tet4_quad4_alex11()
{
  GEO::CUT::MeshIntersection intersection;

  std::map<std::string, int> nodeids;

  {
    std::vector<int> nids( 4 );

    Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

    nxyz1(0.89075394384624373423,0.094851207190588718299,0.31907718246312877231);
    nxyz2(0.91085589872612870987,0.09462691177514795382,0.31915342860180828666);
    nxyz3(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
    nxyz4(0.89100776662223735158,0.11422875682038154121,0.31907009737708769137);

    intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );
  }

  {
    std::vector<int> nids( 3 );

    Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

    nxyz1(0.90127599046596340582,0.12378087657319215842,0.31909880503242982197);
    nxyz2(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
    nxyz3(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);

    intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.90127599046596340582,0.12378087657319215842,0.31909880503242982197);
    nxyz2(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);
    nxyz3(0.89130437399999995485,0.13360949336949076716,0.31906187923446299726);

    intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.90127599046596340582,0.12378087657319215842,0.31909880503242982197);
    nxyz2(0.89130437399999995485,0.13360949336949076716,0.31906187923446299726);
    nxyz3(0.89130437399999984383,0.11422492436052207598,0.31907111909968571828);

    intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.90127599046596340582,0.12378087657319215842,0.31909880503242982197);
    nxyz2(0.89130437399999984383,0.11422492436052207598,0.31907111909968571828);
    nxyz3(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);

    intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::tri3 );
  }

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 0.90127599046596340582;
  tet4_xyze(1,0) = 0.12378087657319215842;
  tet4_xyze(2,0) = 0.31909880503242982197;
  tet4_xyze(0,1) = 0.89130437399999995485;
  tet4_xyze(1,1) = 0.11363636700000000201;
  tet4_xyze(2,1) = 0.31907136403769043032;
  tet4_xyze(0,2) = 0.91110160767503201029;
  tet4_xyze(1,2) = 0.11363636700000001589;
  tet4_xyze(2,2) = 0.31913957190701630617;
  tet4_xyze(0,3) = 0.89130437399999995485;
  tet4_xyze(1,3) = 0.11363636700000000201;
  tet4_xyze(2,3) = 0.34999999399999998095;

  std::vector<int> nids;
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex12()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.89129457552795998687,0.13360963458984104202,0.31906184946907101274);
  nxyz2(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);
  nxyz3(0.91170007309449896393,0.15267896040857797946,0.31910503128749823087);
  nxyz4(0.89160886085045754079,0.15299394383896400273,0.31905289788465629464);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89100776662223735158,0.11422875682038154121,0.31907009737708769137);
  nxyz2(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
  nxyz3(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);
  nxyz4(0.89129457552795998687,0.13360963458984104202,0.31906184946907101274);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89075394384624373423,0.094851207190588718299,0.31907718246312877231);
  nxyz2(0.91085589872612870987,0.09462691177514795382,0.31915342860180828666);
  nxyz3(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
  nxyz4(0.89100776662223735158,0.11422875682038154121,0.31907009737708769137);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91085589872612870987,0.09462691177514795382,0.31915342860180828666);
  nxyz2(0.91092075060668153963,0.094621463232861316439,0.29905679848330324333);
  nxyz3(0.91116510843108489137,0.11396257312124993821,0.29904728417926035311);
  nxyz4(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);
  nxyz2(0.91144151174136411875,0.13331279245508514952,0.29903508638817211107);
  nxyz3(0.91174456614349386196,0.15267165352639497367,0.2990214845365067875);
  nxyz4(0.91170007309449896393,0.15267896040857797946,0.31910503128749823087);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.87121404901980936675,0.13387972670277956122,0.31899881128682672227);
  nxyz2(0.89129457552795998687,0.13360963458984104202,0.31906184946907101274);
  nxyz3(0.89160886085045754079,0.15299394383896400273,0.31905289788465629464);
  nxyz4(0.87152955219270600296,0.15328918881726930068,0.31899878565014916365);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
  nxyz2(0.91116510843108489137,0.11396257312124993821,0.29904728417926035311);
  nxyz3(0.91144151174136411875,0.13331279245508514952,0.29903508638817211107);
  nxyz4(0.91138930580985977326,0.13332001955212607891,0.31912289240735142171);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.91304349900000003615;
  hex8_xyze(1,0) = 0.13636364000000000818;
  hex8_xyze(2,0) = 0.30000001199999998258;
  hex8_xyze(0,1) = 0.89130437399999995485;
  hex8_xyze(1,1) = 0.13636364000000000818;
  hex8_xyze(2,1) = 0.30000001199999998258;
  hex8_xyze(0,2) = 0.89130437399999995485;
  hex8_xyze(1,2) = 0.11363636700000000201;
  hex8_xyze(2,2) = 0.30000001199999998258;
  hex8_xyze(0,3) = 0.91304349900000003615;
  hex8_xyze(1,3) = 0.11363636700000000201;
  hex8_xyze(2,3) = 0.30000001199999998258;
  hex8_xyze(0,4) = 0.91304349900000003615;
  hex8_xyze(1,4) = 0.13636364000000000818;
  hex8_xyze(2,4) = 0.34999999399999998095;
  hex8_xyze(0,5) = 0.89130437399999995485;
  hex8_xyze(1,5) = 0.13636364000000000818;
  hex8_xyze(2,5) = 0.34999999399999998095;
  hex8_xyze(0,6) = 0.89130437399999995485;
  hex8_xyze(1,6) = 0.11363636700000000201;
  hex8_xyze(2,6) = 0.34999999399999998095;
  hex8_xyze(0,7) = 0.91304349900000003615;
  hex8_xyze(1,7) = 0.11363636700000000201;
  hex8_xyze(2,7) = 0.34999999399999998095;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex13()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.90858273762257057982,0.14839508461400718886,0.31969337939658182268);
  nxyz2(0.92841430010507919324,0.1441334754323351508,0.32040076497947073131);
  nxyz3(0.93278815903190737124,0.16251264451130459365,0.3201473781817474884);
  nxyz4(0.91303749122979560582,0.16703848943081595069,0.31957292565401690387);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91771332751466661026,0.18566961672157589192,0.31945566376661693742);
  nxyz2(0.93737772492281612458,0.18093860723878399388,0.31990688123057475778);
  nxyz3(0.94214941835382437496,0.19941159083636084137,0.31968499159999225201);
  nxyz4(0.92255713529109029114,0.20430041973201396033,0.31934439946786041808);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.91303749122979560582,0.16703848943081595069,0.31957292565401690387);
  nxyz2(0.93278815903190737124,0.16251264451130459365,0.3201473781817474884);
  nxyz3(0.93737772492281612458,0.18093860723878399388,0.31990688123057475778);
  nxyz4(0.91771332751466661026,0.18566961672157589192,0.31945566376661693742);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.89336634885174737164,0.17127872118794218448,0.31897854401834946358);
  nxyz2(0.91303749122979560582,0.16703848943081595069,0.31957292565401690387);
  nxyz3(0.91771332751466661026,0.18566961672157589192,0.31945566376661693742);
  nxyz4(0.89809864231398062184,0.19013934730278916896,0.31898286201516340421);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.93278815903190737124,0.16251264451130459365,0.3201473781817474884);
  nxyz2(0.93324183024656148788,0.16230920279488567082,0.29987602045833194886);
  nxyz3(0.93773075623750146157,0.18076200977819448235,0.29968194775221990156);
  nxyz4(0.93737772492281612458,0.18093860723878399388,0.31990688123057475778);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.92841430010507919324,0.1441334754323351508,0.32040076497947073131);
  nxyz2(0.92899096297831118552,0.14390714335069226815,0.30007400109817816691);
  nxyz3(0.93324183024656148788,0.16230920279488567082,0.29987602045833194886);
  nxyz4(0.93278815903190737124,0.16251264451130459365,0.3201473781817474884);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.91304349900000003615;
  hex8_xyze(1,0) = 0.18181818699999999223;
  hex8_xyze(2,0) = 0.30000001199999998258;
  hex8_xyze(0,1) = 0.91304349900000003615;
  hex8_xyze(1,1) = 0.15909090600000000415;
  hex8_xyze(2,1) = 0.30000001199999998258;
  hex8_xyze(0,2) = 0.93478262400000000643;
  hex8_xyze(1,2) = 0.15909090600000000415;
  hex8_xyze(2,2) = 0.30000001199999998258;
  hex8_xyze(0,3) = 0.93478262400000000643;
  hex8_xyze(1,3) = 0.18181818699999999223;
  hex8_xyze(2,3) = 0.30000001199999998258;
  hex8_xyze(0,4) = 0.91304349900000003615;
  hex8_xyze(1,4) = 0.18181818699999999223;
  hex8_xyze(2,4) = 0.34999999399999998095;
  hex8_xyze(0,5) = 0.91304349900000003615;
  hex8_xyze(1,5) = 0.15909090600000000415;
  hex8_xyze(2,5) = 0.34999999399999998095;
  hex8_xyze(0,6) = 0.93478262400000000643;
  hex8_xyze(1,6) = 0.15909090600000000415;
  hex8_xyze(2,6) = 0.34999999399999998095;
  hex8_xyze(0,7) = 0.93478262400000000643;
  hex8_xyze(1,7) = 0.18181818699999999223;
  hex8_xyze(2,7) = 0.34999999399999998095;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex14()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1 (0.8690043, 0.4168291, 0.288951) ;
  nxyz2 (0.8955157, 0.4204433, 0.2888444) ;
  nxyz3 (0.8956636, 0.4204574, 0.3209518) ;
  nxyz4 (0.8691562, 0.4168354, 0.3210609) ;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8956636, 0.4204574, 0.3209518) ;
  nxyz2 (0.8913324, 0.45184, 0.3209449) ;
  nxyz3 (0.8648265, 0.4482161, 0.3210346) ;
  nxyz4 (0.8691562, 0.4168354, 0.3210609) ;

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8648265, 0.4482161, 0.3210346) ;
  nxyz2 (0.8647073, 0.448206, 0.2889356) ;
  nxyz3 (0.8690043, 0.4168291, 0.288951) ;
  nxyz4 (0.8691562, 0.4168354, 0.3210609) ;

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8691559, 0.3831646, 0.321061) ;
  nxyz2 (0.869004, 0.3831709, 0.2889511) ;
  nxyz3 (0.864707, 0.3517939, 0.2889357) ;
  nxyz4 (0.8648262, 0.3517838, 0.3210347) ;

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8691559, 0.3831646, 0.321061) ;
  nxyz2 (0.8956633, 0.3795426, 0.3209518) ;
  nxyz3 (0.8955154, 0.3795567, 0.2888444) ;
  nxyz4 (0.869004, 0.3831709, 0.2889511) ;

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8913322, 0.34816, 0.3209449) ;
  nxyz2 (0.8956633, 0.3795426, 0.3209518) ;
  nxyz3 (0.8691559, 0.3831646, 0.321061) ;
  nxyz4 (0.8648262, 0.3517838, 0.3210347) ;

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 8.846154e-01;
  hex8_xyze(1,0) = 4.190476e-01;
  hex8_xyze(2,0) = 3.333333e-01;
  hex8_xyze(0,1) = 8.846154e-01;
  hex8_xyze(1,1) = 3.809524e-01;
  hex8_xyze(2,1) = 3.333333e-01;
  hex8_xyze(0,2) = 8.461538e-01;
  hex8_xyze(1,2) = 3.809524e-01;
  hex8_xyze(2,2) = 3.333333e-01;
  hex8_xyze(0,3) = 8.461538e-01;
  hex8_xyze(1,3) = 4.190476e-01;
  hex8_xyze(2,3) = 3.333333e-01;
  hex8_xyze(0,4) = 8.846154e-01;
  hex8_xyze(1,4) = 4.190476e-01;
  hex8_xyze(2,4) = 3.055556e-01;
  hex8_xyze(0,5) = 8.846154e-01;
  hex8_xyze(1,5) = 3.809524e-01;
  hex8_xyze(2,5) = 3.055556e-01;
  hex8_xyze(0,6) = 8.461538e-01;
  hex8_xyze(1,6) = 3.809524e-01;
  hex8_xyze(2,6) = 3.055556e-01;
  hex8_xyze(0,7) = 8.461538e-01;
  hex8_xyze(1,7) = 4.190476e-01;
  hex8_xyze(2,7) = 3.055556e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex15()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1 (0.8687621, 0.3830149, 0.06015265) ;
  nxyz2 (0.8687252, 0.3830117, 0.04006952) ;
  nxyz3 (0.865972, 0.3620882, 0.0400713) ;
  nxyz4 (0.866005, 0.3620912, 0.0601553) ;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8660498, 0.4379045, 0.08023742) ;
  nxyz2 (0.866005, 0.4379088, 0.0601553) ;
  nxyz3 (0.8687621, 0.4169851, 0.06015264) ;
  nxyz4 (0.8688126, 0.4169804, 0.08023394) ;

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8688126, 0.4169804, 0.08023394) ;
  nxyz2 (0.8887074, 0.4195985, 0.08018836) ;
  nxyz3 (0.8887702, 0.4195995, 0.1002582) ;
  nxyz4 (0.8688754, 0.4169743, 0.1003129) ;

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8688126, 0.3830196, 0.08023395) ;
  nxyz2 (0.8887074, 0.3804015, 0.08018838) ;
  nxyz3 (0.888657, 0.3804024, 0.06011734) ;
  nxyz4 (0.8687621, 0.3830149, 0.06015265) ;

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8688754, 0.3830257, 0.1003129) ;
  nxyz2 (0.8688126, 0.3830196, 0.08023395) ;
  nxyz3 (0.8660498, 0.3620955, 0.08023743) ;
  nxyz4 (0.8661055, 0.3621014, 0.1003172) ;

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8661055, 0.4378986, 0.1003172) ;
  nxyz2 (0.8660498, 0.4379045, 0.08023742) ;
  nxyz3 (0.8688126, 0.4169804, 0.08023394) ;
  nxyz4 (0.8688754, 0.4169743, 0.1003129) ;

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8688126, 0.3830196, 0.08023395) ;
  nxyz2 (0.8687621, 0.3830149, 0.06015265) ;
  nxyz3 (0.866005, 0.3620912, 0.0601553) ;
  nxyz4 (0.8660498, 0.3620955, 0.08023743) ;

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8687621, 0.3830149, 0.06015265) ;
  nxyz2 (0.888657, 0.3804024, 0.06011734) ;
  nxyz3 (0.88862, 0.3804031, 0.04004542) ;
  nxyz4 (0.8687252, 0.3830117, 0.04006952) ;

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8688754, 0.3830257, 0.1003129) ;
  nxyz2 (0.8887702, 0.3804005, 0.1002582) ;
  nxyz3 (0.8887074, 0.3804015, 0.08018838) ;
  nxyz4 (0.8688126, 0.3830196, 0.08023395) ;

  intersection.AddCutSide( 9, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8687621, 0.4169851, 0.06015264) ;
  nxyz2 (0.888657, 0.4195976, 0.06011733) ;
  nxyz3 (0.8887074, 0.4195985, 0.08018836) ;
  nxyz4 (0.8688126, 0.4169804, 0.08023394) ;

  intersection.AddCutSide( 10, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.866005, 0.4379088, 0.0601553) ;
  nxyz2 (0.865972, 0.4379118, 0.04007129) ;
  nxyz3 (0.8687252, 0.4169883, 0.04006952) ;
  nxyz4 (0.8687621, 0.4169851, 0.06015264) ;

  intersection.AddCutSide( 11, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1 (0.8687252, 0.4169883, 0.04006952) ;
  nxyz2 (0.88862, 0.4195969, 0.04004542) ;
  nxyz3 (0.888657, 0.4195976, 0.06011733) ;
  nxyz4 (0.8687621, 0.4169851, 0.06015264) ;

  intersection.AddCutSide( 12, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 8.846154e-01;
  hex8_xyze(1,0) = 4.190476e-01;
  hex8_xyze(2,0) = 8.333333e-02;
  hex8_xyze(0,1) = 8.846154e-01;
  hex8_xyze(1,1) = 3.809524e-01;
  hex8_xyze(2,1) = 8.333333e-02;
  hex8_xyze(0,2) = 8.461538e-01;
  hex8_xyze(1,2) = 3.809524e-01;
  hex8_xyze(2,2) = 8.333333e-02;
  hex8_xyze(0,3) = 8.461538e-01;
  hex8_xyze(1,3) = 4.190476e-01;
  hex8_xyze(2,3) = 8.333333e-02;
  hex8_xyze(0,4) = 8.846154e-01;
  hex8_xyze(1,4) = 4.190476e-01;
  hex8_xyze(2,4) = 5.555556e-02;
  hex8_xyze(0,5) = 8.846154e-01;
  hex8_xyze(1,5) = 3.809524e-01;
  hex8_xyze(2,5) = 5.555556e-02;
  hex8_xyze(0,6) = 8.461538e-01;
  hex8_xyze(1,6) = 3.809524e-01;
  hex8_xyze(2,6) = 5.555556e-02;
  hex8_xyze(0,7) = 8.461538e-01;
  hex8_xyze(1,7) = 4.190476e-01;
  hex8_xyze(2,7) = 5.555556e-02;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_tet4_quad4_alex16()
{
  GEO::CUT::MeshIntersection intersection;

  std::map<std::string, int> nodeids;

//   {
//     std::vector<int> nids( 4 );

//     Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

//     nxyz1(0.89075394384624373423,0.094851207190588718299,0.31907718246312877231);
//     nxyz2(0.91085589872612870987,0.09462691177514795382,0.31915342860180828666);
//     nxyz3(0.91110590805399416237,0.11396906901062973938,0.31913932938821898411);
//     nxyz4(0.89100776662223735158,0.11422875682038154121,0.31907009737708769137);

//     intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );
//   }

  {
    std::vector<int> nids( 3 );

    Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

    nxyz1(0.899594,0.3799,0.016);
    nxyz2(0.9103,0.3799,-8.26403e-20);
    nxyz3(0.9103,0.3799,0.032);

    intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.376617,0.016);
    nxyz2(0.9103,0.373333,0);
    nxyz3(0.9103,0.373333,0.032);

    intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.376617,0.016);
    nxyz2(0.9103,0.373333,0.032);
    nxyz3(0.9103,0.3799,0.032);

    intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.376617,0.016);
    nxyz2(0.9103,0.3799,0.032);
    nxyz3(0.9103,0.3799,-8.26403e-20);

    intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.376617,0.016);
    nxyz2(0.9103,0.3799,-8.26403e-20);
    nxyz3(0.9103,0.373333,0);

    intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.899594,0.4201,0.016);
    nxyz2(0.9103,0.4201,0.032);
    nxyz3(0.9103,0.4201,-8.26403e-20);

    intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.899594,0.4201,0.016);
    nxyz2(0.9103,0.4201,-8.26403e-20);
    nxyz3(0.888889,0.4201,0);

    intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.423383,0.016);
    nxyz2(0.9103,0.4201,-8.26403e-20);
    nxyz3(0.9103,0.4201,0.032);

    intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.423383,0.016);
    nxyz2(0.9103,0.4201,0.032);
    nxyz3(0.9103,0.426667,0.032);

    intersection.AddCutSide( 9, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.423383,0.016);
    nxyz2(0.9103,0.426667,0.032);
    nxyz3(0.9103,0.426667,0);

    intersection.AddCutSide( 10, nids, quad4_xyze, DRT::Element::tri3 );

    nxyz1(0.9103,0.423383,0.016);
    nxyz2(0.9103,0.426667,0);
    nxyz3(0.9103,0.4201,-8.26403e-20);

    intersection.AddCutSide( 11, nids, quad4_xyze, DRT::Element::tri3 );
  }

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 0.9103;
  tet4_xyze(1,0) = 0.4201;
  tet4_xyze(2,0) = -8.26403e-20;
  tet4_xyze(0,1) = 0.9103;
  tet4_xyze(1,1) = 0.423383;
  tet4_xyze(2,1) = 0.016;
  tet4_xyze(0,2) = 0.9103;
  tet4_xyze(1,2) = 0.376617;
  tet4_xyze(2,2) = 0.016;
  tet4_xyze(0,3) = 0.944444;
  tet4_xyze(1,3) = 0.426667;
  tet4_xyze(2,3) = 0;

  std::vector<int> nids;
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex17()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.883533,0.3799,0.032);
  nxyz2(0.883533,0.3799,0.0641);
  nxyz3(0.9103,0.3799,0.0641);
  nxyz4(0.9103,0.3799,0.032);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.3799,-0.0001);
  nxyz2(0.883533,0.3799,0.032);
  nxyz3(0.9103,0.3799,0.032);
  nxyz4(0.9103,0.3799,-0.0001);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.348233,0.0641);
  nxyz2(0.9103,0.348233,0.032);
  nxyz3(0.9103,0.3799,0.032);
  nxyz4(0.9103,0.3799,0.0641);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.348233,0.032);
  nxyz2(0.9103,0.348233,-0.0001);
  nxyz3(0.9103,0.3799,-0.0001);
  nxyz4(0.9103,0.3799,0.032);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.4201,0.0641);
  nxyz2(0.883533,0.4201,0.032);
  nxyz3(0.9103,0.4201,0.032);
  nxyz4(0.9103,0.4201,0.0641);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.4201,0.032);
  nxyz2(0.883533,0.4201,-0.0001);
  nxyz3(0.9103,0.4201,-0.0001);
  nxyz4(0.9103,0.4201,0.032);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.4201,0.0641);
  nxyz2(0.9103,0.4201,0.032);
  nxyz3(0.9103,0.451767,0.032);
  nxyz4(0.9103,0.451767,0.0641);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.4201,0.032);
  nxyz2(0.9103,0.4201,-0.0001);
  nxyz3(0.9103,0.451767,-0.0001);
  nxyz4(0.9103,0.451767,0.032);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.944444;
  hex8_xyze(1,0) = 0.426667;
  hex8_xyze(2,0) = 0.0416667;
  hex8_xyze(0,1) = 0.944444;
  hex8_xyze(1,1) = 0.373333;
  hex8_xyze(2,1) = 0.0416667;
  hex8_xyze(0,2) = 0.888889;
  hex8_xyze(1,2) = 0.373333;
  hex8_xyze(2,2) = 0.0416667;
  hex8_xyze(0,3) = 0.888889;
  hex8_xyze(1,3) = 0.426667;
  hex8_xyze(2,3) = 0.0416667;
  hex8_xyze(0,4) = 0.944444;
  hex8_xyze(1,4) = 0.426667;
  hex8_xyze(2,4) = 0;
  hex8_xyze(0,5) = 0.944444;
  hex8_xyze(1,5) = 0.373333;
  hex8_xyze(2,5) = 0;
  hex8_xyze(0,6) = 0.888889;
  hex8_xyze(1,6) = 0.373333;
  hex8_xyze(2,6) = 0;
  hex8_xyze(0,7) = 0.888889;
  hex8_xyze(1,7) = 0.426667;
  hex8_xyze(2,7) = 0;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex18()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.302831e-01,3.444803e-01,1.923959e-01);
  nxyz2(9.301966e-01,3.445126e-01,1.603039e-01);
  nxyz3(9.326507e-01,3.760812e-01,1.603040e-01);
  nxyz4(9.327576e-01,3.760461e-01,1.923951e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.303734e-01,3.444429e-01,2.244904e-01);
  nxyz2(9.302831e-01,3.444803e-01,1.923959e-01);
  nxyz3(9.327576e-01,3.760461e-01,1.923951e-01);
  nxyz4(9.328715e-01,3.760068e-01,2.244875e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.059590e-01,3.781486e-01,1.603723e-01);
  nxyz2(9.060679e-01,3.781315e-01,1.924694e-01);
  nxyz3(9.327576e-01,3.760461e-01,1.923951e-01);
  nxyz4(9.326507e-01,3.760812e-01,1.603040e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.060679e-01,3.781315e-01,1.924694e-01);
  nxyz2(9.061837e-01,3.781127e-01,2.245668e-01);
  nxyz3(9.328715e-01,3.760068e-01,2.244875e-01);
  nxyz4(9.327576e-01,3.760461e-01,1.923951e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.253731e-01;
  hex8_xyze(1,0) = 3.703704e-01;
  hex8_xyze(2,0) = 1.764706e-01;
  hex8_xyze(0,1) = 9.328358e-01;
  hex8_xyze(1,1) = 3.703704e-01;
  hex8_xyze(2,1) = 1.764706e-01;
  hex8_xyze(0,2) = 9.328358e-01;
  hex8_xyze(1,2) = 3.777778e-01;
  hex8_xyze(2,2) = 1.764706e-01;
  hex8_xyze(0,3) = 9.253731e-01;
  hex8_xyze(1,3) = 3.777778e-01;
  hex8_xyze(2,3) = 1.764706e-01;
  hex8_xyze(0,4) = 9.253731e-01;
  hex8_xyze(1,4) = 3.703704e-01;
  hex8_xyze(2,4) = 2.058824e-01;
  hex8_xyze(0,5) = 9.328358e-01;
  hex8_xyze(1,5) = 3.703704e-01;
  hex8_xyze(2,5) = 2.058824e-01;
  hex8_xyze(0,6) = 9.328358e-01;
  hex8_xyze(1,6) = 3.777778e-01;
  hex8_xyze(2,6) = 2.058824e-01;
  hex8_xyze(0,7) = 9.253731e-01;
  hex8_xyze(1,7) = 3.777778e-01;
  hex8_xyze(2,7) = 2.058824e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex19()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.700810e-01,4.338294e-01,1.922623e-01);
  nxyz2(9.698407e-01,4.337088e-01,1.601800e-01);
  nxyz3(9.632220e-01,4.646628e-01,1.601806e-01);
  nxyz4(9.634112e-01,4.647694e-01,1.922662e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.441727e-01,4.282853e-01,2.245274e-01);
  nxyz2(9.439085e-01,4.282049e-01,1.924339e-01);
  nxyz3(9.700810e-01,4.338294e-01,1.922623e-01);
  nxyz4(9.703337e-01,4.339620e-01,2.243473e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.703337e-01,4.339620e-01,2.243473e-01);
  nxyz2(9.700810e-01,4.338294e-01,1.922623e-01);
  nxyz3(9.634112e-01,4.647694e-01,1.922662e-01);
  nxyz4(9.636037e-01,4.648893e-01,2.243579e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.439085e-01,4.282049e-01,1.924339e-01);
  nxyz2(9.436571e-01,4.281295e-01,1.603399e-01);
  nxyz3(9.698407e-01,4.337088e-01,1.601800e-01);
  nxyz4(9.700810e-01,4.338294e-01,1.922623e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.701493e-01;
  hex8_xyze(1,0) = 4.296296e-01;
  hex8_xyze(2,0) = 1.764706e-01;
  hex8_xyze(0,1) = 9.701493e-01;
  hex8_xyze(1,1) = 4.370370e-01;
  hex8_xyze(2,1) = 1.764706e-01;
  hex8_xyze(0,2) = 9.626866e-01;
  hex8_xyze(1,2) = 4.370370e-01;
  hex8_xyze(2,2) = 1.764706e-01;
  hex8_xyze(0,3) = 9.626866e-01;
  hex8_xyze(1,3) = 4.296296e-01;
  hex8_xyze(2,3) = 1.764706e-01;
  hex8_xyze(0,4) = 9.701493e-01;
  hex8_xyze(1,4) = 4.296296e-01;
  hex8_xyze(2,4) = 2.058824e-01;
  hex8_xyze(0,5) = 9.701493e-01;
  hex8_xyze(1,5) = 4.370370e-01;
  hex8_xyze(2,5) = 2.058824e-01;
  hex8_xyze(0,6) = 9.626866e-01;
  hex8_xyze(1,6) = 4.370370e-01;
  hex8_xyze(2,6) = 2.058824e-01;
  hex8_xyze(0,7) = 9.626866e-01;
  hex8_xyze(1,7) = 4.296296e-01;
  hex8_xyze(2,7) = 2.058824e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex20()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.856159e-01,3.268669e-01,1.601288e-01);
  nxyz2(9.853922e-01,3.270129e-01,1.280563e-01);
  nxyz3(9.947877e-01,3.572443e-01,1.280577e-01);
  nxyz4(9.950667e-01,3.570740e-01,1.601279e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.853922e-01,3.270129e-01,1.280563e-01);
  nxyz2(9.851978e-01,3.271313e-01,9.599957e-02);
  nxyz3(9.945510e-01,3.573831e-01,9.600163e-02);
  nxyz4(9.947877e-01,3.572443e-01,1.280577e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.761806e-01,2.967183e-01,1.601926e-01);
  nxyz2(9.760040e-01,2.968428e-01,1.281035e-01);
  nxyz3(9.853922e-01,3.270129e-01,1.280563e-01);
  nxyz4(9.856159e-01,3.268669e-01,1.601288e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.760040e-01,2.968428e-01,1.281035e-01);
  nxyz2(9.758444e-01,2.969431e-01,9.603294e-02);
  nxyz3(9.851978e-01,3.271313e-01,9.599957e-02);
  nxyz4(9.853922e-01,3.270129e-01,1.280563e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.850746e-01;
  hex8_xyze(1,0) = 3.259259e-01;
  hex8_xyze(2,0) = 1.176471e-01;
  hex8_xyze(0,1) = 1.000000e+00;
  hex8_xyze(1,1) = 3.259259e-01;
  hex8_xyze(2,1) = 1.176471e-01;
  hex8_xyze(0,2) = 1.000000e+00;
  hex8_xyze(1,2) = 3.407407e-01;
  hex8_xyze(2,2) = 1.176471e-01;
  hex8_xyze(0,3) = 9.850746e-01;
  hex8_xyze(1,3) = 3.407407e-01;
  hex8_xyze(2,3) = 1.176471e-01;
  hex8_xyze(0,4) = 9.850746e-01;
  hex8_xyze(1,4) = 3.259259e-01;
  hex8_xyze(2,4) = 1.470588e-01;
  hex8_xyze(0,5) = 1.000000e+00;
  hex8_xyze(1,5) = 3.259259e-01;
  hex8_xyze(2,5) = 1.470588e-01;
  hex8_xyze(0,6) = 1.000000e+00;
  hex8_xyze(1,6) = 3.407407e-01;
  hex8_xyze(2,6) = 1.470588e-01;
  hex8_xyze(0,7) = 9.850746e-01;
  hex8_xyze(1,7) = 3.407407e-01;
  hex8_xyze(2,7) = 1.470588e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex21()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.971,0.525,0.225);
  nxyz2(0.95,0.55,0.224);
  nxyz3(0.95,0.55,0.193);
  nxyz4(0.971,0.525,0.193);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.95,0.55,0.224);
  nxyz2(0.929,0.576,0.224);
  nxyz3(0.929,0.576,0.192);
  nxyz4(0.95,0.55,0.193);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.971,0.525,0.257);
  nxyz2(0.95,0.549,0.256);
  nxyz3(0.95,0.55,0.224);
  nxyz4(0.971,0.525,0.225);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.992,0.501,0.257);
  nxyz2(0.971,0.525,0.257);
  nxyz3(0.971,0.525,0.225);
  nxyz4(0.992,0.501,0.225);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.95,0.549,0.256);
  nxyz2(0.929,0.575,0.256);
  nxyz3(0.929,0.576,0.224);
  nxyz4(0.95,0.55,0.224);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.992,0.501,0.225);
  nxyz2(0.971,0.525,0.225);
  nxyz3(0.971,0.525,0.193);
  nxyz4(0.992,0.501,0.193);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.975;
  hex8_xyze(1,0) = 0.55;
  hex8_xyze(2,0) = 0.225;
  hex8_xyze(0,1) = 0.975;
  hex8_xyze(1,1) = 0.525;
  hex8_xyze(2,1) = 0.225;
  hex8_xyze(0,2) = 0.95;
  hex8_xyze(1,2) = 0.525;
  hex8_xyze(2,2) = 0.225;
  hex8_xyze(0,3) = 0.95;
  hex8_xyze(1,3) = 0.55;
  hex8_xyze(2,3) = 0.225;
  hex8_xyze(0,4) = 0.975;
  hex8_xyze(1,4) = 0.55;
  hex8_xyze(2,4) = 0.2;
  hex8_xyze(0,5) = 0.975;
  hex8_xyze(1,5) = 0.525;
  hex8_xyze(2,5) = 0.2;
  hex8_xyze(0,6) = 0.95;
  hex8_xyze(1,6) = 0.525;
  hex8_xyze(2,6) = 0.2;
  hex8_xyze(0,7) = 0.95;
  hex8_xyze(1,7) = 0.55;
  hex8_xyze(2,7) = 0.2;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex22()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.631016e-01,3.410436e-01,1.287286e-01);
  nxyz2(9.787952e-01,3.686359e-01,1.287249e-01);
  nxyz3(9.784376e-01,3.687304e-01,9.655058e-02);
  nxyz4(9.628123e-01,3.411170e-01,9.655388e-02);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.471839e-01,3.133966e-01,9.651074e-02);
  nxyz2(9.628123e-01,3.411170e-01,9.655388e-02);
  nxyz3(9.625893e-01,3.411801e-01,6.435218e-02);
  nxyz4(9.469989e-01,3.134496e-01,6.432435e-02);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.628123e-01,3.411170e-01,9.655388e-02);
  nxyz2(9.784376e-01,3.687304e-01,9.655058e-02);
  nxyz3(9.781685e-01,3.688079e-01,6.434987e-02);
  nxyz4(9.625893e-01,3.411801e-01,6.435218e-02);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.474189e-01,3.133388e-01,1.286687e-01);
  nxyz2(9.631016e-01,3.410436e-01,1.287286e-01);
  nxyz3(9.628123e-01,3.411170e-01,9.655388e-02);
  nxyz4(9.471839e-01,3.133966e-01,9.651074e-02);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.701493e-01;
  hex8_xyze(1,0) = 3.407407e-01;
  hex8_xyze(2,0) = 8.823529e-02;
  hex8_xyze(0,1) = 9.701493e-01;
  hex8_xyze(1,1) = 3.481481e-01;
  hex8_xyze(2,1) = 8.823529e-02;
  hex8_xyze(0,2) = 9.626866e-01;
  hex8_xyze(1,2) = 3.481481e-01;
  hex8_xyze(2,2) = 8.823529e-02;
  hex8_xyze(0,3) = 9.626866e-01;
  hex8_xyze(1,3) = 3.407407e-01;
  hex8_xyze(2,3) = 8.823529e-02;
  hex8_xyze(0,4) = 9.701493e-01;
  hex8_xyze(1,4) = 3.407407e-01;
  hex8_xyze(2,4) = 1.176471e-01;
  hex8_xyze(0,5) = 9.701493e-01;
  hex8_xyze(1,5) = 3.481481e-01;
  hex8_xyze(2,5) = 1.176471e-01;
  hex8_xyze(0,6) = 9.626866e-01;
  hex8_xyze(1,6) = 3.481481e-01;
  hex8_xyze(2,6) = 1.176471e-01;
  hex8_xyze(0,7) = 9.626866e-01;
  hex8_xyze(1,7) = 3.407407e-01;
  hex8_xyze(2,7) = 1.176471e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex23()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.883533,0.3799,0.032);
  nxyz2(0.883533,0.3799,0.0641);
  nxyz3(0.9103,0.3799,0.0641);
  nxyz4(0.9103,0.3799,0.032);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.3799,-0.0001);
  nxyz2(0.883533,0.3799,0.032);
  nxyz3(0.9103,0.3799,0.032);
  nxyz4(0.9103,0.3799,-0.0001);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.348233,0.0641);
  nxyz2(0.9103,0.348233,0.032);
  nxyz3(0.9103,0.3799,0.032);
  nxyz4(0.9103,0.3799,0.0641);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.348233,0.032);
  nxyz2(0.9103,0.348233,-0.0001);
  nxyz3(0.9103,0.3799,-0.0001);
  nxyz4(0.9103,0.3799,0.032);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.4201,0.0641);
  nxyz2(0.883533,0.4201,0.032);
  nxyz3(0.9103,0.4201,0.032);
  nxyz4(0.9103,0.4201,0.0641);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.4201,0.032);
  nxyz2(0.883533,0.4201,-0.0001);
  nxyz3(0.9103,0.4201,-0.0001);
  nxyz4(0.9103,0.4201,0.032);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.4201,0.0641);
  nxyz2(0.9103,0.4201,0.032);
  nxyz3(0.9103,0.451767,0.032);
  nxyz4(0.9103,0.451767,0.0641);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.4201,0.032);
  nxyz2(0.9103,0.4201,-0.0001);
  nxyz3(0.9103,0.451767,-0.0001);
  nxyz4(0.9103,0.451767,0.032);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.944444;
  hex8_xyze(1,0) = 0.426667;
  hex8_xyze(2,0) = 0.0416667;
  hex8_xyze(0,1) = 0.944444;
  hex8_xyze(1,1) = 0.373333;
  hex8_xyze(2,1) = 0.0416667;
  hex8_xyze(0,2) = 0.888889;
  hex8_xyze(1,2) = 0.373333;
  hex8_xyze(2,2) = 0.0416667;
  hex8_xyze(0,3) = 0.888889;
  hex8_xyze(1,3) = 0.426667;
  hex8_xyze(2,3) = 0.0416667;
  hex8_xyze(0,4) = 0.944444;
  hex8_xyze(1,4) = 0.426667;
  hex8_xyze(2,4) = 0;
  hex8_xyze(0,5) = 0.944444;
  hex8_xyze(1,5) = 0.373333;
  hex8_xyze(2,5) = 0;
  hex8_xyze(0,6) = 0.888889;
  hex8_xyze(1,6) = 0.373333;
  hex8_xyze(2,6) = 0;
  hex8_xyze(0,7) = 0.888889;
  hex8_xyze(1,7) = 0.426667;
  hex8_xyze(2,7) = 0;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex24()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.883533,0.3799,-0.0001);
  nxyz2(0.883533,0.3799,0.032);
  nxyz3(0.9103,0.3799,0.032);
  nxyz4(0.9103,0.3799,-0.0001);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.348233,0.032);
  nxyz2(0.9103,0.348233,-0.0001);
  nxyz3(0.9103,0.3799,-0.0001);
  nxyz4(0.9103,0.3799,0.032);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.883533,0.4201,0.032);
  nxyz2(0.883533,0.4201,-0.0001);
  nxyz3(0.9103,0.4201,-0.0001);
  nxyz4(0.9103,0.4201,0.032);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.9103,0.4201,0.032);
  nxyz2(0.9103,0.4201,-0.0001);
  nxyz3(0.9103,0.451767,-0.0001);
  nxyz4(0.9103,0.451767,0.032);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.93617;
  hex8_xyze(1,0) = 0.421053;
  hex8_xyze(2,0) = 0.03125;
  hex8_xyze(0,1) = 0.93617;
  hex8_xyze(1,1) = 0.378947;
  hex8_xyze(2,1) = 0.03125;
  hex8_xyze(0,2) = 0.893617;
  hex8_xyze(1,2) = 0.378947;
  hex8_xyze(2,2) = 0.03125;
  hex8_xyze(0,3) = 0.893617;
  hex8_xyze(1,3) = 0.421053;
  hex8_xyze(2,3) = 0.03125;
  hex8_xyze(0,4) = 0.93617;
  hex8_xyze(1,4) = 0.421053;
  hex8_xyze(2,4) = 0;
  hex8_xyze(0,5) = 0.93617;
  hex8_xyze(1,5) = 0.378947;
  hex8_xyze(2,5) = 0;
  hex8_xyze(0,6) = 0.893617;
  hex8_xyze(1,6) = 0.378947;
  hex8_xyze(2,6) = 0;
  hex8_xyze(0,7) = 0.893617;
  hex8_xyze(1,7) = 0.421053;
  hex8_xyze(2,7) = 0;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex25()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.692014e-01,4.167234e-01,2.888914e-01);
  nxyz3(8.956973e-01,4.204100e-01,2.888267e-01);
  nxyz4(8.957910e-01,4.204236e-01,3.209346e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.648857e-01,3.518945e-01,3.209665e-01);
  nxyz2(8.913832e-01,3.482061e-01,3.209214e-01);
  nxyz3(8.957910e-01,3.795764e-01,3.209346e-01);
  nxyz4(8.692972e-01,3.832692e-01,3.210027e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.648857e-01,4.481055e-01,3.209665e-01);
  nxyz3(8.648223e-01,4.480965e-01,2.888696e-01);
  nxyz4(8.692014e-01,4.167234e-01,2.888914e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692014e-01,3.832766e-01,2.888914e-01);
  nxyz2(8.692972e-01,3.832692e-01,3.210027e-01);
  nxyz3(8.957910e-01,3.795764e-01,3.209346e-01);
  nxyz4(8.956973e-01,3.795900e-01,2.888267e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.957910e-01,4.204236e-01,3.209346e-01);
  nxyz3(8.913832e-01,4.517939e-01,3.209214e-01);
  nxyz4(8.648857e-01,4.481055e-01,3.209665e-01);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.648857e-01,3.518945e-01,3.209665e-01);
  nxyz2(8.692972e-01,3.832692e-01,3.210027e-01);
  nxyz3(8.692014e-01,3.832766e-01,2.888914e-01);
  nxyz4(8.648223e-01,3.519035e-01,2.888696e-01);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 8.936170e-01;
  hex8_xyze(1,0) = 4.210526e-01;
  hex8_xyze(2,0) = 3.437500e-01;
  hex8_xyze(0,1) = 8.936170e-01;
  hex8_xyze(1,1) = 3.789474e-01;
  hex8_xyze(2,1) = 3.437500e-01;
  hex8_xyze(0,2) = 8.510638e-01;
  hex8_xyze(1,2) = 3.789474e-01;
  hex8_xyze(2,2) = 3.437500e-01;
  hex8_xyze(0,3) = 8.510638e-01;
  hex8_xyze(1,3) = 4.210526e-01;
  hex8_xyze(2,3) = 3.437500e-01;
  hex8_xyze(0,4) = 8.936170e-01;
  hex8_xyze(1,4) = 4.210526e-01;
  hex8_xyze(2,4) = 3.125000e-01;
  hex8_xyze(0,5) = 8.936170e-01;
  hex8_xyze(1,5) = 3.789474e-01;
  hex8_xyze(2,5) = 3.125000e-01;
  hex8_xyze(0,6) = 8.510638e-01;
  hex8_xyze(1,6) = 3.789474e-01;
  hex8_xyze(2,6) = 3.125000e-01;
  hex8_xyze(0,7) = 8.510638e-01;
  hex8_xyze(1,7) = 4.210526e-01;
  hex8_xyze(2,7) = 3.125000e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex26()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.692014e-01,4.167234e-01,2.888914e-01);
  nxyz3(8.956973e-01,4.204100e-01,2.888267e-01);
  nxyz4(8.957910e-01,4.204236e-01,3.209346e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.648857e-01,3.518945e-01,3.209665e-01);
  nxyz2(8.913832e-01,3.482061e-01,3.209214e-01);
  nxyz3(8.957910e-01,3.795764e-01,3.209346e-01);
  nxyz4(8.692972e-01,3.832692e-01,3.210027e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.648857e-01,4.481055e-01,3.209665e-01);
  nxyz3(8.648223e-01,4.480965e-01,2.888696e-01);
  nxyz4(8.692014e-01,4.167234e-01,2.888914e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692014e-01,3.832766e-01,2.888914e-01);
  nxyz2(8.692972e-01,3.832692e-01,3.210027e-01);
  nxyz3(8.957910e-01,3.795764e-01,3.209346e-01);
  nxyz4(8.956973e-01,3.795900e-01,2.888267e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.692972e-01,4.167308e-01,3.210027e-01);
  nxyz2(8.957910e-01,4.204236e-01,3.209346e-01);
  nxyz3(8.913832e-01,4.517939e-01,3.209214e-01);
  nxyz4(8.648857e-01,4.481055e-01,3.209665e-01);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.648857e-01,3.518945e-01,3.209665e-01);
  nxyz2(8.692972e-01,3.832692e-01,3.210027e-01);
  nxyz3(8.692014e-01,3.832766e-01,2.888914e-01);
  nxyz4(8.648223e-01,3.519035e-01,2.888696e-01);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 8.936170e-01;
  hex8_xyze(1,0) = 4.210526e-01;
  hex8_xyze(2,0) = 3.437500e-01;
  hex8_xyze(0,1) = 8.936170e-01;
  hex8_xyze(1,1) = 3.789474e-01;
  hex8_xyze(2,1) = 3.437500e-01;
  hex8_xyze(0,2) = 8.510638e-01;
  hex8_xyze(1,2) = 3.789474e-01;
  hex8_xyze(2,2) = 3.437500e-01;
  hex8_xyze(0,3) = 8.510638e-01;
  hex8_xyze(1,3) = 4.210526e-01;
  hex8_xyze(2,3) = 3.437500e-01;
  hex8_xyze(0,4) = 8.936170e-01;
  hex8_xyze(1,4) = 4.210526e-01;
  hex8_xyze(2,4) = 3.125000e-01;
  hex8_xyze(0,5) = 8.936170e-01;
  hex8_xyze(1,5) = 3.789474e-01;
  hex8_xyze(2,5) = 3.125000e-01;
  hex8_xyze(0,6) = 8.510638e-01;
  hex8_xyze(1,6) = 3.789474e-01;
  hex8_xyze(2,6) = 3.125000e-01;
  hex8_xyze(0,7) = 8.510638e-01;
  hex8_xyze(1,7) = 4.210526e-01;
  hex8_xyze(2,7) = 3.125000e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex27()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.759851e-01,2.968514e-01,1.281036e-01);
  nxyz2(9.758256e-01,2.969516e-01,9.603311e-02);
  nxyz3(9.851760e-01,3.271408e-01,9.599976e-02);
  nxyz4(9.853703e-01,3.270224e-01,1.280565e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.761614e-01,2.967270e-01,1.601928e-01);
  nxyz2(9.759851e-01,2.968514e-01,1.281036e-01);
  nxyz3(9.853703e-01,3.270224e-01,1.280565e-01);
  nxyz4(9.855938e-01,3.268765e-01,1.601290e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.853703e-01,3.270224e-01,1.280565e-01);
  nxyz2(9.851760e-01,3.271408e-01,9.599976e-02);
  nxyz3(9.945263e-01,3.573934e-01,9.600185e-02);
  nxyz4(9.947628e-01,3.572547e-01,1.280579e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.855938e-01,3.268765e-01,1.601290e-01);
  nxyz2(9.853703e-01,3.270224e-01,1.280565e-01);
  nxyz3(9.947628e-01,3.572547e-01,1.280579e-01);
  nxyz4(9.950415e-01,3.570846e-01,1.601281e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.850746e-01;
  hex8_xyze(1,0) = 3.259259e-01;
  hex8_xyze(2,0) = 1.176471e-01;
  hex8_xyze(0,1) = 1.000000e+00;
  hex8_xyze(1,1) = 3.259259e-01;
  hex8_xyze(2,1) = 1.176471e-01;
  hex8_xyze(0,2) = 1.000000e+00;
  hex8_xyze(1,2) = 3.407407e-01;
  hex8_xyze(2,2) = 1.176471e-01;
  hex8_xyze(0,3) = 9.850746e-01;
  hex8_xyze(1,3) = 3.407407e-01;
  hex8_xyze(2,3) = 1.176471e-01;
  hex8_xyze(0,4) = 9.850746e-01;
  hex8_xyze(1,4) = 3.259259e-01;
  hex8_xyze(2,4) = 1.470588e-01;
  hex8_xyze(0,5) = 1.000000e+00;
  hex8_xyze(1,5) = 3.259259e-01;
  hex8_xyze(2,5) = 1.470588e-01;
  hex8_xyze(0,6) = 1.000000e+00;
  hex8_xyze(1,6) = 3.407407e-01;
  hex8_xyze(2,6) = 1.470588e-01;
  hex8_xyze(0,7) = 9.850746e-01;
  hex8_xyze(1,7) = 3.407407e-01;
  hex8_xyze(2,7) = 1.470588e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex28()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.865,0.0973,0.32);
  nxyz2(0.891,0.0928,0.322);
  nxyz3(0.898,0.123,0.322);
  nxyz4(0.872,0.129,0.32);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.898,0.123,0.322);
  nxyz2(0.925,0.117,0.323);
  nxyz3(0.933,0.146,0.323);
  nxyz4(0.907,0.153,0.322);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.925,0.117,0.323);
  nxyz2(0.926,0.116,0.291);
  nxyz3(0.934,0.145,0.29);
  nxyz4(0.933,0.146,0.323);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.891,0.0928,0.322);
  nxyz2(0.918,0.0877,0.324);
  nxyz3(0.925,0.117,0.323);
  nxyz4(0.898,0.123,0.322);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.918,0.0877,0.324);
  nxyz2(0.92,0.0873,0.291);
  nxyz3(0.926,0.116,0.291);
  nxyz4(0.925,0.117,0.323);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.886,0.0619,0.322);
  nxyz2(0.913,0.0586,0.324);
  nxyz3(0.918,0.0877,0.324);
  nxyz4(0.891,0.0928,0.322);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.925;
  hex8_xyze(1,0) = 0.0889;
  hex8_xyze(2,0) = 0.294;
  hex8_xyze(0,1) = 0.925;
  hex8_xyze(1,1) = 0.119;
  hex8_xyze(2,1) = 0.294;
  hex8_xyze(0,2) = 0.896;
  hex8_xyze(1,2) = 0.119;
  hex8_xyze(2,2) = 0.294;
  hex8_xyze(0,3) = 0.896;
  hex8_xyze(1,3) = 0.0889;
  hex8_xyze(2,3) = 0.294;
  hex8_xyze(0,4) = 0.925;
  hex8_xyze(1,4) = 0.0889;
  hex8_xyze(2,4) = 0.324;
  hex8_xyze(0,5) = 0.925;
  hex8_xyze(1,5) = 0.119;
  hex8_xyze(2,5) = 0.324;
  hex8_xyze(0,6) = 0.896;
  hex8_xyze(1,6) = 0.119;
  hex8_xyze(2,6) = 0.324;
  hex8_xyze(0,7) = 0.896;
  hex8_xyze(1,7) = 0.0889;
  hex8_xyze(2,7) = 0.324;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex29()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.632214e-01,4.646626e-01,1.601807e-01);
  nxyz2(9.630460e-01,4.645722e-01,1.281052e-01);
  nxyz3(9.564701e-01,4.954965e-01,1.281392e-01);
  nxyz4(9.566099e-01,4.955746e-01,1.602263e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.630460e-01,4.645722e-01,1.281052e-01);
  nxyz2(9.628950e-01,4.644997e-01,9.604120e-02);
  nxyz3(9.563455e-01,4.954348e-01,9.606543e-02);
  nxyz4(9.564701e-01,4.954965e-01,1.281392e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.698401e-01,4.337086e-01,1.601801e-01);
  nxyz2(9.696236e-01,4.336045e-01,1.281060e-01);
  nxyz3(9.630460e-01,4.645722e-01,1.281052e-01);
  nxyz4(9.632214e-01,4.646626e-01,1.601807e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.696236e-01,4.336045e-01,1.281060e-01);
  nxyz2(9.694408e-01,4.335201e-01,9.604240e-02);
  nxyz3(9.628950e-01,4.644997e-01,9.604120e-02);
  nxyz4(9.630460e-01,4.645722e-01,1.281052e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.701493e-01;
  hex8_xyze(1,0) = 4.592593e-01;
  hex8_xyze(2,0) = 1.176471e-01;
  hex8_xyze(0,1) = 9.701493e-01;
  hex8_xyze(1,1) = 4.666667e-01;
  hex8_xyze(2,1) = 1.176471e-01;
  hex8_xyze(0,2) = 9.626866e-01;
  hex8_xyze(1,2) = 4.666667e-01;
  hex8_xyze(2,2) = 1.176471e-01;
  hex8_xyze(0,3) = 9.626866e-01;
  hex8_xyze(1,3) = 4.592593e-01;
  hex8_xyze(2,3) = 1.176471e-01;
  hex8_xyze(0,4) = 9.701493e-01;
  hex8_xyze(1,4) = 4.592593e-01;
  hex8_xyze(2,4) = 1.470588e-01;
  hex8_xyze(0,5) = 9.701493e-01;
  hex8_xyze(1,5) = 4.666667e-01;
  hex8_xyze(2,5) = 1.470588e-01;
  hex8_xyze(0,6) = 9.626866e-01;
  hex8_xyze(1,6) = 4.666667e-01;
  hex8_xyze(2,6) = 1.470588e-01;
  hex8_xyze(0,7) = 9.626866e-01;
  hex8_xyze(1,7) = 4.592593e-01;
  hex8_xyze(2,7) = 1.470588e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex30()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.697729e-01,4.336404e-01,1.601536e-01);
  nxyz2(9.695358e-01,4.335307e-01,1.280851e-01);
  nxyz3(9.629934e-01,4.645134e-01,1.280847e-01);
  nxyz4(9.631879e-01,4.646092e-01,1.601547e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.631879e-01,4.646092e-01,1.601547e-01);
  nxyz2(9.629934e-01,4.645134e-01,1.280847e-01);
  nxyz3(9.564486e-01,4.954533e-01,1.281203e-01);
  nxyz4(9.566062e-01,4.955365e-01,1.602023e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.634008e-01,4.647223e-01,1.922347e-01);
  nxyz2(9.631879e-01,4.646092e-01,1.601547e-01);
  nxyz3(9.566062e-01,4.955365e-01,1.602023e-01);
  nxyz4(9.567716e-01,4.956366e-01,1.922968e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.700387e-01,4.337679e-01,1.922301e-01);
  nxyz2(9.697729e-01,4.336404e-01,1.601536e-01);
  nxyz3(9.631879e-01,4.646092e-01,1.601547e-01);
  nxyz4(9.634008e-01,4.647223e-01,1.922347e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.701493e-01;
  hex8_xyze(1,0) = 4.592593e-01;
  hex8_xyze(2,0) = 1.470588e-01;
  hex8_xyze(0,1) = 9.701493e-01;
  hex8_xyze(1,1) = 4.666667e-01;
  hex8_xyze(2,1) = 1.470588e-01;
  hex8_xyze(0,2) = 9.626866e-01;
  hex8_xyze(1,2) = 4.666667e-01;
  hex8_xyze(2,2) = 1.470588e-01;
  hex8_xyze(0,3) = 9.626866e-01;
  hex8_xyze(1,3) = 4.592593e-01;
  hex8_xyze(2,3) = 1.470588e-01;
  hex8_xyze(0,4) = 9.701493e-01;
  hex8_xyze(1,4) = 4.592593e-01;
  hex8_xyze(2,4) = 1.764706e-01;
  hex8_xyze(0,5) = 9.701493e-01;
  hex8_xyze(1,5) = 4.666667e-01;
  hex8_xyze(2,5) = 1.764706e-01;
  hex8_xyze(0,6) = 9.626866e-01;
  hex8_xyze(1,6) = 4.666667e-01;
  hex8_xyze(2,6) = 1.764706e-01;
  hex8_xyze(0,7) = 9.626866e-01;
  hex8_xyze(1,7) = 4.592593e-01;
  hex8_xyze(2,7) = 1.764706e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true,INPAR::CUT::VCellGaussPts_Tessellation);
  intersection.Status();
}

void test_hex8_quad4_alex31()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.007e-01,3.514e-01,2.249e-01);
  nxyz2(9.093e-01,3.819e-01,2.249e-01);
  nxyz3(9.089e-01,3.819e-01,1.928e-01);
  nxyz4(9.004e-01,3.514e-01,1.928e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.004e-01,3.514e-01,1.928e-01);
  nxyz2(9.089e-01,3.819e-01,1.928e-01);
  nxyz3(9.086e-01,3.819e-01,1.607e-01);
  nxyz4(9.002e-01,3.514e-01,1.607e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.921e-01,3.208e-01,2.248e-01);
  nxyz2(9.007e-01,3.514e-01,2.249e-01);
  nxyz3(9.004e-01,3.514e-01,1.928e-01);
  nxyz4(8.919e-01,3.208e-01,1.927e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.919e-01,3.208e-01,1.927e-01);
  nxyz2(9.004e-01,3.514e-01,1.928e-01);
  nxyz3(9.002e-01,3.514e-01,1.607e-01);
  nxyz4(8.917e-01,3.208e-01,1.606e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.250e-01;
  hex8_xyze(1,0) = 3.750e-01;
  hex8_xyze(2,0) = 2.000e-01;
  hex8_xyze(0,1) = 9.250e-01;
  hex8_xyze(1,1) = 3.500e-01;
  hex8_xyze(2,1) = 2.000e-01;
  hex8_xyze(0,2) = 9.000e-01;
  hex8_xyze(1,2) = 3.500e-01;
  hex8_xyze(2,2) = 2.000e-01;
  hex8_xyze(0,3) = 9.000e-01;
  hex8_xyze(1,3) = 3.750e-01;
  hex8_xyze(2,3) = 2.000e-01;
  hex8_xyze(0,4) = 9.250e-01;
  hex8_xyze(1,4) = 3.750e-01;
  hex8_xyze(2,4) = 1.750e-01;
  hex8_xyze(0,5) = 9.250e-01;
  hex8_xyze(1,5) = 3.500e-01;
  hex8_xyze(2,5) = 1.750e-01;
  hex8_xyze(0,6) = 9.000e-01;
  hex8_xyze(1,6) = 3.500e-01;
  hex8_xyze(2,6) = 1.750e-01;
  hex8_xyze(0,7) = 9.000e-01;
  hex8_xyze(1,7) = 3.750e-01;
  hex8_xyze(2,7) = 1.750e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex32()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.006951e-01,3.514136e-01,2.249071e-01);
  nxyz2(9.092586e-01,3.819322e-01,2.249058e-01);
  nxyz3(9.089236e-01,3.819199e-01,1.927974e-01);
  nxyz4(9.004321e-01,3.513832e-01,1.928019e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.921272e-01,3.208394e-01,2.248295e-01);
  nxyz2(9.006951e-01,3.514136e-01,2.249071e-01);
  nxyz3(9.004321e-01,3.513832e-01,1.928019e-01);
  nxyz4(8.919317e-01,3.207909e-01,1.927414e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(8.919317e-01,3.207909e-01,1.927414e-01);
  nxyz2(9.004321e-01,3.513832e-01,1.928019e-01);
  nxyz3(9.001781e-01,3.513680e-01,1.606864e-01);
  nxyz4(8.917339e-01,3.207654e-01,1.606396e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.004321e-01,3.513832e-01,1.928019e-01);
  nxyz2(9.089236e-01,3.819199e-01,1.927974e-01);
  nxyz3(9.086066e-01,3.819159e-01,1.606810e-01);
  nxyz4(9.001781e-01,3.513680e-01,1.606864e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.250000e-01;
  hex8_xyze(1,0) = 3.750000e-01;
  hex8_xyze(2,0) = 2.000000e-01;
  hex8_xyze(0,1) = 9.250000e-01;
  hex8_xyze(1,1) = 3.500000e-01;
  hex8_xyze(2,1) = 2.000000e-01;
  hex8_xyze(0,2) = 9.000000e-01;
  hex8_xyze(1,2) = 3.500000e-01;
  hex8_xyze(2,2) = 2.000000e-01;
  hex8_xyze(0,3) = 9.000000e-01;
  hex8_xyze(1,3) = 3.750000e-01;
  hex8_xyze(2,3) = 2.000000e-01;
  hex8_xyze(0,4) = 9.250000e-01;
  hex8_xyze(1,4) = 3.750000e-01;
  hex8_xyze(2,4) = 1.750000e-01;
  hex8_xyze(0,5) = 9.250000e-01;
  hex8_xyze(1,5) = 3.500000e-01;
  hex8_xyze(2,5) = 1.750000e-01;
  hex8_xyze(0,6) = 9.000000e-01;
  hex8_xyze(1,6) = 3.500000e-01;
  hex8_xyze(2,6) = 1.750000e-01;
  hex8_xyze(0,7) = 9.000000e-01;
  hex8_xyze(1,7) = 3.750000e-01;
  hex8_xyze(2,7) = 1.750000e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex33()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.498530e-01,4.279690e-01,2.565073e-01);
  nxyz2(9.496624e-01,4.278901e-01,2.244228e-01);
  nxyz3(9.452894e-01,4.592464e-01,2.244288e-01);
  nxyz4(9.454342e-01,4.593227e-01,2.565221e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.452894e-01,4.592464e-01,2.244288e-01);
  nxyz2(9.451478e-01,4.591754e-01,1.923355e-01);
  nxyz3(9.408188e-01,4.905113e-01,1.923740e-01);
  nxyz4(9.409217e-01,4.905776e-01,2.244782e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.496624e-01,4.278901e-01,2.244228e-01);
  nxyz2(9.494806e-01,4.278134e-01,1.923335e-01);
  nxyz3(9.451478e-01,4.591754e-01,1.923355e-01);
  nxyz4(9.452894e-01,4.592464e-01,2.244288e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.454342e-01,4.593227e-01,2.565221e-01);
  nxyz2(9.452894e-01,4.592464e-01,2.244288e-01);
  nxyz3(9.409217e-01,4.905776e-01,2.244782e-01);
  nxyz4(9.410197e-01,4.906515e-01,2.565843e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.552239e-01;
  hex8_xyze(1,0) = 4.444444e-01;
  hex8_xyze(2,0) = 2.058824e-01;
  hex8_xyze(0,1) = 9.552239e-01;
  hex8_xyze(1,1) = 4.592593e-01;
  hex8_xyze(2,1) = 2.058824e-01;
  hex8_xyze(0,2) = 9.402985e-01;
  hex8_xyze(1,2) = 4.592593e-01;
  hex8_xyze(2,2) = 2.058824e-01;
  hex8_xyze(0,3) = 9.402985e-01;
  hex8_xyze(1,3) = 4.444444e-01;
  hex8_xyze(2,3) = 2.058824e-01;
  hex8_xyze(0,4) = 9.552239e-01;
  hex8_xyze(1,4) = 4.444444e-01;
  hex8_xyze(2,4) = 2.352941e-01;
  hex8_xyze(0,5) = 9.552239e-01;
  hex8_xyze(1,5) = 4.592593e-01;
  hex8_xyze(2,5) = 2.352941e-01;
  hex8_xyze(0,6) = 9.402985e-01;
  hex8_xyze(1,6) = 4.592593e-01;
  hex8_xyze(2,6) = 2.352941e-01;
  hex8_xyze(0,7) = 9.402985e-01;
  hex8_xyze(1,7) = 4.444444e-01;
  hex8_xyze(2,7) = 2.352941e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true ,INPAR::CUT::VCellGaussPts_DirectDivergence, INPAR::CUT::BCellGaussPts_Tessellation);
  intersection.Status();
}

void test_hex8_quad4_alex34()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.654395e-01,3.724018e-01,1.929025e-01);
  nxyz2(9.658465e-01,3.723620e-01,2.249979e-01);
  nxyz3(9.896019e-01,3.601832e-01,2.247724e-01);
  nxyz4(9.892369e-01,3.603100e-01,1.926747e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.513290e-01,3.441290e-01,2.249951e-01);
  nxyz2(9.658465e-01,3.723620e-01,2.249979e-01);
  nxyz3(9.654395e-01,3.724018e-01,1.929025e-01);
  nxyz4(9.510218e-01,3.441215e-01,1.929042e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.658465e-01,3.723620e-01,2.249979e-01);
  nxyz2(9.662598e-01,3.723388e-01,2.571007e-01);
  nxyz3(9.899602e-01,3.600738e-01,2.568700e-01);
  nxyz4(9.896019e-01,3.601832e-01,2.247724e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.516350e-01,3.441686e-01,2.570851e-01);
  nxyz2(9.662598e-01,3.723388e-01,2.571007e-01);
  nxyz3(9.658465e-01,3.723620e-01,2.249979e-01);
  nxyz4(9.513290e-01,3.441290e-01,2.249951e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.750000e-01;
  hex8_xyze(1,0) = 3.750000e-01;
  hex8_xyze(2,0) = 2.250000e-01;
  hex8_xyze(0,1) = 9.750000e-01;
  hex8_xyze(1,1) = 3.500000e-01;
  hex8_xyze(2,1) = 2.250000e-01;
  hex8_xyze(0,2) = 9.500000e-01;
  hex8_xyze(1,2) = 3.500000e-01;
  hex8_xyze(2,2) = 2.250000e-01;
  hex8_xyze(0,3) = 9.500000e-01;
  hex8_xyze(1,3) = 3.750000e-01;
  hex8_xyze(2,3) = 2.250000e-01;
  hex8_xyze(0,4) = 9.750000e-01;
  hex8_xyze(1,4) = 3.750000e-01;
  hex8_xyze(2,4) = 2.000000e-01;
  hex8_xyze(0,5) = 9.750000e-01;
  hex8_xyze(1,5) = 3.500000e-01;
  hex8_xyze(2,5) = 2.000000e-01;
  hex8_xyze(0,6) = 9.500000e-01;
  hex8_xyze(1,6) = 3.500000e-01;
  hex8_xyze(2,6) = 2.000000e-01;
  hex8_xyze(0,7) = 9.500000e-01;
  hex8_xyze(1,7) = 3.750000e-01;
  hex8_xyze(2,7) = 2.000000e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex35()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.557365e-01,3.428832e-01,1.286608e-01);
  nxyz2(9.705318e-01,3.709260e-01,1.286612e-01);
  nxyz3(9.702512e-01,3.709883e-01,9.650093e-02);
  nxyz4(9.555130e-01,3.429299e-01,9.650103e-02);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.409015e-01,3.147291e-01,1.286048e-01);
  nxyz2(9.557365e-01,3.428832e-01,1.286608e-01);
  nxyz3(9.555130e-01,3.429299e-01,9.650103e-02);
  nxyz4(9.407246e-01,3.147650e-01,9.646085e-02);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.555130e-01,3.429299e-01,9.650103e-02);
  nxyz2(9.702512e-01,3.709883e-01,9.650093e-02);
  nxyz3(9.700342e-01,3.710408e-01,6.431592e-02);
  nxyz4(9.553374e-01,3.429723e-01,6.431625e-02);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.407246e-01,3.147650e-01,9.646085e-02);
  nxyz2(9.555130e-01,3.429299e-01,9.650103e-02);
  nxyz3(9.553374e-01,3.429723e-01,6.431625e-02);
  nxyz4(9.405818e-01,3.148011e-01,6.429014e-02);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.615385e-01;
  hex8_xyze(1,0) = 3.809524e-01;
  hex8_xyze(2,0) = 1.111111e-01;
  hex8_xyze(0,1) = 9.615385e-01;
  hex8_xyze(1,1) = 3.428571e-01;
  hex8_xyze(2,1) = 1.111111e-01;
  hex8_xyze(0,2) = 9.230769e-01;
  hex8_xyze(1,2) = 3.428571e-01;
  hex8_xyze(2,2) = 1.111111e-01;
  hex8_xyze(0,3) = 9.230769e-01;
  hex8_xyze(1,3) = 3.809524e-01;
  hex8_xyze(2,3) = 1.111111e-01;
  hex8_xyze(0,4) = 9.615385e-01;
  hex8_xyze(1,4) = 3.809524e-01;
  hex8_xyze(2,4) = 8.333333e-02;
  hex8_xyze(0,5) = 9.615385e-01;
  hex8_xyze(1,5) = 3.428571e-01;
  hex8_xyze(2,5) = 8.333333e-02;
  hex8_xyze(0,6) = 9.230769e-01;
  hex8_xyze(1,6) = 3.428571e-01;
  hex8_xyze(2,6) = 8.333333e-02;
  hex8_xyze(0,7) = 9.230769e-01;
  hex8_xyze(1,7) = 3.809524e-01;
  hex8_xyze(2,7) = 8.333333e-02;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex36()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(1.013108e+00,4.491975e-01,1.927528e-01);
  nxyz2(1.012861e+00,4.490732e-01,1.606404e-01);
  nxyz3(1.035538e+00,4.631656e-01,1.605541e-01);
  nxyz4(1.035744e+00,4.633647e-01,1.926876e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.012861e+00,4.490732e-01,1.606404e-01);
  nxyz2(1.012618e+00,4.489440e-01,1.285139e-01);
  nxyz3(1.035330e+00,4.629754e-01,1.284178e-01);
  nxyz4(1.035538e+00,4.631656e-01,1.605541e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.058235e+00,4.772966e-01,1.604059e-01);
  nxyz2(1.058066e+00,4.770501e-01,1.282692e-01);
  nxyz3(1.041395e+00,5.038824e-01,1.282505e-01);
  nxyz4(1.041505e+00,5.040821e-01,1.603877e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.041505e+00,5.040821e-01,1.603877e-01);
  nxyz2(1.041395e+00,5.038824e-01,1.282505e-01);
  nxyz3(1.024790e+00,5.306275e-01,1.282941e-01);
  nxyz4(1.024853e+00,5.307872e-01,1.604537e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.035744e+00,4.633647e-01,1.926876e-01);
  nxyz2(1.035538e+00,4.631656e-01,1.605541e-01);
  nxyz3(1.058235e+00,4.772966e-01,1.604059e-01);
  nxyz4(1.058396e+00,4.775672e-01,1.925550e-01);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.035538e+00,4.631656e-01,1.605541e-01);
  nxyz2(1.035330e+00,4.629754e-01,1.284178e-01);
  nxyz3(1.058066e+00,4.770501e-01,1.282692e-01);
  nxyz4(1.058235e+00,4.772966e-01,1.604059e-01);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.058396e+00,4.775672e-01,1.925550e-01);
  nxyz2(1.058235e+00,4.772966e-01,1.604059e-01);
  nxyz3(1.041505e+00,5.040821e-01,1.603877e-01);
  nxyz4(1.041586e+00,5.042994e-01,1.925434e-01);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(1.041586e+00,5.042994e-01,1.925434e-01);
  nxyz2(1.041505e+00,5.040821e-01,1.603877e-01);
  nxyz3(1.024853e+00,5.307872e-01,1.604537e-01);
  nxyz4(1.024866e+00,5.309578e-01,1.926398e-01);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1.063830e+00;
  hex8_xyze(1,0) = 5.052632e-01;
  hex8_xyze(2,0) = 1.875000e-01;
  hex8_xyze(0,1) = 1.063830e+00;
  hex8_xyze(1,1) = 4.631579e-01;
  hex8_xyze(2,1) = 1.875000e-01;
  hex8_xyze(0,2) = 1.021277e+00;
  hex8_xyze(1,2) = 4.631579e-01;
  hex8_xyze(2,2) = 1.875000e-01;
  hex8_xyze(0,3) = 1.021277e+00;
  hex8_xyze(1,3) = 5.052632e-01;
  hex8_xyze(2,3) = 1.875000e-01;
  hex8_xyze(0,4) = 1.063830e+00;
  hex8_xyze(1,4) = 5.052632e-01;
  hex8_xyze(2,4) = 1.562500e-01;
  hex8_xyze(0,5) = 1.063830e+00;
  hex8_xyze(1,5) = 4.631579e-01;
  hex8_xyze(2,5) = 1.562500e-01;
  hex8_xyze(0,6) = 1.021277e+00;
  hex8_xyze(1,6) = 4.631579e-01;
  hex8_xyze(2,6) = 1.562500e-01;
  hex8_xyze(0,7) = 1.021277e+00;
  hex8_xyze(1,7) = 5.052632e-01;
  hex8_xyze(2,7) = 1.562500e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex37()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 3 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

  std::map<std::string, int> nodeids;

  nxyz1(9.164261e-01,2.517463e-01,3.209446e-01);
  nxyz2(9.152678e-01,2.201901e-01,3.210346e-01);
  nxyz3(9.158694e-01,2.361631e-01,2.933258e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(8.899081e-01,2.626309e-01,3.209111e-01);
  nxyz2(8.888432e-01,2.311533e-01,3.209298e-01);
  nxyz3(9.164261e-01,2.517463e-01,3.209446e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(8.888432e-01,2.311533e-01,3.209298e-01);
  nxyz2(9.152678e-01,2.201901e-01,3.210346e-01);
  nxyz3(9.164261e-01,2.517463e-01,3.209446e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(9.164261e-01,2.517463e-01,3.209446e-01);
  nxyz2(9.158694e-01,2.361631e-01,2.933258e-01);
  nxyz3(9.170605e-01,2.682050e-01,2.929160e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(8.899081e-01,2.626309e-01,3.209111e-01);
  nxyz2(9.164261e-01,2.517463e-01,3.209446e-01);
  nxyz3(9.176455e-01,2.833532e-01,3.209024e-01);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(9.176455e-01,2.833532e-01,3.209024e-01);
  nxyz2(9.164261e-01,2.517463e-01,3.209446e-01);
  nxyz3(9.170605e-01,2.682050e-01,2.929160e-01);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::tri3 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 9.104478e-01;
  hex8_xyze(1,0) = 2.518519e-01;
  hex8_xyze(2,0) = 2.941176e-01;
  hex8_xyze(0,1) = 9.104478e-01;
  hex8_xyze(1,1) = 2.444444e-01;
  hex8_xyze(2,1) = 2.941176e-01;
  hex8_xyze(0,2) = 9.179104e-01;
  hex8_xyze(1,2) = 2.444444e-01;
  hex8_xyze(2,2) = 2.941176e-01;
  hex8_xyze(0,3) = 9.179104e-01;
  hex8_xyze(1,3) = 2.518519e-01;
  hex8_xyze(2,3) = 2.941176e-01;
  hex8_xyze(0,4) = 9.104478e-01;
  hex8_xyze(1,4) = 2.518519e-01;
  hex8_xyze(2,4) = 3.235294e-01;
  hex8_xyze(0,5) = 9.104478e-01;
  hex8_xyze(1,5) = 2.444444e-01;
  hex8_xyze(2,5) = 3.235294e-01;
  hex8_xyze(0,6) = 9.179104e-01;
  hex8_xyze(1,6) = 2.444444e-01;
  hex8_xyze(2,6) = 3.235294e-01;
  hex8_xyze(0,7) = 9.179104e-01;
  hex8_xyze(1,7) = 2.518519e-01;
  hex8_xyze(2,7) = 3.235294e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_alex38()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(9.177310e-01,3.765734e-01,2.887735e-01);
  nxyz2(9.178660e-01,3.765640e-01,3.208629e-01);
  nxyz3(9.444573e-01,3.733318e-01,3.207611e-01);
  nxyz4(9.443162e-01,3.733500e-01,2.886787e-01);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.405980e-01,3.419116e-01,3.207990e-01);
  nxyz2(9.404867e-01,3.419193e-01,2.887009e-01);
  nxyz3(9.443162e-01,3.733500e-01,2.886787e-01);
  nxyz4(9.444573e-01,3.733318e-01,3.207611e-01);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.444571e-01,4.266682e-01,3.207612e-01);
  nxyz2(9.443160e-01,4.266500e-01,2.886787e-01);
  nxyz3(9.404865e-01,4.580807e-01,2.887009e-01);
  nxyz4(9.405979e-01,4.580884e-01,3.207990e-01);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.178659e-01,4.234360e-01,3.208629e-01);
  nxyz2(9.444571e-01,4.266682e-01,3.207612e-01);
  nxyz3(9.405979e-01,4.580884e-01,3.207990e-01);
  nxyz4(9.140159e-01,4.548570e-01,3.208757e-01);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.178659e-01,4.234360e-01,3.208629e-01);
  nxyz2(9.177309e-01,4.234266e-01,2.887735e-01);
  nxyz3(9.443160e-01,4.266500e-01,2.886787e-01);
  nxyz4(9.444571e-01,4.266682e-01,3.207612e-01);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(9.140160e-01,3.451430e-01,3.208757e-01);
  nxyz2(9.405980e-01,3.419116e-01,3.207990e-01);
  nxyz3(9.444573e-01,3.733318e-01,3.207611e-01);
  nxyz4(9.178660e-01,3.765640e-01,3.208629e-01);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1.000000e+00;
  hex8_xyze(1,0) = 4.266667e-01;
  hex8_xyze(2,0) = 3.333333e-01;
  hex8_xyze(0,1) = 1.000000e+00;
  hex8_xyze(1,1) = 3.733333e-01;
  hex8_xyze(2,1) = 3.333333e-01;
  hex8_xyze(0,2) = 9.444444e-01;
  hex8_xyze(1,2) = 3.733333e-01;
  hex8_xyze(2,2) = 3.333333e-01;
  hex8_xyze(0,3) = 9.444444e-01;
  hex8_xyze(1,3) = 4.266667e-01;
  hex8_xyze(2,3) = 3.333333e-01;
  hex8_xyze(0,4) = 1.000000e+00;
  hex8_xyze(1,4) = 4.266667e-01;
  hex8_xyze(2,4) = 2.916667e-01;
  hex8_xyze(0,5) = 1.000000e+00;
  hex8_xyze(1,5) = 3.733333e-01;
  hex8_xyze(2,5) = 2.916667e-01;
  hex8_xyze(0,6) = 9.444444e-01;
  hex8_xyze(1,6) = 3.733333e-01;
  hex8_xyze(2,6) = 2.916667e-01;
  hex8_xyze(0,7) = 9.444444e-01;
  hex8_xyze(1,7) = 4.266667e-01;
  hex8_xyze(2,7) = 2.916667e-01;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_tri3_ursula1()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 3 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

  std::map<std::string, int> nodeids;

  nxyz1( -1, -1, 0.99999999999999933 );
  nxyz2( 1, -1, 0.99999999999999933 );
  nxyz3( 0, -0.99999999999999989, 0.99999999999999967 );

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1( 1, -1, 0.99999999999999933 );
  nxyz2( 1, -0.99999999999999978, 1 );
  nxyz3( 0, -0.99999999999999989, 0.99999999999999967 );

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1( 1, -0.99999999999999978, 1 );
  nxyz2( -1, -0.99999999999999978, 1 );
  nxyz3( 0, -0.99999999999999989, 0.99999999999999967 );

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1( -1, -0.99999999999999978, 1 );
  nxyz2( -1, -1, 0.99999999999999933 );
  nxyz3( 0, -0.99999999999999989, 0.99999999999999967 );

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::tri3 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = -1;
  hex8_xyze(1,0) = -1;
  hex8_xyze(2,0) = -1;
  hex8_xyze(0,1) = 1;
  hex8_xyze(1,1) = -1;
  hex8_xyze(2,1) = -1;
  hex8_xyze(0,2) = 1;
  hex8_xyze(1,2) = 1;
  hex8_xyze(2,2) = -1;
  hex8_xyze(0,3) = -1;
  hex8_xyze(1,3) = 1;
  hex8_xyze(2,3) = -1;
  hex8_xyze(0,4) = -1;
  hex8_xyze(1,4) = -1;
  hex8_xyze(2,4) = 1;
  hex8_xyze(0,5) = 1;
  hex8_xyze(1,5) = -1;
  hex8_xyze(2,5) = 1;
  hex8_xyze(0,6) = 1;
  hex8_xyze(1,6) = 1;
  hex8_xyze(2,6) = 1;
  hex8_xyze(0,7) = -1;
  hex8_xyze(1,7) = 1;
  hex8_xyze(2,7) = 1;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel7()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 3 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

  std::map<std::string, int> nodeids;

  int sideid = 0;

  nxyz1(1.0992740341033049312,0.23760272423331274538,0.0050000250500000006232);
  nxyz2(1.0909091200000002875,0.21889414848398383584,0.010000050099999999512);
  nxyz3(1.107676579371152048,0.25615896479243038808,0.010000050100000001246);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341033049312,0.23760272423331274538,0.0050000250500000006232);
  nxyz2(1.107676579371152048,0.25615896479243038808,0.010000050100000001246);
  nxyz3(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341033049312,0.23760272423331274538,0.0050000250500000006232);
  nxyz2(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);
  nxyz3(1.0909091200000000654,0.21915380243406717975,4.4629817099252685055e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341033049312,0.23760272423331274538,0.0050000250500000006232);
  nxyz2(1.0909091200000000654,0.21915380243406717975,4.4629817099252685055e-19);
  nxyz3(1.0909091200000002875,0.21889414848398383584,0.010000050099999999512);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341020581507,0.23760272423108014239,-0.0050000250499999997558);
  nxyz2(1.0909091200000000654,0.21915380243406717975,4.4629817099252685055e-19);
  nxyz3(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341020581507,0.23760272423108014239,-0.0050000250499999997558);
  nxyz2(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);
  nxyz3(1.1076765793661651482,0.25615896478304178707,-0.010000050099999999512);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341020581507,0.23760272423108014239,-0.0050000250499999997558);
  nxyz2(1.1076765793661651482,0.25615896478304178707,-0.010000050099999999512);
  nxyz3(1.0909091200000000654,0.2188941484844418861,-0.010000050099999999512);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.0992740341020581507,0.23760272423108014239,-0.0050000250499999997558);
  nxyz2(1.0909091200000000654,0.2188941484844418861,-0.010000050099999999512);
  nxyz3(1.0909091200000000654,0.21915380243406717975,4.4629817099252685055e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.1116124183961666461,0.2644543770038000341,0.0050000250500000006232);
  nxyz2(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);
  nxyz3(1.107676579371152048,0.25615896479243038808,0.010000050100000001246);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.1116124183961666461,0.2644543770038000341,0.0050000250500000006232);
  nxyz2(1.1155384155882137609,0.27272728099999998808,3.4966861014134080861e-19);
  nxyz3(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.111612418394896773,0.26445437700145285609,-0.0050000250499999997558);
  nxyz2(1.1076765793661651482,0.25615896478304178707,-0.010000050099999999512);
  nxyz3(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(1.111612418394896773,0.26445437700145285609,-0.0050000250499999997558);
  nxyz2(1.1076013170420673237,0.25620398122276971664,3.8158834299999999373e-19);
  nxyz3(1.1155384155882137609,0.27272728099999998808,3.4966861014134080861e-19);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::tri3 );

  Epetra_SerialDenseMatrix tet4_xyze( 3, 8 );

  tet4_xyze(0,0) = 1.0992740341033049312;
  tet4_xyze(1,0) = 0.23760272423331274538;
  tet4_xyze(2,0) = 0.0050000250500000006232;
  tet4_xyze(0,1) = 1.1076013170420673237;
  tet4_xyze(1,1) = 0.25620398122276971664;
  tet4_xyze(2,1) = 3.8158834299999999373e-19;
  tet4_xyze(0,2) = 1.0992740341020581507;
  tet4_xyze(1,2) = 0.23760272423108014239;
  tet4_xyze(2,2) = -0.0050000250499999997558;
  tet4_xyze(0,3) = 1.0909091200000000654;
  tet4_xyze(1,3) = 0.27272728099999998808;
  tet4_xyze(2,3) = -0.010000050099999999512;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel6()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0.72658695971228337562,0.46225225295206362475,0.015000075099999999467);
  nxyz2(0.75257459164385520634,0.54643521666375061141,0.015000075099999999467);
  nxyz3(0.75301582121389776248,0.54642941846884562906,3.1223752000000001688e-19);
  nxyz4(0.7269042669886265351,0.46207429315320791563,1.9079417100000000732e-19);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.70070698096704708835,0.37606951015675532668,0.015000075099999999467);
  nxyz2(0.72658695971228337562,0.46225225295206362475,0.015000075099999999467);
  nxyz3(0.7269042669886265351,0.46207429315320791563,1.9079417100000000732e-19);
  nxyz4(0.70107331434079600552,0.37602569613385122826,5.0313757100000004007e-19);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.70107331434079600552,0.37602569613385122826,5.0313757100000004007e-19);
  nxyz2(0.7269042669886265351,0.46207429315320791563,1.9079417100000000732e-19);
  nxyz3(0.72658695972323428247,0.4622522529408981673,-0.015000075099999999467);
  nxyz4(0.70070698096458938764,0.37606951016009454447,-0.015000075099999999467);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0.7269042669886265351,0.46207429315320791563,1.9079417100000000732e-19);
  nxyz2(0.75301582121389776248,0.54642941846884562906,3.1223752000000001688e-19);
  nxyz3(0.75257459166386664329,0.5464352166907603392,-0.015000075099999999467);
  nxyz4(0.72658695972323428247,0.4622522529408981673,-0.015000075099999999467);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.72727274900000005164;
  hex8_xyze(1,0) = 0.54545456199999997615;
  hex8_xyze(2,0) = 0.010000050099999999512;
  hex8_xyze(0,1) = 0.72727274900000005164;
  hex8_xyze(1,1) = 0.45454546800000000806;
  hex8_xyze(2,1) = 0.010000050099999999512;
  hex8_xyze(0,2) = 0.63636362599999996004;
  hex8_xyze(1,2) = 0.45454546800000000806;
  hex8_xyze(2,2) = 0.010000050099999999512;
  hex8_xyze(0,3) = 0.63636362599999996004;
  hex8_xyze(1,3) = 0.54545456199999997615;
  hex8_xyze(2,3) = 0.010000050099999999512;
  hex8_xyze(0,4) = 0.72727274900000005164;
  hex8_xyze(1,4) = 0.54545456199999997615;
  hex8_xyze(2,4) = -0.010000050099999999512;
  hex8_xyze(0,5) = 0.72727274900000005164;
  hex8_xyze(1,5) = 0.45454546800000000806;
  hex8_xyze(2,5) = -0.010000050099999999512;
  hex8_xyze(0,6) = 0.63636362599999996004;
  hex8_xyze(1,6) = 0.45454546800000000806;
  hex8_xyze(2,6) = -0.010000050099999999512;
  hex8_xyze(0,7) = 0.63636362599999996004;
  hex8_xyze(1,7) = 0.54545456199999997615;
  hex8_xyze(2,7) = -0.010000050099999999512;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel5()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.70107881777883, 0.376288005548294, 0.0150000751);
  nxyz2( 0.727082769140639, 0.462499964659154, 0.0150000751);
  nxyz3( 0.727403031742602, 0.462321135243137, 1.90794171e-19);
  nxyz4( 0.701446387019269, 0.376242724613422, 5.03137571e-19);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.727082769140639, 0.462499964659154, 0.0150000751);
  nxyz2( 0.75320276810349, 0.546724469276745, 0.0150000751);
  nxyz3( 0.753645226950538, 0.546716705138847, 3.1223752e-19);
  nxyz4( 0.727403031742602, 0.462321135243137, 1.90794171e-19);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.701446387019269, 0.376242724613422, 5.03137571e-19);
  nxyz2( 0.727403031742602, 0.462321135243137, 1.90794171e-19);
  nxyz3( 0.727082769041564, 0.462499964711039,-0.0150000751);
  nxyz4( 0.701078817808441, 0.376288005536189,-0.0150000751);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.727403031742602, 0.462321135243137, 1.90794171e-19);
  nxyz2( 0.753645226950538, 0.546716705138847, 3.1223752e-19);
  nxyz3( 0.753202768114659, 0.546724469252098,-0.0150000751);
  nxyz4( 0.727082769041564, 0.462499964711039,-0.0150000751);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.818181813;
  hex8_xyze(1,0) = 0.545454562;
  hex8_xyze(2,0) = 0.0100000501;
  hex8_xyze(0,1) = 0.818181813;
  hex8_xyze(1,1) = 0.454545468;
  hex8_xyze(2,1) = 0.0100000501;
  hex8_xyze(0,2) = 0.727272749;
  hex8_xyze(1,2) = 0.454545468;
  hex8_xyze(2,2) = 0.0100000501;
  hex8_xyze(0,3) = 0.727272749;
  hex8_xyze(1,3) = 0.545454562;
  hex8_xyze(2,3) = 0.0100000501;
  hex8_xyze(0,4) = 0.818181813;
  hex8_xyze(1,4) = 0.545454562;
  hex8_xyze(2,4) = -0.0100000501;
  hex8_xyze(0,5) = 0.818181813;
  hex8_xyze(1,5) = 0.454545468;
  hex8_xyze(2,5) = -0.0100000501;
  hex8_xyze(0,6) = 0.727272749;
  hex8_xyze(1,6) = 0.454545468;
  hex8_xyze(2,6) = -0.0100000501;
  hex8_xyze(0,7) = 0.727272749;
  hex8_xyze(1,7) = 0.545454562;
  hex8_xyze(2,7) = -0.0100000501;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true , INPAR::CUT::VCellGaussPts_Tessellation);
  intersection.Status();
}

void test_hex8_quad4_axel4()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.701073314340809, 0.376025696133856, 5.03137571e-19);
  nxyz2( 0.726904266988644, 0.462074293153213, 1.90794171e-19);
  nxyz3( 0.726586959713632, 0.462252252949456,-0.0150000751);
  nxyz4( 0.700706980965754, 0.376069510157904,-0.0150000751);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.726904266988644, 0.462074293153213, 1.90794171e-19);
  nxyz2( 0.753015821213919, 0.546429418468849, 3.1223752e-19);
  nxyz3( 0.752574591649044, 0.5464352166732,-0.0150000751);
  nxyz4( 0.726586959713632, 0.462252252949456,-0.0150000751);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.700706980965909, 0.376069510158957, 0.0150000751);
  nxyz2( 0.726586959721921, 0.462252252943515, 0.0150000751);
  nxyz3( 0.726904266988644, 0.462074293153213, 1.90794171e-19);
  nxyz4( 0.701073314340809, 0.376025696133856, 5.03137571e-19);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.726586959721921, 0.462252252943515, 0.0150000751);
  nxyz2( 0.75257459165872, 0.546435216681318, 0.0150000751);
  nxyz3( 0.753015821213919, 0.546429418468849, 3.1223752e-19);
  nxyz4( 0.726904266988644, 0.462074293153213, 1.90794171e-19);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.727272749;
  hex8_xyze(1,0) = 0.545454562;
  hex8_xyze(2,0) = 0.0100000501;
  hex8_xyze(0,1) = 0.727272749;
  hex8_xyze(1,1) = 0.454545468;
  hex8_xyze(2,1) = 0.0100000501;
  hex8_xyze(0,2) = 0.636363626;
  hex8_xyze(1,2) = 0.454545468;
  hex8_xyze(2,2) = 0.0100000501;
  hex8_xyze(0,3) = 0.636363626;
  hex8_xyze(1,3) = 0.545454562;
  hex8_xyze(2,3) = 0.0100000501;
  hex8_xyze(0,4) = 0.727272749;
  hex8_xyze(1,4) = 0.545454562;
  hex8_xyze(2,4) = -0.0100000501;
  hex8_xyze(0,5) = 0.727272749;
  hex8_xyze(1,5) = 0.454545468;
  hex8_xyze(2,5) = -0.0100000501;
  hex8_xyze(0,6) = 0.636363626;
  hex8_xyze(1,6) = 0.454545468;
  hex8_xyze(2,6) = -0.0100000501;
  hex8_xyze(0,7) = 0.636363626;
  hex8_xyze(1,7) = 0.545454562;
  hex8_xyze(2,7) = -0.0100000501;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel3()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.727052975423276, 0.462358773566033, 0.0150000751);
  nxyz2( 0.75313893222712, 0.546607019429788, 0.0150000751);
  nxyz3( 0.753577731984463, 0.546598870509341, 3.1223752e-19);
  nxyz4( 0.727374824073424, 0.462177278436167, 1.90794171e-19);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.701081603684843, 0.37616332593259, 0.0150000751);
  nxyz2( 0.727052975423276, 0.462358773566033, 0.0150000751);
  nxyz3( 0.727374824073424, 0.462177278436167, 1.90794171e-19);
  nxyz4( 0.701449546070791, 0.376118628569801, 5.03137571e-19);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.701449546070791, 0.376118628569801, 5.03137571e-19);
  nxyz2( 0.727374824073424, 0.462177278436167, 1.90794171e-19);
  nxyz3( 0.727052975413278, 0.462358773546029,-0.0150000751);
  nxyz4( 0.701081603687067, 0.376163325942043,-0.0150000751);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.727374824073424, 0.462177278436167, 1.90794171e-19);
  nxyz2( 0.753577731984463, 0.546598870509341, 3.1223752e-19);
  nxyz3( 0.75313893225359, 0.546607019469259,-0.0150000751);
  nxyz4( 0.727052975413278, 0.462358773546029,-0.0150000751);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.818181813;
  hex8_xyze(1,0) = 0.545454562;
  hex8_xyze(2,0) = 0.0100000501;
  hex8_xyze(0,1) = 0.818181813;
  hex8_xyze(1,1) = 0.454545468;
  hex8_xyze(2,1) = 0.0100000501;
  hex8_xyze(0,2) = 0.727272749;
  hex8_xyze(1,2) = 0.454545468;
  hex8_xyze(2,2) = 0.0100000501;
  hex8_xyze(0,3) = 0.727272749;
  hex8_xyze(1,3) = 0.545454562;
  hex8_xyze(2,3) = 0.0100000501;
  hex8_xyze(0,4) = 0.818181813;
  hex8_xyze(1,4) = 0.545454562;
  hex8_xyze(2,4) = -0.0100000501;
  hex8_xyze(0,5) = 0.818181813;
  hex8_xyze(1,5) = 0.454545468;
  hex8_xyze(2,5) = -0.0100000501;
  hex8_xyze(0,6) = 0.727272749;
  hex8_xyze(1,6) = 0.454545468;
  hex8_xyze(2,6) = -0.0100000501;
  hex8_xyze(0,7) = 0.727272749;
  hex8_xyze(1,7) = 0.545454562;
  hex8_xyze(2,7) = -0.0100000501;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel2()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 1.13739619241744, 0.4378327878848, 0.0150000751);
  nxyz2( 1.13725580113275, 0.437863771860061, 1.38913403e-19);
  nxyz3( 1.16330399014196, 0.518279763846018, 0);
  nxyz4( 1.16335802489377, 0.518256842898359, 0.0150000751);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.13725580113275, 0.437863771860061, 1.38913403e-19);
  nxyz2( 1.13739619242228, 0.43783278790288,-0.0150000751);
  nxyz3( 1.16335802490459, 0.518256842856954,-0.0150000751);
  nxyz4( 1.16330399014196, 0.518279763846018, 0);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.08635704859813, 0.541649676999336, 5.16901856e-19);
  nxyz2( 1.08625184048568, 0.541799217272295, 0.0150000751);
  nxyz3( 1.16335802489377, 0.518256842898359, 0.0150000751);
  nxyz4( 1.16330399014196, 0.518279763846018, 0);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.08625184047561, 0.541799217272842,-0.0150000751);
  nxyz2( 1.08635704859813, 0.541649676999336, 5.16901856e-19);
  nxyz3( 1.16330399014196, 0.518279763846018, 0);
  nxyz4( 1.16335802490459, 0.518256842856954,-0.0150000751);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1.18181813;
  hex8_xyze(1,0) = 0.545454562;
  hex8_xyze(2,0) = 0.0100000501;
  hex8_xyze(0,1) = 1.18181813;
  hex8_xyze(1,1) = 0.454545468;
  hex8_xyze(2,1) = 0.0100000501;
  hex8_xyze(0,2) = 1.09090912;
  hex8_xyze(1,2) = 0.454545468;
  hex8_xyze(2,2) = 0.0100000501;
  hex8_xyze(0,3) = 1.09090912;
  hex8_xyze(1,3) = 0.545454562;
  hex8_xyze(2,3) = 0.0100000501;
  hex8_xyze(0,4) = 1.18181813;
  hex8_xyze(1,4) = 0.545454562;
  hex8_xyze(2,4) = -0.0100000501;
  hex8_xyze(0,5) = 1.18181813;
  hex8_xyze(1,5) = 0.454545468;
  hex8_xyze(2,5) = -0.0100000501;
  hex8_xyze(0,6) = 1.09090912;
  hex8_xyze(1,6) = 0.454545468;
  hex8_xyze(2,6) = -0.0100000501;
  hex8_xyze(0,7) = 1.09090912;
  hex8_xyze(1,7) = 0.545454562;
  hex8_xyze(2,7) = -0.0100000501;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_axel1()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 1.14226661309652, 0.511653661871148, 5.16901856e-19);
  nxyz2( 1.14217876340323, 0.511845001037696, 0.0150000751);
  nxyz3( 1.2159743756806, 0.477514442783813, 0.0150000751);
  nxyz4( 1.21595128604308, 0.477503887225635, 0);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.14217876336518, 0.511845001050777,-0.0150000751);
  nxyz2( 1.14226661309652, 0.511653661871148, 5.16901856e-19);
  nxyz3( 1.21595128604308, 0.477503887225635, 0);
  nxyz4( 1.21597437572149, 0.477514442626296,-0.0150000751);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.17834292908302, 0.402093972352961, 0.0150000751);
  nxyz2( 1.17825626213532, 0.402150428778217, 1.38913403e-19);
  nxyz3( 1.21595128604308, 0.477503887225635, 0);
  nxyz4( 1.2159743756806, 0.477514442783813, 0.0150000751);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.17825626213532, 0.402150428778217, 1.38913403e-19);
  nxyz2( 1.17834292911377, 0.402093972415847,-0.0150000751);
  nxyz3( 1.21597437572149, 0.477514442626296,-0.0150000751);
  nxyz4( 1.21595128604308, 0.477503887225635, 0);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1.27272725;
  hex8_xyze(1,0) = 0.545454562;
  hex8_xyze(2,0) = 0.0100000501;
  hex8_xyze(0,1) = 1.27272725;
  hex8_xyze(1,1) = 0.454545468;
  hex8_xyze(2,1) = 0.0100000501;
  hex8_xyze(0,2) = 1.18181813;
  hex8_xyze(1,2) = 0.454545468;
  hex8_xyze(2,2) = 0.0100000501;
  hex8_xyze(0,3) = 1.18181813;
  hex8_xyze(1,3) = 0.545454562;
  hex8_xyze(2,3) = 0.0100000501;
  hex8_xyze(0,4) = 1.27272725;
  hex8_xyze(1,4) = 0.545454562;
  hex8_xyze(2,4) = -0.0100000501;
  hex8_xyze(0,5) = 1.27272725;
  hex8_xyze(1,5) = 0.454545468;
  hex8_xyze(2,5) = -0.0100000501;
  hex8_xyze(0,6) = 1.18181813;
  hex8_xyze(1,6) = 0.454545468;
  hex8_xyze(2,6) = -0.0100000501;
  hex8_xyze(0,7) = 1.18181813;
  hex8_xyze(1,7) = 0.545454562;
  hex8_xyze(2,7) = -0.0100000501;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true, INPAR::CUT::VCellGaussPts_DirectDivergence );
  intersection.Status();
}

void test_hex8_quad4_shadan5()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;
  int sideid = 0;

  nxyz1(-0.086426593400000001344,-0.14840835299999999286,0.1050907819999999937);
  nxyz2(-0.073279127499999999196,-0.09286977350000000242,0.13056352700000001255);
  nxyz3(-0.098751873774999984756,-0.064204173160000002629,0.081210988075000004049);
  nxyz4(-0.11189933900000000044,-0.1197427509999999945,0.05573824420000000196);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.073279127499999999196,-0.09286977350000000242,0.13056352700000001255);
  nxyz2(-0.060131665299999997876,-0.037331193700000001034,0.15603627300000000311);
  nxyz3(-0.085604406899999999081,-0.00866559706999999943,0.10668373100000000386);
  nxyz4(-0.098751873774999984756,-0.064204173160000002629,0.081210988075000004049);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.086426593400000001344,-0.14840835299999999286,0.1050907819999999937);
  nxyz2(-0.030888017300000000787,-0.14840835299999999286,0.076425179800000006547);
  nxyz3(-0.0177405521200000027,-0.09286977350000000242,0.10189792465000002686);
  nxyz4(-0.073279127499999999196,-0.09286977350000000242,0.13056352700000001255);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0177405521200000027,-0.09286977350000000242,0.10189792465000002686);
  nxyz2(0.037798024700000001919,-0.09286977350000000242,0.073232330400000003601);
  nxyz3(0.050945486900000003239,-0.037331193700000001034,0.098705075700000005634);
  nxyz4(-0.0045930896899999996974,-0.037331193700000001034,0.1273706699999999914);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.073279127499999999196,-0.09286977350000000242,0.13056352700000001255);
  nxyz2(-0.0177405521200000027,-0.09286977350000000242,0.10189792465000002686);
  nxyz3(-0.0045930896899999996974,-0.037331193700000001034,0.1273706699999999914);
  nxyz4(-0.060131665299999997876,-0.037331193700000001034,0.15603627300000000311);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.060131665299999997876,-0.037331193700000001034,0.15603627300000000311);
  nxyz2(-0.046984203199999997891,0.018207382399999999523,0.1815090179999999942);
  nxyz3(-0.072456946699999996708,0.046872977115000005743,0.13215647934999999458);
  nxyz4(-0.085604406899999999081,-0.00866559706999999943,0.10668373100000000386);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.060131665299999997876,-0.037331193700000001034,0.15603627300000000311);
  nxyz2(-0.0045930896899999996974,-0.037331193700000001034,0.1273706699999999914);
  nxyz3(0.0085543726800000018107,0.018207382399999999523,0.1528434175749999735);
  nxyz4(-0.046984203199999997891,0.018207382399999999523,0.1815090179999999942);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0045930896899999996974,-0.037331193700000001034,0.1273706699999999914);
  nxyz2(0.050945486900000003239,-0.037331193700000001034,0.098705075700000005634);
  nxyz3(0.064092948999999996285,0.018207382399999999523,0.12417781399999999747);
  nxyz4(0.0085543726800000018107,0.018207382399999999523,0.1528434175749999735);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = -0.125;
  hex8_xyze(1,0) = -0.125;
  hex8_xyze(2,0) = 0.25;
  hex8_xyze(0,1) = -0.125;
  hex8_xyze(1,1) = -0.125;
  hex8_xyze(2,1) = 0.125;
  hex8_xyze(0,2) = -0.125;
  hex8_xyze(1,2) = 0;
  hex8_xyze(2,2) = 0.125;
  hex8_xyze(0,3) = -0.125;
  hex8_xyze(1,3) = 0;
  hex8_xyze(2,3) = 0.25;
  hex8_xyze(0,4) = 0;
  hex8_xyze(1,4) = -0.125;
  hex8_xyze(2,4) = 0.25;
  hex8_xyze(0,5) = 0;
  hex8_xyze(1,5) = -0.125;
  hex8_xyze(2,5) = 0.125;
  hex8_xyze(0,6) = 0;
  hex8_xyze(1,6) = 0;
  hex8_xyze(2,6) = 0.125;
  hex8_xyze(0,7) = 0;
  hex8_xyze(1,7) = 0;
  hex8_xyze(2,7) = 0.25;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_shadan4()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;
  int sideid = 0;

  nxyz1(-0.0308880173,-0.148408353, 0.0764251798);
  nxyz2(-0.056360762772,-0.119742751, 0.0270726447);
  nxyz3(-0.000822183094,-0.119742751,-0.00159295055);
  nxyz4( 0.0246505607,-0.148408353, 0.0477595851);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0308880173,-0.148408353, 0.0764251798);
  nxyz2( 0.0246505607,-0.148408353, 0.0477595851);
  nxyz3( 0.0377980247,-0.0928697735, 0.0732323304);
  nxyz4(-0.01774055212,-0.0928697735, 0.10189792465);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.01774055212,-0.0928697735, 0.10189792465);
  nxyz2( 0.0377980247,-0.0928697735, 0.0732323304);
  nxyz3( 0.0509454869,-0.0373311937, 0.0987050757);
  nxyz4(-0.00459308969,-0.0373311937, 0.12737067);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.00459308969,-0.0373311937, 0.12737067);
  nxyz2( 0.0509454869,-0.0373311937, 0.0987050757);
  nxyz3( 0.064092949, 0.0182073824, 0.124177814);
  nxyz4( 0.00855437268, 0.0182073824, 0.152843417575);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.0246505607,-0.148408353, 0.0477595851);
  nxyz2(-0.000822183094,-0.119742751,-0.00159295055);
  nxyz3( 0.054716394153,-0.119742751,-0.030258549235);
  nxyz4( 0.0801891387,-0.148408353, 0.0190939885);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.0246505607,-0.148408353, 0.0477595851);
  nxyz2( 0.0801891387,-0.148408353, 0.0190939885);
  nxyz3( 0.09333660155,-0.0928697735, 0.044566731015);
  nxyz4( 0.0377980247,-0.0928697735, 0.0732323304);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.0801891387,-0.148408353, 0.0190939885);
  nxyz2( 0.135727718,-0.148408353,-0.00957160816);
  nxyz3( 0.148875177,-0.0928697735, 0.0159011353);
  nxyz4( 0.09333660155,-0.0928697735, 0.044566731015);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.09333660155,-0.0928697735, 0.044566731015);
  nxyz2( 0.148875177,-0.0928697735, 0.0159011353);
  nxyz3( 0.162022635,-0.0373311937, 0.0413738787);
  nxyz4( 0.106484063,-0.0373311937, 0.0700394735);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.0377980247,-0.0928697735, 0.0732323304);
  nxyz2( 0.09333660155,-0.0928697735, 0.044566731015);
  nxyz3( 0.106484063,-0.0373311937, 0.0700394735);
  nxyz4( 0.0509454869,-0.0373311937, 0.0987050757);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.0509454869,-0.0373311937, 0.0987050757);
  nxyz2( 0.106484063,-0.0373311937, 0.0700394735);
  nxyz3( 0.11963152805, 0.0182073824, 0.095512217075);
  nxyz4( 0.064092949, 0.0182073824, 0.124177814);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.106484063,-0.0373311937, 0.0700394735);
  nxyz2( 0.162022635,-0.0373311937, 0.0413738787);
  nxyz3( 0.175170109, 0.0182073824, 0.066846624);
  nxyz4( 0.11963152805, 0.0182073824, 0.095512217075);

  intersection.AddCutSide( ++sideid, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0;
  hex8_xyze(1,0) = -0.125;
  hex8_xyze(2,0) = 0.125;
  hex8_xyze(0,1) = 0;
  hex8_xyze(1,1) = -0.125;
  hex8_xyze(2,1) = 0;
  hex8_xyze(0,2) = 0;
  hex8_xyze(1,2) = 0;
  hex8_xyze(2,2) = 0;
  hex8_xyze(0,3) = 0;
  hex8_xyze(1,3) = 0;
  hex8_xyze(2,3) = 0.125;
  hex8_xyze(0,4) = 0.125;
  hex8_xyze(1,4) = -0.125;
  hex8_xyze(2,4) = 0.125;
  hex8_xyze(0,5) = 0.125;
  hex8_xyze(1,5) = -0.125;
  hex8_xyze(2,5) = 0;
  hex8_xyze(0,6) = 0.125;
  hex8_xyze(1,6) = 0;
  hex8_xyze(2,6) = 0;
  hex8_xyze(0,7) = 0.125;
  hex8_xyze(1,7) = 0;
  hex8_xyze(2,7) = 0.125;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_shadan3()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.952719331, 0.328301698, 0.49193272);
  nxyz2( 0.968041003, 0.328301698, 0.553384781);
  nxyz3( 0.926130831, 0.374620765, 0.56383419);
  nxyz4( 0.9108091, 0.374620765, 0.5023821);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.952719331, 0.328301698, 0.49193272);
  nxyz2( 1.00678062, 0.302910298, 0.478453726);
  nxyz3( 1.02210236, 0.302910298, 0.539905787);
  nxyz4( 0.968041003, 0.328301698, 0.553384781);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.937397599, 0.328301698, 0.430480659);
  nxyz2( 0.952719331, 0.328301698, 0.49193272);
  nxyz3( 0.9108091, 0.374620765, 0.5023821);
  nxyz4( 0.895487368, 0.374620765, 0.440930039);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.937397599, 0.328301698, 0.430480659);
  nxyz2( 0.991458893, 0.302910298, 0.417001665);
  nxyz3( 1.00678062, 0.302910298, 0.478453726);
  nxyz4( 0.952719331, 0.328301698, 0.49193272);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.940340042, 0.344281912, 0.371446699);
  nxyz2( 0.937397599, 0.328301698, 0.430480659);
  nxyz3( 0.895487368, 0.374620765, 0.440930039);
  nxyz4( 0.898429871, 0.390600979, 0.381896079);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.937397599, 0.328301698, 0.430480659);
  nxyz2( 0.940340042, 0.344281912, 0.371446699);
  nxyz3( 0.980782151, 0.317581624, 0.37206167);
  nxyz4( 0.991458893, 0.302910298, 0.417001665);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.991458893, 0.302910298, 0.417001665);
  nxyz2( 1.0471071, 0.324353248, 0.403126985);
  nxyz3( 1.06242883, 0.324353248, 0.464579046);
  nxyz4( 1.00678062, 0.302910298, 0.478453726);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.0471071, 0.324353248, 0.403126985);
  nxyz2( 0.991458893, 0.302910298, 0.417001665);
  nxyz3( 0.980782151, 0.317581624, 0.37206167);
  nxyz4( 1.01791644, 0.341489941, 0.352104753);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.916666687;
  hex8_xyze(1,0) = 0.25;
  hex8_xyze(2,0) = 0.5;
  hex8_xyze(0,1) = 0.916666687;
  hex8_xyze(1,1) = 0.25;
  hex8_xyze(2,1) = 0.416666657;
  hex8_xyze(0,2) = 0.916666687;
  hex8_xyze(1,2) = 0.333333343;
  hex8_xyze(2,2) = 0.416666657;
  hex8_xyze(0,3) = 0.916666687;
  hex8_xyze(1,3) = 0.333333343;
  hex8_xyze(2,3) = 0.5;
  hex8_xyze(0,4) = 1;
  hex8_xyze(1,4) = 0.25;
  hex8_xyze(2,4) = 0.5;
  hex8_xyze(0,5) = 1;
  hex8_xyze(1,5) = 0.25;
  hex8_xyze(2,5) = 0.416666657;
  hex8_xyze(0,6) = 1;
  hex8_xyze(1,6) = 0.333333343;
  hex8_xyze(2,6) = 0.416666657;
  hex8_xyze(0,7) = 1;
  hex8_xyze(1,7) = 0.333333343;
  hex8_xyze(2,7) = 0.5;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_shadan2()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1( 0.983362734, 0.328301698, 0.614836872);
  nxyz2( 1.01367557, 0.344281912, 0.665579319);
  nxyz3( 0.971765339, 0.390600979, 0.676028728);
  nxyz4( 0.941452563, 0.374620765, 0.625286222);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 0.983362734, 0.328301698, 0.614836872);
  nxyz2( 1.03742409, 0.302910298, 0.601357818);
  nxyz3( 1.04909503, 0.317581624, 0.646049917);
  nxyz4( 1.01367557, 0.344281912, 0.665579319);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.0930723, 0.324353248, 0.587483168);
  nxyz2( 1.09125197, 0.341489941, 0.646237373);
  nxyz3( 1.04909503, 0.317581624, 0.646049917);
  nxyz4( 1.03742409, 0.302910298, 0.601357818);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.01367557, 0.344281912, 0.665579319);
  nxyz2( 1.05948687, 0.382861555, 0.678306043);
  nxyz3( 1.01757669, 0.429180622, 0.688755453);
  nxyz4( 0.971765339, 0.390600979, 0.676028728);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.05948687, 0.382861555, 0.678306043);
  nxyz2( 1.01367557, 0.344281912, 0.665579319);
  nxyz3( 1.04909503, 0.317581624, 0.646049917);
  nxyz4( 1.09125197, 0.341489941, 0.646237373);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1( 1.13619518, 0.384683162, 0.63503176);
  nxyz2( 1.10443008, 0.426054776, 0.667100489);
  nxyz3( 1.05948687, 0.382861555, 0.678306043);
  nxyz4( 1.09125197, 0.341489941, 0.646237373);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 1;
  hex8_xyze(1,0) = 0.333333343;
  hex8_xyze(2,0) = 0.666666687;
  hex8_xyze(0,1) = 1;
  hex8_xyze(1,1) = 0.333333343;
  hex8_xyze(2,1) = 0.583333313;
  hex8_xyze(0,2) = 1;
  hex8_xyze(1,2) = 0.416666657;
  hex8_xyze(2,2) = 0.583333313;
  hex8_xyze(0,3) = 1;
  hex8_xyze(1,3) = 0.416666657;
  hex8_xyze(2,3) = 0.666666687;
  hex8_xyze(0,4) = 1.08000004;
  hex8_xyze(1,4) = 0.333333343;
  hex8_xyze(2,4) = 0.666666687;
  hex8_xyze(0,5) = 1.08000004;
  hex8_xyze(1,5) = 0.333333343;
  hex8_xyze(2,5) = 0.583333313;
  hex8_xyze(0,6) = 1.08000004;
  hex8_xyze(1,6) = 0.416666657;
  hex8_xyze(2,6) = 0.583333313;
  hex8_xyze(0,7) = 1.08000004;
  hex8_xyze(1,7) = 0.416666657;
  hex8_xyze(2,7) = 0.666666687;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_hex8_quad4_shadan1()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

#if 1
#if 1
  nxyz1(-0.0864265934,-0.148408353, 0.105090782);
  nxyz2(-0.0732791275,-0.0928697735, 0.130563527);
  nxyz3(-0.098751873775,-0.06420417316, 0.081210988075);
  nxyz4(-0.111899339,-0.119742751, 0.0557382442);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0732791275,-0.0928697735, 0.130563527);
  nxyz2(-0.0601316653,-0.0373311937, 0.156036273);
  nxyz3(-0.0856044069,-0.00866559707, 0.106683731);
  nxyz4(-0.098751873775,-0.06420417316, 0.081210988075);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.098751873775,-0.06420417316, 0.081210988075);
  nxyz2(-0.0856044069,-0.00866559707, 0.106683731);
  nxyz3(-0.111077152, 0.0199999996, 0.0573311932);
  nxyz4(-0.124224618,-0.0355385765, 0.0318584517);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.111899339,-0.119742751, 0.0557382442);
  nxyz2(-0.098751873775,-0.06420417316, 0.081210988075);
  nxyz3(-0.124224618,-0.0355385765, 0.0318584517);
  nxyz4(-0.137372077,-0.0910771564, 0.0063857073);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0864265934,-0.148408353, 0.105090782);
  nxyz2(-0.111899339,-0.119742751, 0.0557382442);
  nxyz3(-0.056360762772,-0.119742751, 0.0270726447);
  nxyz4(-0.0308880173,-0.148408353, 0.0764251798);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.111899339,-0.119742751, 0.0557382442);
  nxyz2(-0.137372077,-0.0910771564, 0.0063857073);
  nxyz3(-0.0818335041,-0.0910771564,-0.0222798903);
  nxyz4(-0.056360762772,-0.119742751, 0.0270726447);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.056360762772,-0.119742751, 0.0270726447);
  nxyz2(-0.0818335041,-0.0910771564,-0.0222798903);
  nxyz3(-0.0262949262,-0.0910771564,-0.0509454869);
  nxyz4(-0.000822183094,-0.119742751,-0.00159295055);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0308880173,-0.148408353, 0.0764251798);
  nxyz2(-0.056360762772,-0.119742751, 0.0270726447);
  nxyz3(-0.000822183094,-0.119742751,-0.00159295055);
  nxyz4( 0.0246505607,-0.148408353, 0.0477595851);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0864265934,-0.148408353, 0.105090782);
  nxyz2(-0.0308880173,-0.148408353, 0.0764251798);
  nxyz3(-0.01774055212,-0.0928697735, 0.10189792465);
  nxyz4(-0.0732791275,-0.0928697735, 0.130563527);

  intersection.AddCutSide( 9, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0308880173,-0.148408353, 0.0764251798);
  nxyz2( 0.0246505607,-0.148408353, 0.0477595851);
  nxyz3( 0.0377980247,-0.0928697735, 0.0732323304);
  nxyz4(-0.01774055212,-0.0928697735, 0.10189792465);

  intersection.AddCutSide( 10, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.01774055212,-0.0928697735, 0.10189792465);
  nxyz2( 0.0377980247,-0.0928697735, 0.0732323304);
  nxyz3( 0.0509454869,-0.0373311937, 0.0987050757);
  nxyz4(-0.00459308969,-0.0373311937, 0.12737067);

  intersection.AddCutSide( 11, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0732791275,-0.0928697735, 0.130563527);
  nxyz2(-0.01774055212,-0.0928697735, 0.10189792465);
  nxyz3(-0.00459308969,-0.0373311937, 0.12737067);
  nxyz4(-0.0601316653,-0.0373311937, 0.156036273);

  intersection.AddCutSide( 12, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.137372077,-0.0910771564, 0.0063857073);
  nxyz2(-0.124224618,-0.0355385765, 0.0318584517);
  nxyz3(-0.14969736925,-0.006872979,-0.017494084545);
  nxyz4(-0.162844822,-0.0624115579,-0.0429668278);

  intersection.AddCutSide( 13, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.124224618,-0.0355385765, 0.0318584517);
  nxyz2(-0.111077152, 0.0199999996, 0.0573311932);
  nxyz3(-0.136549905, 0.048665598, 0.00797865726);
  nxyz4(-0.14969736925,-0.006872979,-0.017494084545);

  intersection.AddCutSide( 14, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.137372077,-0.0910771564, 0.0063857073);
  nxyz2(-0.162844822,-0.0624115579,-0.0429668278);
  nxyz3(-0.107306245075,-0.0624115579,-0.071632426575);
  nxyz4(-0.0818335041,-0.0910771564,-0.0222798903);

  intersection.AddCutSide( 15, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0601316653,-0.0373311937, 0.156036273);
  nxyz2(-0.0469842032, 0.0182073824, 0.181509018);
  nxyz3(-0.0724569467, 0.046872977115, 0.13215647935);
  nxyz4(-0.0856044069,-0.00866559707, 0.106683731);

  intersection.AddCutSide( 16, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0856044069,-0.00866559707, 0.106683731);
  nxyz2(-0.0724569467, 0.046872977115, 0.13215647935);
  nxyz3(-0.0979296938, 0.0755385756, 0.0828039348);
  nxyz4(-0.111077152, 0.0199999996, 0.0573311932);

  intersection.AddCutSide( 17, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.00459308969,-0.0373311937, 0.12737067);
  nxyz2( 0.0509454869,-0.0373311937, 0.0987050757);
  nxyz3( 0.064092949, 0.0182073824, 0.124177814);
  nxyz4( 0.00855437268, 0.0182073824, 0.152843417575);

  intersection.AddCutSide( 18, nids, quad4_xyze, DRT::Element::quad4 );

#else

  nxyz1(-0.0864266,-0.148408,0.105091);
  nxyz2(-0.0732791,-0.0928698,0.130564);
  nxyz3(-0.0987519,-0.0642042,0.081211);
  nxyz4(-0.111899,-0.119743,0.0557382);

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0732791,-0.0928698,0.130564);
  nxyz2(-0.0601317,-0.0373312,0.156036);
  nxyz3(-0.0856044,-0.0086656,0.106684);
  nxyz4(-0.0987519,-0.0642042,0.081211);

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0987519,-0.0642042,0.081211);
  nxyz2(-0.0856044,-0.0086656,0.106684);
  nxyz3(-0.111077,0.02,0.0573312);
  nxyz4(-0.124225,-0.0355386,0.0318585);

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.111899,-0.119743,0.0557382);
  nxyz2(-0.0987519,-0.0642042,0.081211);
  nxyz3(-0.124225,-0.0355386,0.0318585);
  nxyz4(-0.137372,-0.0910772,0.00638571);

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0864266,-0.148408,0.105091);
  nxyz2(-0.111899,-0.119743,0.0557382);
  nxyz3(-0.0563608,-0.119743,0.0270726);
  nxyz4(-0.030888,-0.148408,0.0764252);

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.111899,-0.119743,0.0557382);
  nxyz2(-0.137372,-0.0910772,0.00638571);
  nxyz3(-0.0818335,-0.0910772,-0.0222799);
  nxyz4(-0.0563608,-0.119743,0.0270726);

  intersection.AddCutSide( 6, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0563608,-0.119743,0.0270726);
  nxyz2(-0.0818335,-0.0910772,-0.0222799);
  nxyz3(-0.0262949,-0.0910772,-0.0509455);
  nxyz4(-0.000822183,-0.119743,-0.00159295);

  intersection.AddCutSide( 7, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.030888,-0.148408,0.0764252);
  nxyz2(-0.0563608,-0.119743,0.0270726);
  nxyz3(-0.000822183,-0.119743,-0.00159295);
  nxyz4(0.0246506,-0.148408,0.0477596);

  intersection.AddCutSide( 8, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0864266,-0.148408,0.105091);
  nxyz2(-0.030888,-0.148408,0.0764252);
  nxyz3(-0.0177406,-0.0928698,0.101898);
  nxyz4(-0.0732791,-0.0928698,0.130564);

  intersection.AddCutSide( 9, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.030888,-0.148408,0.0764252);
  nxyz2(0.0246506,-0.148408,0.0477596);
  nxyz3(0.037798,-0.0928698,0.0732323);
  nxyz4(-0.0177406,-0.0928698,0.101898);

  intersection.AddCutSide( 10, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0177406,-0.0928698,0.101898);
  nxyz2(0.037798,-0.0928698,0.0732323);
  nxyz3(0.0509455,-0.0373312,0.0987051);
  nxyz4(-0.00459309,-0.0373312,0.127371);

  intersection.AddCutSide( 11, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0732791,-0.0928698,0.130564);
  nxyz2(-0.0177406,-0.0928698,0.101898);
  nxyz3(-0.00459309,-0.0373312,0.127371);
  nxyz4(-0.0601317,-0.0373312,0.156036);

  intersection.AddCutSide( 12, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.137372,-0.0910772,0.00638571);
  nxyz2(-0.124225,-0.0355386,0.0318585);
  nxyz3(-0.149697,-0.00687298,-0.0174941);
  nxyz4(-0.162845,-0.0624116,-0.0429668);

  intersection.AddCutSide( 13, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.124225,-0.0355386,0.0318585);
  nxyz2(-0.111077,0.02,0.0573312);
  nxyz3(-0.13655,0.0486656,0.00797866);
  nxyz4(-0.149697,-0.00687298,-0.0174941);

  intersection.AddCutSide( 14, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.137372,-0.0910772,0.00638571);
  nxyz2(-0.162845,-0.0624116,-0.0429668);
  nxyz3(-0.107306,-0.0624116,-0.0716324);
  nxyz4(-0.0818335,-0.0910772,-0.0222799);

  intersection.AddCutSide( 15, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0601317,-0.0373312,0.156036);
  nxyz2(-0.0469842,0.0182074,0.181509);
  nxyz3(-0.0724569,0.046873,0.132156);
  nxyz4(-0.0856044,-0.0086656,0.106684);

  intersection.AddCutSide( 16, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.0856044,-0.0086656,0.106684);
  nxyz2(-0.0724569,0.046873,0.132156);
  nxyz3(-0.0979297,0.0755386,0.0828039);
  nxyz4(-0.111077,0.02,0.0573312);

  intersection.AddCutSide( 17, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(-0.00459309,-0.0373312,0.127371);
  nxyz2(0.0509455,-0.0373312,0.0987051);
  nxyz3(0.0640929,0.0182074,0.124178);
  nxyz4(0.00855437,0.0182074,0.152843);

  intersection.AddCutSide( 18, nids, quad4_xyze, DRT::Element::quad4 );

#endif
#endif

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = -0.125;
  hex8_xyze(1,0) = -0.125;
  hex8_xyze(2,0) = 0.125;
  hex8_xyze(0,1) = -0.125;
  hex8_xyze(1,1) = -0.125;
  hex8_xyze(2,1) = 0;
  hex8_xyze(0,2) = -0.125;
  hex8_xyze(1,2) = 0;
  hex8_xyze(2,2) = 0;
  hex8_xyze(0,3) = -0.125;
  hex8_xyze(1,3) = 0;
  hex8_xyze(2,3) = 0.125;
  hex8_xyze(0,4) = 0;
  hex8_xyze(1,4) = -0.125;
  hex8_xyze(2,4) = 0.125;
  hex8_xyze(0,5) = 0;
  hex8_xyze(1,5) = -0.125;
  hex8_xyze(2,5) = 0;
  hex8_xyze(0,6) = 0;
  hex8_xyze(1,6) = 0;
  hex8_xyze(2,6) = 0;
  hex8_xyze(0,7) = 0;
  hex8_xyze(1,7) = 0;
  hex8_xyze(2,7) = 0.125;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

//INFO: This Cut test will not work with DD in LOCAL, as the hex8 does not fullfill the restrictions of this method!
void test_hex8_quad4_mesh_many()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  nids.clear();
  nids.push_back(0);
  quad4_xyze(0,0) = 0.666276;
  quad4_xyze(1,0) = 0.828602;
  quad4_xyze(2,0) = 0.02525;
  nids.push_back(1);
  quad4_xyze(0,1) = 0.673325;
  quad4_xyze(1,1) = 0.847231;
  quad4_xyze(2,1) = 0.02525;
  nids.push_back(2);
  quad4_xyze(0,2) = 0.673325;
  quad4_xyze(1,2) = 0.847231;
  quad4_xyze(2,2) = -0.02525;
  nids.push_back(3);
  quad4_xyze(0,3) = 0.666276;
  quad4_xyze(1,3) = 0.828602;
  quad4_xyze(2,3) = -0.02525;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(4);
  quad4_xyze(0,0) = 0.680311;
  quad4_xyze(1,0) = 0.86595;
  quad4_xyze(2,0) = 0.02525;
  nids.push_back(5);
  quad4_xyze(0,1) = 0.68731;
  quad4_xyze(1,1) = 0.884695;
  quad4_xyze(2,1) = 0.02525;
  nids.push_back(6);
  quad4_xyze(0,2) = 0.68731;
  quad4_xyze(1,2) = 0.884695;
  quad4_xyze(2,2) = -0.02525;
  nids.push_back(7);
  quad4_xyze(0,3) = 0.680311;
  quad4_xyze(1,3) = 0.86595;
  quad4_xyze(2,3) = -0.02525;

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(1);
  quad4_xyze(0,0) = 0.673325;
  quad4_xyze(1,0) = 0.847231;
  quad4_xyze(2,0) = 0.02525;
  nids.push_back(4);
  quad4_xyze(0,1) = 0.680311;
  quad4_xyze(1,1) = 0.86595;
  quad4_xyze(2,1) = 0.02525;
  nids.push_back(7);
  quad4_xyze(0,2) = 0.680311;
  quad4_xyze(1,2) = 0.86595;
  quad4_xyze(2,2) = -0.02525;
  nids.push_back(2);
  quad4_xyze(0,3) = 0.673325;
  quad4_xyze(1,3) = 0.847231;
  quad4_xyze(2,3) = -0.02525;

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(6);
  quad4_xyze(0,0) = 0.68731;
  quad4_xyze(1,0) = 0.884695;
  quad4_xyze(2,0) = -0.02525;
  nids.push_back(5);
  quad4_xyze(0,1) = 0.68731;
  quad4_xyze(1,1) = 0.884695;
  quad4_xyze(2,1) = 0.02525;
  nids.push_back(8);
  quad4_xyze(0,2) = 0.702923;
  quad4_xyze(1,2) = 0.878831;
  quad4_xyze(2,2) = 0.02525;
  nids.push_back(9);
  quad4_xyze(0,3) = 0.702923;
  quad4_xyze(1,3) = 0.878831;
  quad4_xyze(2,3) = -0.02525;

  intersection.AddCutSide( 4, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(9);
  quad4_xyze(0,0) = 0.702923;
  quad4_xyze(1,0) = 0.878831;
  quad4_xyze(2,0) = -0.02525;
  nids.push_back(8);
  quad4_xyze(0,1) = 0.702923;
  quad4_xyze(1,1) = 0.878831;
  quad4_xyze(2,1) = 0.02525;
  nids.push_back(10);
  quad4_xyze(0,2) = 0.718519;
  quad4_xyze(1,2) = 0.87295;
  quad4_xyze(2,2) = 0.02525;
  nids.push_back(11);
  quad4_xyze(0,3) = 0.718519;
  quad4_xyze(1,3) = 0.87295;
  quad4_xyze(2,3) = -0.02525;

  intersection.AddCutSide( 5, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.707126021;
  hex8_xyze(1,0) = 0.885084629;
  hex8_xyze(2,0) = -0.0250000004;
  hex8_xyze(0,1) = 0.669295728;
  hex8_xyze(1,1) = 0.88072747;
  hex8_xyze(2,1) = -0.0250000004;
  hex8_xyze(0,2) = 0.669837534;
  hex8_xyze(1,2) = 0.831363618;
  hex8_xyze(2,2) = -0.0250000004;
  hex8_xyze(0,3) = 0.713419497;
  hex8_xyze(1,3) = 0.846297681;
  hex8_xyze(2,3) = -0.0250000004;
  hex8_xyze(0,4) = 0.707126021;
  hex8_xyze(1,4) = 0.885084629;
  hex8_xyze(2,4) = 0.0250000004;
  hex8_xyze(0,5) = 0.669295728;
  hex8_xyze(1,5) = 0.88072747;
  hex8_xyze(2,5) = 0.0250000004;
  hex8_xyze(0,6) = 0.669837534;
  hex8_xyze(1,6) = 0.831363618;
  hex8_xyze(2,6) = 0.0250000004;
  hex8_xyze(0,7) = 0.713419497;
  hex8_xyze(1,7) = 0.846297681;
  hex8_xyze(2,7) = 0.0250000004;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.CutTest_Cut( true );
}

//INFO: This Cut test will not work with DD in LOCAL, as the hex8 does not fullfill the restrictions of this method!
void test_hex8_quad4_mesh_edgecut()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  nids.clear();
  nids.push_back(0);
  quad4_xyze(0,0) = 0.5399141482156006866;
  quad4_xyze(1,0) = 0.77952300801931839747;
  quad4_xyze(2,0) = 0.025250000899999999748;
  nids.push_back(1);
  quad4_xyze(0,1) = 0.53991414821296357385;
  quad4_xyze(1,1) = 0.77952300803036944643;
  quad4_xyze(2,1) = -0.025250000899999999748;
  nids.push_back(2);
  quad4_xyze(0,2) = 0.5400635612166613253;
  quad4_xyze(1,2) = 0.79951197365630055636;
  quad4_xyze(2,2) = -0.025250000899999999748;
  nids.push_back(3);
  quad4_xyze(0,3) = 0.54006356122364373995;
  quad4_xyze(1,3) = 0.79951197365214676793;
  quad4_xyze(2,3) = 0.025250000899999999748;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(4);
  quad4_xyze(0,0) = 0.53955076912075083939;
  quad4_xyze(1,0) = 0.73952115635249415782;
  quad4_xyze(2,0) = 0.025250000899999999748;
  nids.push_back(5);
  quad4_xyze(0,1) = 0.53955076912314114956;
  quad4_xyze(1,1) = 0.73952115637886628452;
  quad4_xyze(2,1) = -0.025250000899999999748;
  nids.push_back(6);
  quad4_xyze(0,2) = 0.53974580830286555955;
  quad4_xyze(1,2) = 0.75952575976383629452;
  quad4_xyze(2,2) = -0.025250000899999999748;
  nids.push_back(7);
  quad4_xyze(0,3) = 0.53974580831400054137;
  quad4_xyze(1,3) = 0.75952575975287395238;
  quad4_xyze(2,3) = 0.025250000899999999748;

  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  nids.clear();
  nids.push_back(7);
  quad4_xyze(0,0) = 0.53974580831400054137;
  quad4_xyze(1,0) = 0.75952575975287395238;
  quad4_xyze(2,0) = 0.025250000899999999748;
  nids.push_back(6);
  quad4_xyze(0,1) = 0.53974580830286555955;
  quad4_xyze(1,1) = 0.75952575976383629452;
  quad4_xyze(2,1) = -0.025250000899999999748;
  nids.push_back(1);
  quad4_xyze(0,2) = 0.53991414821296357385;
  quad4_xyze(1,2) = 0.77952300803036944643;
  quad4_xyze(2,2) = -0.025250000899999999748;
  nids.push_back(0);
  quad4_xyze(0,3) = 0.5399141482156006866;
  quad4_xyze(1,3) = 0.77952300801931839747;
  quad4_xyze(2,3) = 0.025250000899999999748;

  intersection.AddCutSide( 3, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.53135645399999997807;
  hex8_xyze(1,0) = 0.75096273400000002063;
  hex8_xyze(2,0) = -0.025000000399999999789;
  hex8_xyze(0,1) = 0.56228893999999995934;
  hex8_xyze(1,1) = 0.74424672099999999997;
  hex8_xyze(2,1) = -0.025000000399999999789;
  hex8_xyze(0,2) = 0.57124853099999994832;
  hex8_xyze(1,2) = 0.7783957719999999858;
  hex8_xyze(2,2) = -0.025000000399999999789;
  hex8_xyze(0,3) = 0.53994804600000001482;
  hex8_xyze(1,3) = 0.78398263499999998327;
  hex8_xyze(2,3) = -0.025000000399999999789;
  hex8_xyze(0,4) = 0.53135645399999997807;
  hex8_xyze(1,4) = 0.75096273400000002063;
  hex8_xyze(2,4) = 0.025000000399999999789;
  hex8_xyze(0,5) = 0.56228893999999995934;
  hex8_xyze(1,5) = 0.74424672099999999997;
  hex8_xyze(2,5) = 0.025000000399999999789;
  hex8_xyze(0,6) = 0.57124853099999994832;
  hex8_xyze(1,6) = 0.7783957719999999858;
  hex8_xyze(2,6) = 0.025000000399999999789;
  hex8_xyze(0,7) = 0.53994804600000001482;
  hex8_xyze(1,7) = 0.78398263499999998327;
  hex8_xyze(2,7) = 0.025000000399999999789;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.CutTest_Cut( true );
}

//INFO: This Cut test will not work with DD in LOCAL, as the hex8 does not fullfill the restrictions of this method!
void test_hex8_quad4_mesh_edgecut2()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  nids.clear();
  nids.push_back(0);
  quad4_xyze(0,0) = 0.53991414821573446847;
  quad4_xyze(1,0) = 0.7795230080202141254;
  quad4_xyze(2,0) = 0.025250000899999999748;
  nids.push_back(1);
  quad4_xyze(0,1) = 0.53991414821143801639;
  quad4_xyze(1,1) = 0.77952300802897134258;
  quad4_xyze(2,1) = -0.025250000899999999748;
  nids.push_back(2);
  quad4_xyze(0,2) = 0.54006356121681575733;
  quad4_xyze(1,2) = 0.79951197365563708708;
  quad4_xyze(2,2) = -0.025250000899999999748;
  nids.push_back(3);
  quad4_xyze(0,3) = 0.54006356122155441124;
  quad4_xyze(1,3) = 0.79951197365235482373;
  quad4_xyze(2,3) = 0.025250000899999999748;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.53135645399999997807;
  hex8_xyze(1,0) = 0.75096273400000002063;
  hex8_xyze(2,0) = -0.025000000399999999789;
  hex8_xyze(0,1) = 0.53994804600000001482;
  hex8_xyze(1,1) = 0.78398263499999998327;
  hex8_xyze(2,1) = -0.025000000399999999789;
  hex8_xyze(0,2) = 0.50934117999999994897;
  hex8_xyze(1,2) = 0.78892862799999996515;
  hex8_xyze(2,2) = -0.025000000399999999789;
  hex8_xyze(0,3) = 0.50265735400000000066;
  hex8_xyze(1,3) = 0.75846225000000000449;
  hex8_xyze(2,3) = -0.025000000399999999789;
  hex8_xyze(0,4) = 0.53135645399999997807;
  hex8_xyze(1,4) = 0.75096273400000002063;
  hex8_xyze(2,4) = 0.025000000399999999789;
  hex8_xyze(0,5) = 0.53994804600000001482;
  hex8_xyze(1,5) = 0.78398263499999998327;
  hex8_xyze(2,5) = 0.025000000399999999789;
  hex8_xyze(0,6) = 0.50934117999999994897;
  hex8_xyze(1,6) = 0.78892862799999996515;
  hex8_xyze(2,6) = 0.025000000399999999789;
  hex8_xyze(0,7) = 0.50265735400000000066;
  hex8_xyze(1,7) = 0.75846225000000000449;
  hex8_xyze(2,7) = 0.025000000399999999789;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.CutTest_Cut( true );
}

void test_hex8_quad4_mesh_inner()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids;

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  nids.clear();
  nids.push_back(0);
  quad4_xyze(0,0) = 0.162023;
  quad4_xyze(1,0) = -0.0373312;
  quad4_xyze(2,0) = 0.0413739;
  nids.push_back(1);
  quad4_xyze(0,1) = 0.13655;
  quad4_xyze(1,1) = -0.0086656;
  quad4_xyze(2,1) = -0.00797866;
  nids.push_back(2);
  quad4_xyze(0,2) = 0.149697;
  quad4_xyze(1,2) = 0.046873;
  quad4_xyze(2,2) = 0.0174941;
  nids.push_back(3);
  quad4_xyze(0,3) = 0.17517;
  quad4_xyze(1,3) = 0.0182074;
  quad4_xyze(2,3) = 0.0668466;

  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix hex8_xyze( 3, 8 );

  hex8_xyze(0,0) = 0.0833333;
  hex8_xyze(1,0) = -0.0833333;
  hex8_xyze(2,0) = 0.0833333;
  hex8_xyze(0,1) = 0.0833333;
  hex8_xyze(1,1) = -0.0833333;
  hex8_xyze(2,1) = 1.38778e-17;
  hex8_xyze(0,2) = 0.0833333;
  hex8_xyze(1,2) = -3.46945e-18;
  hex8_xyze(2,2) = 1.38778e-17;
  hex8_xyze(0,3) = 0.0833333;
  hex8_xyze(1,3) = -6.93889e-18;
  hex8_xyze(2,3) = 0.0833333;
  hex8_xyze(0,4) = 0.166667;
  hex8_xyze(1,4) = -0.0833333;
  hex8_xyze(2,4) = 0.0833333;
  hex8_xyze(0,5) = 0.166667;
  hex8_xyze(1,5) = -0.0833333;
  hex8_xyze(2,5) = 1.38778e-17;
  hex8_xyze(0,6) = 0.166667;
  hex8_xyze(1,6) = -1.73472e-18;
  hex8_xyze(2,6) = 1.38778e-17;
  hex8_xyze(0,7) = 0.166667;
  hex8_xyze(1,7) = 0;
  hex8_xyze(2,7) = 0.0833333;

  nids.clear();
  for ( int i=0; i<8; ++i )
    nids.push_back( i );

  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );

  intersection.CutTest_Cut( true );
}

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

  intersection.CutTest_Cut( true );
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

  intersection.CutTest_Cut( true );
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

  intersection.CutTest_Cut( true ,INPAR::CUT::VCellGaussPts_Tessellation,INPAR::CUT::BCellGaussPts_Tessellation,true,true,true);
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

  intersection.CutTest_Cut( true );
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

  intersection.CutTest_Cut( true );
}

void test_tet4_quad4_double()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 4 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 4 );

  std::map<std::string, int> nodeids;

  nxyz1(0,0,0.9);
  nxyz2(1,0,0.9);
  nxyz3(1,0.5,0.9);
  nxyz4(0,0.5,0.9);

  //GEO::CUT::SideHandle * s1 =
  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::quad4 );

  nxyz1(0,0.5,0.9);
  nxyz2(1,0.5,0.9);
  nxyz3(1,1,0.9);
  nxyz4(0,1,0.9);

  //GEO::CUT::SideHandle * s2 =
  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::quad4 );

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 0;
  tet4_xyze(1,0) = 0;
  tet4_xyze(2,0) = 0;
  tet4_xyze(0,1) = 1;
  tet4_xyze(1,1) = 0.5;
  tet4_xyze(2,1) = 0;
  tet4_xyze(0,2) = 0;
  tet4_xyze(1,2) = 1;
  tet4_xyze(2,2) = 0;
  tet4_xyze(0,3) = 0.5;
  tet4_xyze(1,3) = 0.52;
  tet4_xyze(2,3) = 1;

  nids.clear();
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  //GEO::CUT::ElementHandle * e =
  intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
}

void test_tet4_tri3_double()
{
  GEO::CUT::MeshIntersection intersection;

  std::vector<int> nids( 3 );

  Epetra_SerialDenseMatrix quad4_xyze( 3, 3 );

  std::map<std::string, int> nodeids;

  nxyz1(0,0,0);
  nxyz2(1,0.5,0.9);
  nxyz3(0,0.5,0.9);

  //GEO::CUT::SideHandle * s1 =
  intersection.AddCutSide( 1, nids, quad4_xyze, DRT::Element::tri3 );

  nxyz1(0,0.5,0.9);
  nxyz2(1,0.5,0.9);
  nxyz3(0,1,0);

  //GEO::CUT::SideHandle * s2 =
  intersection.AddCutSide( 2, nids, quad4_xyze, DRT::Element::tri3 );

  Epetra_SerialDenseMatrix tet4_xyze( 3, 4 );

  tet4_xyze(0,0) = 0;
  tet4_xyze(1,0) = 0;
  tet4_xyze(2,0) = 0;
  tet4_xyze(0,1) = 1;
  tet4_xyze(1,1) = 0.5;
  tet4_xyze(2,1) = 0;
  tet4_xyze(0,2) = 0;
  tet4_xyze(1,2) = 1;
  tet4_xyze(2,2) = 0;
  tet4_xyze(0,3) = 0.5;
  tet4_xyze(1,3) = 0.52;
  tet4_xyze(2,3) = 1;

  nids.clear();
  for ( int i=0; i<4; ++i )
    nids.push_back( i );

  //GEO::CUT::ElementHandle * e =
  intersection.AddElement( 1, nids, tet4_xyze, DRT::Element::tet4 );

  intersection.Status();
  intersection.CutTest_Cut( true );
  intersection.Status();
  std::cout << __LINE__ << std::endl;
  std::vector<double> tessVol,momFitVol,dirDivVol;

  GEO::CUT::Mesh mesh = intersection.NormalMesh();
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell> > & other_cells = mesh.VolumeCells();
  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
        i!=other_cells.end();
        ++i )
  {
    GEO::CUT::VolumeCell * vc = &**i;
    tessVol.push_back(vc->Volume());
  }

  std::cout<<"the volumes predicted by\n tessellation ";
  for(unsigned i=0;i<tessVol.size();i++)
  {
    std::cout<<tessVol[i]<<"\n";
  }
}
