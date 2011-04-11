
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

#define nxyz( x, y, z ) { if ( nodeids.count( #x#y#z )==0 ) { nodeids[#x#y#z] = nodeids.size(); }}
#define nxyz1( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 0 )=x; quad4_xyze( 1, 0 )=y; quad4_xyze( 2, 0 )=z; nids[0] = nodeids[#x#y#z];
#define nxyz2( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 1 )=x; quad4_xyze( 1, 1 )=y; quad4_xyze( 2, 1 )=z; nids[1] = nodeids[#x#y#z];
#define nxyz3( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 2 )=x; quad4_xyze( 1, 2 )=y; quad4_xyze( 2, 2 )=z; nids[2] = nodeids[#x#y#z];
#define nxyz4( x, y, z ) nxyz( x, y, z ); quad4_xyze( 0, 3 )=x; quad4_xyze( 1, 3 )=y; quad4_xyze( 2, 3 )=z; nids[3] = nodeids[#x#y#z];


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
  intersection.Cut( true );
  intersection.Status();
}

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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
  intersection.Status();
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
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
  intersection.Cut( true );
  intersection.Status();
}

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

  intersection.Cut( true );
}

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

  intersection.Cut( true );
}

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

  intersection.Cut( true );
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

  intersection.Cut( true );
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

  intersection.Cut( true );
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

  intersection.Cut( true );
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

  intersection.Cut( true );
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

  intersection.Cut( true );
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

  intersection.Cut( true );
}
