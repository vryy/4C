
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
  intersection.Cut();
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
  intersection.Cut();
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
  intersection.Cut();
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
  intersection.Cut();
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
  intersection.Cut();
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
  intersection.Cut();
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
  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
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

  intersection.Cut();
}
