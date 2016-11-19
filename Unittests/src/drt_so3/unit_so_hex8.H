/*!---------------------------------------------------------------------*
\file unit_so_hex8.H

\brief Unittests for the So_hex8 class

\level 3

\maintainer Michael Gee
*-----------------------------------------------------------------------*/

#ifndef UNIT_SO_HEX8_H_
#define UNIT_SO_HEX8_H_

#include <cxxtest/TestSuite.h>
#include <Epetra_SerialComm.h>
#include "src/drt_so3/so_hex8.H"

namespace DRT { namespace ELEMENTS { class So_hex8_TestSuite; } }

class DRT::ELEMENTS::So_hex8_TestSuite : public CxxTest::TestSuite
{
 private:
  Teuchos::RCP<Discretization> testdis_;
  Teuchos::RCP<So_hex8> testele_;
  Teuchos::RCP<So_hex8> copytestele_;

 public:
  void setUp()
  {
    // create a discretization, that creates node to element pointers and keeps the nodes alive
    testdis_ = Teuchos::rcp(new Discretization("dummy",Teuchos::rcp(new Epetra_SerialComm)));

    // create 8 nodes
    const int nodeids[8] = {0,1,2,3,4,5,6,7};
    double coords[24] = {-0.1,-0.2,-0.5, 1.25,0.23,0.66, 1.20,0.99,0.5, -0.11,1.20,0.66, -0.10,-0.2,1.9, 1.00,0.00,1.90, 1.20,0.99,1.50, -0.11,-0.20,1.66};
    for (int lid = 0; lid < 8; ++lid)
      testdis_->AddNode(Teuchos::rcp(new Node(lid,&coords[3*lid],0)));

    // create 1 element
    testele_ = Teuchos::rcp(new So_hex8(0, 0));
    testele_->SetNodeIds(8, nodeids);
    testdis_->AddElement(testele_);
    testdis_->FillComplete(false, false, false);

    copytestele_ = Teuchos::rcp(new So_hex8(*testele_));
  }

  void tearDown()
  {
    copytestele_ = Teuchos::null;
    testele_ = Teuchos::null;
    testdis_ = Teuchos::null;
  }

  void test_NumDofPerElement()
  {
    TS_ASSERT_EQUALS(testele_->NumDofPerElement(), 0);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerElement(), 0);
  }

  void test_NumDofPerNode()
  {
    double pd[3] = {1,2,3};
    DRT::Node node_dummy(0,pd,false);
    TS_ASSERT_EQUALS(testele_->NumDofPerNode(node_dummy), 3);
    TS_ASSERT_EQUALS(copytestele_->NumDofPerNode(node_dummy), 3);
  }

  void test_Degree()
  {
    TS_ASSERT_EQUALS(testele_->Degree(), 1);
    TS_ASSERT_EQUALS(copytestele_->Degree(), 1);
  }

  void test_NumVolume()
  {
    TS_ASSERT_EQUALS(testele_->NumVolume(), 1);
    TS_ASSERT_EQUALS(copytestele_->NumVolume(), 1);
  }

  void test_NumSurface()
  {
    TS_ASSERT_EQUALS(testele_->NumSurface(), 6);
    TS_ASSERT_EQUALS(copytestele_->NumSurface(), 6);
  }

  void test_NumLine()
  {
    TS_ASSERT_EQUALS(testele_->NumLine(), 12);
    TS_ASSERT_EQUALS(copytestele_->NumLine(), 12);
  }

  void test_ElementCenterRefeCoords()
  {
    double midpoint[3] = {0.528750000000000,0.351250000000000,1.035000000000000};
    for (int i = 0; i < 3; ++i)
    {
      TS_ASSERT_DELTA(testele_->ElementCenterRefeCoords()[i],midpoint[i],1e-14);
      TS_ASSERT_DELTA(copytestele_->ElementCenterRefeCoords()[i],midpoint[i],1e-14);
    }
  }

  // this is a protected method and should not be tested, but until more testcases
  // have been acumulated it is left in as an example
  //void DISABLE_test_soh8_expol()
  //{
  //  // assign stresses at gauss points
  //  double gp_stresses[48] =   {1,1.1,1.2,1.3,1.4,1.5,1.6,1.7, 1,1,1,1,1,1,1,1,
  //    100,100,101,102,105,108,107.9,1, 0,0,0,0,0,0,0,0,
  //    0.154,0.464,0.534,0.654,0.652,0.132,0.789,0.561, 1,100,1000,12,56,800,160,165};
  //    // {sigma11_gp1,sigma11_gp2,...,sigma22_gp1,sigma22_gp2, ...}
  //  LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D> stresses(gp_stresses);
  //
  //  // assign extrapolated stress targets at nodes
  //  double node_stresses[48] = {
  //    0.6803847577293368,0.9803847577293369,1.0267949192431121,1.3267949192431121,1.3732050807568876,1.6732050807568877,1.7196152422706634,2.0196152422706630,
  //    1.0000000000000000,1.0000000000000000,1.0000000000000000,0.9999999999999999,0.9999999999999999,1.0000000000000000,0.9999999999999999,1.0000000000000000,
  //    78.3603668816091243,101.4345474031532319,78.8878659563684721,174.0018800057871999,176.6621340436315961,93.5481199942128399,182.8896331183908615,-160.8845474031532206,
  //    0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
  //    -0.6305150942321529,0.9785844986619437,0.1817200395244018,1.1561542236045113,1.4062799604755982,-0.7681542236045114,1.4895150942321531,0.1264155013380563,
  //    239.5767664977294942,-941.5664569454437469,2537.3027868719659637,-747.2028238815928489,-516.3027868719652815,2036.2028238815928489,-903.5767664977294089,589.5664569454437469};
  //
  //  // build Epetra_MultiVector
  //  Epetra_SerialComm SerComm;
  //  Epetra_Map Map(NUMGPT_SOH8,0,SerComm);
  //  Epetra_MultiVector expolstresses(Map,MAT::NUM_STRESS_3D);
  //
  //  // test soh8_expol()
  //  testele_->soh8_expol(stresses,expolstresses);
  //
  //  int i = 0;
  //  for (int sid = 0; sid < MAT::NUM_STRESS_3D; ++sid)
  //  {
  //    for (int lid = 0; lid < NUMGPT_SOH8; ++lid)
  //    {
  //      TS_ASSERT_DELTA((*(expolstresses(sid)))[lid],node_stresses[i],1e-14);
  //      i += 1;
  //    }
  //  }
  //}

}; // class So_hex8_TestSuite


#endif /* UNIT_SO_Hex8_H_ */
