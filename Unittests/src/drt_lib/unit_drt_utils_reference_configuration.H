/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests to check utility functions for the reference configuration

\level 2
*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                              cschmidt 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_DRT_UTILS_REFERENCE_CONFIGURATION_H_
#define UNIT_DRT_UTILS_REFERENCE_CONFIGURATION_H_

/*---------------------------------------------------------------------------*
 | headers                                                  cschmidt 11/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_lib/drt_utils_reference_configuration.H"
#include "src/drt_so3/so_tet4.H"
#include "src/drt_so3/so_hex8.H"
#include "src/drt_contact/contact_element.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                      cschmidt 11/2018 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace UTILS
  {
    class UtilsRefConfig_TestSuite;
  }
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | utils reference configuration test suite                 cschmidt 11/2018 |
 *---------------------------------------------------------------------------*/
class DRT::UTILS::UtilsRefConfig_TestSuite : public CxxTest::TestSuite
{
 private:
  Teuchos::RCP<Discretization> testdis_;

 public:
  void setUp()
  {
    // create a discretization, to store the created elements and nodes
    testdis_ = Teuchos::rcp(new DRT::Discretization("dummy", Teuchos::rcp(new Epetra_SerialComm)));

    // create hex8 element and store it in the test discretization
    const int nodeidshex8[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const double coordshex8[24] = {-0.10, -0.20, -0.50, 1.25, 0.23, 0.66, 1.20, 0.99, 0.50, -0.11,
        1.20, 0.66, -0.10, -0.20, 1.90, 1.00, 0.00, 1.90, 1.20, 0.99, 1.50, -0.11, -0.20, 1.66};
    for (int i = 0; i < 8; ++i)
    {
      testdis_->AddNode(Teuchos::rcp(new DRT::Node(nodeidshex8[i], &coordshex8[3 * i], 0)));
    }
    Teuchos::RCP<DRT::ELEMENTS::So_hex8> testhex8ele =
        Teuchos::rcp(new DRT::ELEMENTS::So_hex8(0, 0));
    testhex8ele->SetNodeIds(8, nodeidshex8);
    testdis_->AddElement(testhex8ele);

    // create corresponding quad4 surface contact element and store it
    Teuchos::RCP<CONTACT::CoElement> testcontactquad4ele =
        Teuchos::rcp(new CONTACT::CoElement(testhex8ele->Id() + 1, testhex8ele->Owner(),
            testhex8ele->Shape(), testhex8ele->NumNode(), testhex8ele->NodeIds(), false, false));
    testdis_->AddElement(testcontactquad4ele);

    // create tet4 element and store it in the test discretization
    const int nodeidstet4[4] = {8, 9, 10, 11};
    const double coordstet4[12] = {2.5, -0.5, 0.0, 1.0, -1.1, 0.1, 1.1, 0.11, 0.15, 1.5, -0.5, 2.0};
    for (int j = 0; j < 4; ++j)
    {
      testdis_->AddNode(Teuchos::rcp(new DRT::Node(nodeidstet4[j], &coordstet4[3 * j], 0)));
    }
    Teuchos::RCP<DRT::ELEMENTS::So_tet4> testtet4ele =
        Teuchos::rcp(new DRT::ELEMENTS::So_tet4(2, 0));
    testtet4ele->SetNodeIds(4, nodeidstet4);
    testdis_->AddElement(testtet4ele);

    // create corresponding tri3 surface contact element and store it
    Teuchos::RCP<CONTACT::CoElement> testcontacttri3ele =
        Teuchos::rcp(new CONTACT::CoElement(testtet4ele->Id() + 1, testtet4ele->Owner(),
            testtet4ele->Shape(), testtet4ele->NumNode(), testtet4ele->NodeIds(), false, false));
    testdis_->AddElement(testcontacttri3ele);
    testdis_->FillComplete(false, false, false);
  }

  void tearDown() { testdis_ = Teuchos::null; }

  void test_LocalToGlobalPositionAtXiRefConfig()
  {
    // get hex8 element and test it
    const DRT::Element* hex8ele = testdis_->gElement(0);
    LINALG::Matrix<3, 1> xicenterhex8ele(true);
    LINALG::Matrix<3, 1> hex8elecoords(true);
    const double hex8refsolution[3] = {423.0 / 800.0, 281.0 / 800.0, 207.0 / 200.0};
    LocalToGlobalPositionAtXiRefConfig<3, DRT::Element::hex8>(
        hex8ele, xicenterhex8ele, hex8elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(hex8elecoords(i, 0), hex8refsolution[i], 1e-14);

    // get quad4 element and test it
    const DRT::Element* quad4ele = testdis_->gElement(1);
    LINALG::Matrix<2, 1> xicenterquad4ele(true);
    LINALG::Matrix<3, 1> quad4elecoords(true);
    const double quad4refsolution[3] = {14.0 / 25.0, 111.0 / 200.0, 33.0 / 100.0};
    LocalToGlobalPositionAtXiRefConfig<3, DRT::Element::quad4>(
        quad4ele, xicenterquad4ele, quad4elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(quad4elecoords(i, 0), quad4refsolution[i], 1e-14);

    // get tet4 element stuff and test it
    const DRT::Element* tet4ele = testdis_->gElement(2);
    LINALG::Matrix<3, 1> xicentertet4ele(true);
    LINALG::Matrix<3, 1> tet4elecoords(true);
    xicentertet4ele.PutScalar(1.0 / 4.0);
    const double tet4refsolution[3] = {61.0 / 40.0, -199.0 / 400.0, 9.0 / 16.0};
    LocalToGlobalPositionAtXiRefConfig<3, DRT::Element::tet4>(
        tet4ele, xicentertet4ele, tet4elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(tet4elecoords(i, 0), tet4refsolution[i], 1e-14);

    // get tri3 element and test it
    const DRT::Element* tri3ele = testdis_->gElement(3);
    LINALG::Matrix<2, 1> xicentertri3ele(true);
    LINALG::Matrix<3, 1> tri3elecoords(true);
    xicentertri3ele.PutScalar(1.0 / 3.0);
    const double tri3refsolution[3] = {23.0 / 15.0, -149.0 / 300.0, 1.0 / 12.0};
    LocalToGlobalPositionAtXiRefConfig<3, DRT::Element::tri3>(
        tri3ele, xicentertri3ele, tri3elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(tri3elecoords(i, 0), tri3refsolution[i], 1e-14);
  }

  void test_ComputeUnitNormalAtXiRefConfig()
  {
    // get quad4 element and test it
    const DRT::Element* quad4ele = testdis_->gElement(1);
    LINALG::Matrix<2, 1> xicenterquad4ele(true);
    LINALG::Matrix<3, 1> quad4elecoords(true);
    const double quad4refsolution[3] = {-0.29138926578643, -0.40854577471087, 0.86497551742829};
    ComputeUnitNormalAtXiRefConfig<DRT::Element::quad4>(quad4ele, xicenterquad4ele, quad4elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(quad4elecoords(i, 0), quad4refsolution[i], 1e-14);

    // get tri3 element and test it
    const DRT::Element* tri3ele = testdis_->gElement(3);
    LINALG::Matrix<2, 1> xicentertri3ele(true);
    LINALG::Matrix<3, 1> tri3elecoords(true);
    xicentertri3ele.PutScalar(1.0 / 3.0);
    const double tri3refsolution[3] = {-0.085623542490578, 0.048198682858935, -0.995161040205065};
    ComputeUnitNormalAtXiRefConfig<DRT::Element::tri3>(tri3ele, xicentertri3ele, tri3elecoords);
    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(tri3elecoords(i, 0), tri3refsolution[i], 1e-14);
  }
};

#endif  // UNIT_DRT_UTILS_REFERENCE_CONFIGURATION_H_
