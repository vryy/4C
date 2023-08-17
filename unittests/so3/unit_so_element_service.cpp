/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the methods in so_element_service

\level 3

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_lib_element.H"
#include "baci_so3_element_service.H"

#include <vector>

namespace
{
  TEST(ElementServiceTest, TestProjectNodalQuantityToXiHex8)
  {
    CORE::LINALG::Matrix<3, 1> xi(true);
    xi(0) = 0.01, xi(1) = 0.25, xi(2) = 0.115;
    std::vector<double> nodal_quantity = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7};
    std::vector<double> ref_val{1.397875};
    auto test_val = DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::hex8>(xi, nodal_quantity);
    for (std::size_t i = 0; i < ref_val.size(); ++i) EXPECT_NEAR(test_val[i], ref_val[i], 1.0e-10);
  }

  TEST(ElementServiceTest, TestProjectNodalQuantityToXiHex27)
  {
    CORE::LINALG::Matrix<3, 1> xi(true);
    xi(0) = 0.01, xi(1) = 0.25, xi(2) = 0.115;
    std::vector<double> nodal_quantity = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6};
    std::vector<double> ref_val = {3.623649611383};
    auto test_val =
        DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::hex27>(xi, nodal_quantity);
    for (std::size_t i = 0; i < ref_val.size(); ++i) EXPECT_NEAR(test_val[i], ref_val[i], 1.0e-10);
  }

  TEST(ElementServiceTest, TestProjectNodalQuantityToXiTet4)
  {
    CORE::LINALG::Matrix<3, 1> xi(true);
    xi(0) = 0.01, xi(1) = 0.25, xi(2) = 0.115;
    std::vector<double> nodal_quantity = {1.0, 2.0, 1.1, 2.1, 1.2, 2.2, 1.3, 2.3};
    std::vector<double> ref_val = {1.0855, 2.0855};
    auto test_val = DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::tet4>(xi, nodal_quantity);
    for (std::size_t i = 0; i < ref_val.size(); ++i) EXPECT_NEAR(test_val[i], ref_val[i], 1.0e-10);
  }

  TEST(ElementServiceTest, TestProjectNodalQuantityToXiTet10)
  {
    CORE::LINALG::Matrix<3, 1> xi(true);
    xi(0) = 0.01, xi(1) = 0.25, xi(2) = 0.115;
    std::vector<double> nodal_quantity = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
    std::vector<double> ref_val = {1.645885};
    auto test_val =
        DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::tet10>(xi, nodal_quantity);
    for (std::size_t i = 0; i < ref_val.size(); ++i) EXPECT_NEAR(test_val[i], ref_val[i], 1.0e-10);
  }

  TEST(ElementServiceTest, TestProjectNodalQuantityToXiWedge6)
  {
    CORE::LINALG::Matrix<3, 1> xi(true);
    xi(0) = 0.01, xi(1) = 0.25, xi(2) = 0.115;
    std::vector<double> nodal_quantity = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
    std::vector<double> ref_val = {1.34025};
    auto test_val =
        DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::wedge6>(xi, nodal_quantity);
    for (std::size_t i = 0; i < ref_val.size(); ++i) EXPECT_NEAR(test_val[i], ref_val[i], 1.0e-10);
  }

  TEST(ElementServiceTest, TestGaussPointProjectionMatrixHex8)
  {
    constexpr DRT::Element::DiscretizationType distype = DRT::Element::DiscretizationType::hex8;
    constexpr int nsd = 3;

    CORE::DRT::UTILS::IntPointsAndWeights<nsd> intpoints(
        DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<CORE::DRT::UTILS::CollectedGaussPoints> gp =
        Teuchos::rcp(new CORE::DRT::UTILS::CollectedGaussPoints);

    std::array<double, nsd> xi{};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < nsd; ++d) xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    // save default integration rule
    CORE::DRT::UTILS::GaussIntegration integration(gp);
    LINALG::SerialDenseMatrix m =
        CORE::DRT::UTILS::EvaluateGaussPointsToNodesExtrapolationMatrix<distype>(integration);
  }

}  // namespace
