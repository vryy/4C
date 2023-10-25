/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the element calculation library

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_solid_ele_calc_lib.H"

namespace
{
  TEST(EvaluateParameterCoordinateCentroid, DisTypeHex)
  {
    // only tested for hex8, but equivalent for hex18, hex27, ...
    const auto distype = CORE::FE::CellType::hex8;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        DRT::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeTet)
  {
    // only tested for tet4, but equivalent for tet10
    const auto distype = CORE::FE::CellType::tet4;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(0) = 0.25;
    xi_centroid_ref(1) = 0.25;
    xi_centroid_ref(2) = 0.25;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        DRT::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = CORE::FE::CellType::pyramid5;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(2) = 0.25;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        DRT::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = CORE::FE::CellType::wedge6;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(0) = 1.0 / 3.0;
    xi_centroid_ref(1) = 1.0 / 3.0;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        DRT::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeHex)
  {
    const auto distype = CORE::FE::CellType::hex8;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> reference_coords_centroid_ref(
        true);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_nodes<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        nodal_coordinates_reference(true);

    nodal_coordinates_reference(0, 0) = 0;
    nodal_coordinates_reference(0, 1) = 0;
    nodal_coordinates_reference(0, 2) = 0;

    nodal_coordinates_reference(1, 0) = 4;
    nodal_coordinates_reference(1, 1) = 0;
    nodal_coordinates_reference(1, 2) = 0;

    nodal_coordinates_reference(2, 0) = 4;
    nodal_coordinates_reference(2, 1) = 1;
    nodal_coordinates_reference(2, 2) = 0;

    nodal_coordinates_reference(3, 0) = 0;
    nodal_coordinates_reference(3, 1) = 1;
    nodal_coordinates_reference(3, 2) = 0;

    nodal_coordinates_reference(4, 0) = 0;
    nodal_coordinates_reference(4, 1) = 0;
    nodal_coordinates_reference(4, 2) = 2;

    nodal_coordinates_reference(5, 0) = 4;
    nodal_coordinates_reference(5, 1) = 0;
    nodal_coordinates_reference(5, 2) = 2;

    nodal_coordinates_reference(6, 0) = 4;
    nodal_coordinates_reference(6, 1) = 1;
    nodal_coordinates_reference(6, 2) = 2;

    nodal_coordinates_reference(7, 0) = 0;
    nodal_coordinates_reference(7, 1) = 1;
    nodal_coordinates_reference(7, 2) = 2;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid_ref(true);
    x_centroid_ref(0) = 2;
    x_centroid_ref(1) = 0.5;
    x_centroid_ref(2) = 1.0;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid =
        DRT::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates_reference);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeTet)
  {
    const auto distype = CORE::FE::CellType::tet4;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> reference_coords_centroid_ref(
        true);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_nodes<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        nodal_coordinates_reference(true);

    nodal_coordinates_reference(0, 0) = 0;
    nodal_coordinates_reference(0, 1) = 0;
    nodal_coordinates_reference(0, 2) = 0;

    nodal_coordinates_reference(1, 0) = 1;
    nodal_coordinates_reference(1, 1) = 0;
    nodal_coordinates_reference(1, 2) = 0;

    nodal_coordinates_reference(2, 0) = 0;
    nodal_coordinates_reference(2, 1) = 2;
    nodal_coordinates_reference(2, 2) = 0;

    nodal_coordinates_reference(3, 0) = 0;
    nodal_coordinates_reference(3, 1) = 0;
    nodal_coordinates_reference(3, 2) = 4;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid_ref(true);
    x_centroid_ref(0) = 0.25;
    x_centroid_ref(1) = 0.5;
    x_centroid_ref(2) = 1.0;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid =
        DRT::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates_reference);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = CORE::FE::CellType::pyramid5;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> reference_coords_centroid_ref(
        true);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_nodes<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        nodal_coordinates_reference(true);

    nodal_coordinates_reference(0, 0) = -2;
    nodal_coordinates_reference(0, 1) = -1;
    nodal_coordinates_reference(0, 2) = 0;

    nodal_coordinates_reference(1, 0) = 4;
    nodal_coordinates_reference(1, 1) = -1;
    nodal_coordinates_reference(1, 2) = 0;

    nodal_coordinates_reference(2, 0) = 4;
    nodal_coordinates_reference(2, 1) = 1;
    nodal_coordinates_reference(2, 2) = 0;

    nodal_coordinates_reference(3, 0) = -2;
    nodal_coordinates_reference(3, 1) = 1;
    nodal_coordinates_reference(3, 2) = 0;

    nodal_coordinates_reference(4, 0) = 1;
    nodal_coordinates_reference(4, 1) = 0;
    nodal_coordinates_reference(4, 2) = 4;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid_ref(true);
    x_centroid_ref(0) = 1.0;
    x_centroid_ref(1) = 0.0;
    x_centroid_ref(2) = 1.0;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid =
        DRT::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates_reference);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = CORE::FE::CellType::wedge6;

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_dim<distype>, 1> reference_coords_centroid_ref(
        true);

    CORE::LINALG::Matrix<DRT::ELEMENTS::DETAIL::num_nodes<distype>,
        DRT::ELEMENTS::DETAIL::num_dim<distype>>
        nodal_coordinates_reference(true);

    nodal_coordinates_reference(0, 0) = 0;
    nodal_coordinates_reference(0, 1) = 0;
    nodal_coordinates_reference(0, 2) = 0;

    nodal_coordinates_reference(1, 0) = 3;
    nodal_coordinates_reference(1, 1) = 0;
    nodal_coordinates_reference(1, 2) = 0;

    nodal_coordinates_reference(2, 0) = 0;
    nodal_coordinates_reference(2, 1) = 6;
    nodal_coordinates_reference(2, 2) = 0;

    nodal_coordinates_reference(3, 0) = 0;
    nodal_coordinates_reference(3, 1) = 0;
    nodal_coordinates_reference(3, 2) = 1;

    nodal_coordinates_reference(4, 0) = 3;
    nodal_coordinates_reference(4, 1) = 0;
    nodal_coordinates_reference(4, 2) = 1;

    nodal_coordinates_reference(5, 0) = 0;
    nodal_coordinates_reference(5, 1) = 6;
    nodal_coordinates_reference(5, 2) = 1;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid_ref(true);
    x_centroid_ref(0) = 1;
    x_centroid_ref(1) = 2;
    x_centroid_ref(2) = 0.5;

    CORE::LINALG::Matrix<1, DRT::ELEMENTS::DETAIL::num_dim<distype>> x_centroid =
        DRT::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates_reference);

    for (int j = 0; j < DRT::ELEMENTS::DETAIL::num_dim<distype>; j++)
    {
      EXPECT_NEAR(x_centroid(j), x_centroid_ref(j), 1e-14);
    }
  }
}  // namespace