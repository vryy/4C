/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the element calculation library

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_solid_3D_ele_calc_lib.hpp"

namespace
{
  using namespace FourC;

  TEST(EvaluateParameterCoordinateCentroid, DisTypeHex)
  {
    // only tested for hex8, but equivalent for hex18, hex27, ...
    const auto distype = Core::FE::CellType::hex8;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        Discret::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeTet)
  {
    // only tested for tet4, but equivalent for tet10
    const auto distype = Core::FE::CellType::tet4;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(0) = 0.25;
    xi_centroid_ref(1) = 0.25;
    xi_centroid_ref(2) = 0.25;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        Discret::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = Core::FE::CellType::pyramid5;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(2) = 0.25;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        Discret::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = Core::FE::CellType::wedge6;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid_ref(true);
    xi_centroid_ref(0) = 1.0 / 3.0;
    xi_centroid_ref(1) = 1.0 / 3.0;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> xi_centroid =
        Discret::ELEMENTS::EvaluateParameterCoordinateCentroid<distype>();

    EXPECT_EQ(xi_centroid, xi_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeHex)
  {
    const auto distype = Core::FE::CellType::hex8;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1>
        reference_coords_centroid_ref(true);

    Discret::ELEMENTS::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(0, 1) = 0;
    nodal_coordinates.reference_coordinates(0, 2) = 0;

    nodal_coordinates.reference_coordinates(1, 0) = 4;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(1, 2) = 0;

    nodal_coordinates.reference_coordinates(2, 0) = 4;
    nodal_coordinates.reference_coordinates(2, 1) = 1;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(3, 0) = 0;
    nodal_coordinates.reference_coordinates(3, 1) = 1;
    nodal_coordinates.reference_coordinates(3, 2) = 0;

    nodal_coordinates.reference_coordinates(4, 0) = 0;
    nodal_coordinates.reference_coordinates(4, 1) = 0;
    nodal_coordinates.reference_coordinates(4, 2) = 2;

    nodal_coordinates.reference_coordinates(5, 0) = 4;
    nodal_coordinates.reference_coordinates(5, 1) = 0;
    nodal_coordinates.reference_coordinates(5, 2) = 2;

    nodal_coordinates.reference_coordinates(6, 0) = 4;
    nodal_coordinates.reference_coordinates(6, 1) = 1;
    nodal_coordinates.reference_coordinates(6, 2) = 2;

    nodal_coordinates.reference_coordinates(7, 0) = 0;
    nodal_coordinates.reference_coordinates(7, 1) = 1;
    nodal_coordinates.reference_coordinates(7, 2) = 2;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid_ref(true);
    x_centroid_ref(0) = 2;
    x_centroid_ref(1) = 0.5;
    x_centroid_ref(2) = 1.0;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid =
        Discret::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeTet)
  {
    const auto distype = Core::FE::CellType::tet4;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1>
        reference_coords_centroid_ref(true);

    Discret::ELEMENTS::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(0, 1) = 0;
    nodal_coordinates.reference_coordinates(0, 2) = 0;

    nodal_coordinates.reference_coordinates(1, 0) = 1;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(1, 2) = 0;

    nodal_coordinates.reference_coordinates(2, 0) = 0;
    nodal_coordinates.reference_coordinates(2, 1) = 2;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(3, 0) = 0;
    nodal_coordinates.reference_coordinates(3, 1) = 0;
    nodal_coordinates.reference_coordinates(3, 2) = 4;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid_ref(true);
    x_centroid_ref(0) = 0.25;
    x_centroid_ref(1) = 0.5;
    x_centroid_ref(2) = 1.0;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid =
        Discret::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = Core::FE::CellType::pyramid5;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1>
        reference_coords_centroid_ref(true);

    Discret::ELEMENTS::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = -2;
    nodal_coordinates.reference_coordinates(0, 1) = -1;
    nodal_coordinates.reference_coordinates(0, 2) = 0;

    nodal_coordinates.reference_coordinates(1, 0) = 4;
    nodal_coordinates.reference_coordinates(1, 1) = -1;
    nodal_coordinates.reference_coordinates(1, 2) = 0;

    nodal_coordinates.reference_coordinates(2, 0) = 4;
    nodal_coordinates.reference_coordinates(2, 1) = 1;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(3, 0) = -2;
    nodal_coordinates.reference_coordinates(3, 1) = 1;
    nodal_coordinates.reference_coordinates(3, 2) = 0;

    nodal_coordinates.reference_coordinates(4, 0) = 1;
    nodal_coordinates.reference_coordinates(4, 1) = 0;
    nodal_coordinates.reference_coordinates(4, 2) = 4;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid_ref(true);
    x_centroid_ref(0) = 1.0;
    x_centroid_ref(1) = 0.0;
    x_centroid_ref(2) = 1.0;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid =
        Discret::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates);

    EXPECT_EQ(x_centroid, x_centroid_ref);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = Core::FE::CellType::wedge6;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1>
        reference_coords_centroid_ref(true);

    Discret::ELEMENTS::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(0, 1) = 0;
    nodal_coordinates.reference_coordinates(0, 2) = 0;

    nodal_coordinates.reference_coordinates(1, 0) = 3;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(1, 2) = 0;

    nodal_coordinates.reference_coordinates(2, 0) = 0;
    nodal_coordinates.reference_coordinates(2, 1) = 6;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(3, 0) = 0;
    nodal_coordinates.reference_coordinates(3, 1) = 0;
    nodal_coordinates.reference_coordinates(3, 2) = 1;

    nodal_coordinates.reference_coordinates(4, 0) = 3;
    nodal_coordinates.reference_coordinates(4, 1) = 0;
    nodal_coordinates.reference_coordinates(4, 2) = 1;

    nodal_coordinates.reference_coordinates(5, 0) = 0;
    nodal_coordinates.reference_coordinates(5, 1) = 6;
    nodal_coordinates.reference_coordinates(5, 2) = 1;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid_ref(true);
    x_centroid_ref(0) = 1;
    x_centroid_ref(1) = 2;
    x_centroid_ref(2) = 0.5;

    Core::LinAlg::Matrix<Discret::ELEMENTS::DETAIL::num_dim<distype>, 1> x_centroid =
        Discret::ELEMENTS::EvaluateReferenceCoordinateCentroid<distype>(nodal_coordinates);

    for (int j = 0; j < Discret::ELEMENTS::DETAIL::num_dim<distype>; j++)
    {
      EXPECT_NEAR(x_centroid(j), x_centroid_ref(j), 1e-14);
    }
  }
}  // namespace