/*----------------------------------------------------------------------*/
/*! \file

\brief Service routines for Solid Hex8 element

\level 1


*----------------------------------------------------------------------*/
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_so3_hex8.hpp"

FOUR_C_NAMESPACE_OPEN


void Discret::ELEMENTS::SoHex8::soh8_element_center_refe_coords(
    Core::LinAlg::Matrix<NUMDIM_SOH8, 1>& centercoord,
    Core::LinAlg::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe) const
{
  const Core::FE::CellType distype = Shape();
  Core::LinAlg::Matrix<NUMNOD_SOH8, 1> funct;
  Core::FE::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  centercoord.multiply_tn(xrefe, funct);
  return;
}


void Discret::ELEMENTS::SoHex8::soh8_gauss_point_refe_coords(
    Core::LinAlg::Matrix<NUMDIM_SOH8, 1>& gpcoord,
    Core::LinAlg::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe, int const gp) const
{
  const static std::vector<Core::LinAlg::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  Core::LinAlg::Matrix<NUMNOD_SOH8, 1> funct(true);
  funct = shapefcts[gp];
  gpcoord.multiply_tn(xrefe, funct);

  return;
}

FOUR_C_NAMESPACE_CLOSE
