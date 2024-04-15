/*----------------------------------------------------------------------*/
/*! \file

\brief Service routines for Solid Hex8 element

\level 1


*----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_lib_node.hpp"
#include "baci_so3_hex8.hpp"

FOUR_C_NAMESPACE_OPEN


void DRT::ELEMENTS::SoHex8::soh8_ElementCenterRefeCoords(
    CORE::LINALG::Matrix<NUMDIM_SOH8, 1>& centercoord,
    CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe) const
{
  const CORE::FE::CellType distype = Shape();
  CORE::LINALG::Matrix<NUMNOD_SOH8, 1> funct;
  CORE::FE::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  centercoord.MultiplyTN(xrefe, funct);
  return;
}


void DRT::ELEMENTS::SoHex8::soh8_GaussPointRefeCoords(CORE::LINALG::Matrix<NUMDIM_SOH8, 1>& gpcoord,
    CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe, int const gp) const
{
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  CORE::LINALG::Matrix<NUMNOD_SOH8, 1> funct(true);
  funct = shapefcts[gp];
  gpcoord.MultiplyTN(xrefe, funct);

  return;
}

FOUR_C_NAMESPACE_CLOSE
