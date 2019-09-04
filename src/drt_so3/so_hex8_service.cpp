/*----------------------------------------------------------------------*/
/*! \file

\brief Service routines for Solid Hex8 element

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/
#include "so_hex8.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_node.H"


void DRT::ELEMENTS::So_hex8::soh8_ElementCenterRefeCoords(
    LINALG::Matrix<1, NUMDIM_SOH8>& centercoord,
    LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe) const
{
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOH8, 1> funct;
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  centercoord.MultiplyTN(funct, xrefe);
  return;
}


void DRT::ELEMENTS::So_hex8::soh8_GaussPointRefeCoords(LINALG::Matrix<1, NUMDIM_SOH8>& gpcoord,
    LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe, int const gp) const
{
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  LINALG::Matrix<NUMNOD_SOH8, 1> funct(true);
  funct = shapefcts[gp];
  gpcoord.MultiplyTN(funct, xrefe);

  return;
}
