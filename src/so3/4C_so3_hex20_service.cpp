/*----------------------------------------------------------------------*/
/*! \file
\brief 3D quadratic serendipity element
\level 1

*----------------------------------------------------------------------*/
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_so3_hex20.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  return Center Coords in Reference System                            |
 *----------------------------------------------------------------------*/
std::vector<double> Discret::ELEMENTS::SoHex20::soh20_element_center_refe_coords()
{
  // update element geometry
  Core::Nodes::Node** nodes = Nodes();
  Core::LinAlg::Matrix<NUMNOD_SOH20, NUMDIM_SOH20> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH20; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const Core::FE::CellType distype = Shape();
  Core::LinAlg::Matrix<NUMNOD_SOH20, 1> funct;
  // Element midpoint at r=s=t=0.0
  Core::FE::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  Core::LinAlg::Matrix<1, NUMDIM_SOH20> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}

FOUR_C_NAMESPACE_CLOSE
