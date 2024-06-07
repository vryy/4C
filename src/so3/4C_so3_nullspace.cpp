/*! \file
\brief Helpers for solid element nullspace computation
\level 0
*/

#include "4C_so3_nullspace.hpp"

#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  Core::LinAlg::SerialDenseMatrix ComputeSolid3DNullSpace(
      const Core::Nodes::Node& node, const double* x0)
  {
    /* the rigid body modes for structures are:

          xtrans   ytrans  ztrans   xrot       yrot       zrot
          mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
      -----------------------------------------------------------
    x   |    1       0       0       0          z-z0      -y+y0
    y   |    0       1       0      -z+z0       0          x-x0
    z   |    0       0       1       y-y0      -x+x0       0

    valid element types: ale3, so_hex8, so_hex20, so_hex27, sosh8, so_tet4,
                         so_tet10, so_weg6, sodisp, so_shw6, truss3, torsion3
    */

    const auto& x = node.X();

    Core::LinAlg::SerialDenseMatrix nullspace(3, 6);
    // x-modes
    nullspace(0, 0) = 1.0;
    nullspace(0, 1) = 0.0;
    nullspace(0, 2) = 0.0;
    nullspace(0, 3) = 0.0;
    nullspace(0, 4) = x[2] - x0[2];
    nullspace(0, 5) = -x[1] + x0[1];
    // y-modes
    nullspace(1, 0) = 0.0;
    nullspace(1, 1) = 1.0;
    nullspace(1, 2) = 0.0;
    nullspace(1, 3) = -x[2] + x0[2];
    nullspace(1, 4) = 0.0;
    nullspace(1, 5) = x[0] - x0[0];
    // z-modes
    nullspace(2, 0) = 0.0;
    nullspace(2, 1) = 0.0;
    nullspace(2, 2) = 1.0;
    nullspace(2, 3) = x[1] - x0[1];
    nullspace(2, 4) = -x[0] + x0[0];
    nullspace(2, 5) = 0.0;

    return nullspace;
  }

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  Core::LinAlg::SerialDenseMatrix ComputeSolid2DNullSpace(
      const Core::Nodes::Node& node, const double* x0)
  {
    /* the rigid body modes for structures are:

          xtrans   ytrans   zrot
          mode[0]  mode[1]  mode[3]
        ----------------------------
    x   |    1       0       -y+y0
    y   |    0       1       x-x0

    valid element types: wall1, ale2, torsion2

     */

    const auto& x = node.X();

    Core::LinAlg::SerialDenseMatrix nullspace(2, 3);
    // x-modes
    nullspace(0, 0) = 1.0;
    nullspace(0, 1) = 0.0;
    nullspace(0, 2) = -x[1] + x0[1];
    // y-modes
    nullspace(1, 0) = 0.0;
    nullspace(1, 1) = 1.0;
    nullspace(1, 2) = x[0] - x0[0];

    return nullspace;
  }
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
