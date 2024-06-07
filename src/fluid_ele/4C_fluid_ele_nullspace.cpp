/*! \file
\brief Helpers for fluid element nullspace computation
\level 0
*/
#include "4C_fluid_ele_nullspace.hpp"

#include "4C_fem_general_node.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  Core::LinAlg::SerialDenseMatrix ComputeFluidNullSpace(
      const Core::Nodes::Node& node, const int numdof, const int dimnsp)
  {
    /* the rigid body modes for fluids are:

    xtrans   ytrans  ztrans   pressure
            mode[0]  mode[1] mode[2]  mode[3]
        ----------------------------------------
        x   |    1       0       0       0
        y   |    0       1       0       0
        z   |    0       0       1       0
        p   |    0       0       0       1

        valid element types: fluid3, xfluid3
                */

    if (numdof > 10) FOUR_C_THROW("Cannot define more than 10 degrees of freedom!");

    Core::LinAlg::SerialDenseMatrix nullspace(numdof, dimnsp);
    for (int i = 0; i < numdof; i++)
    {
      for (int j = 0; j < dimnsp; j++)
      {
        if (i == j)
          nullspace(i, j) = 1.0;
        else
          nullspace(i, j) = 0.0;
      }
    }

    return nullspace;
  }
}  // namespace FLD
FOUR_C_NAMESPACE_CLOSE
