/*! \file
\brief Helpers for fluid element nullspace computation
\level 0
*/
#include "baci_fluid_ele_nullspace.H"

#include "baci_lib_node.H"
#include "baci_utils_exceptions.H"

namespace FLD
{
  Teuchos::SerialDenseMatrix<int, double> ComputeFluidNullSpace(
      const DRT::Node& node, const int numdof, const int dimnsp)
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

    if (numdof > 10) dserror("Cannot define more than 10 degrees of freedom!");

    Teuchos::SerialDenseMatrix<int, double> nullspace(numdof, dimnsp);
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