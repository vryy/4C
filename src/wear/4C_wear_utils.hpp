/*----------------------------------------------------------------------*/
/*! \file

\brief  utils for wear algorithm

\level 2

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_WEAR_UTILS_HPP
#define FOUR_C_WEAR_UTILS_HPP

/*----------------------------------------------------------------------*
 | headers                                                  farah 12/13 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_wear_defines.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | forward declarations                                     farah 12/13 |
 *----------------------------------------------------------------------*/
namespace Discret
{
  class LocationArray;
}
/*----------------------------------------------------------------------*
 |                                                          farah 12/13 |
 *----------------------------------------------------------------------*/
namespace Wear
{
  namespace UTILS
  {
    // advection map for elements
    template <Core::FE::CellType distype>
    void av(Core::Elements::Element* ele,               // in
        double* Xtarget,                                // out
        double* Xsource,                                // in
        Teuchos::RCP<const Epetra_Vector> disp_source,  // in
        Teuchos::RCP<const Epetra_Vector> disp_target,  // in
        const std::vector<int>& lm,                     // in
        bool& found, double* e)
    {
      static constexpr int numnod = Core::FE::num_nodes<distype>;
      static constexpr int ndim = Core::FE::dim<distype>;

      Core::LinAlg::Matrix<numnod, 1> funct;
      Core::LinAlg::Matrix<ndim, numnod> xcure;
      Core::LinAlg::Matrix<ndim, ndim> xjm;
      Core::LinAlg::Matrix<ndim, numnod> deriv;

      // spatial displacements
      std::vector<double> mydisp_source(lm.size());
      Core::FE::ExtractMyValues(*disp_source, mydisp_source, lm);

      // material displacements
      std::vector<double> mydisp_target(lm.size());
      Core::FE::ExtractMyValues(*disp_target, mydisp_target, lm);

      // spatial configuration of this element!
      for (int k = 0; k < numnod; ++k)
        for (int j = 0; j < ndim; ++j)
          xcure(j, k) = ele->Nodes()[k]->X()[j] + mydisp_source[k * ndim + j];

      // first estimation for parameter space coordinates
      for (int p = 0; p < 3; ++p) e[p] = 0.0;

      double rhs[ndim];

      // converged
      bool converged = false;
      int j = 0;

      //************************************************
      // loop
      while (!converged and j < 10)
      {
        // reset matriced
        xjm.Clear();
        deriv.Clear();

        if (ndim == 2)
        {
          Core::FE::shape_function_2D(funct, e[0], e[1], distype);
          Core::FE::shape_function_2D_deriv1(deriv, e[0], e[1], distype);
        }
        else if (ndim == 3)
        {
          Core::FE::shape_function_3D(funct, e[0], e[1], e[2], distype);
          Core::FE::shape_function_3D_deriv1(deriv, e[0], e[1], e[2], distype);
        }
        else
          FOUR_C_THROW("Wrong dimension!");

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p)
            for (int l = 0; l < ndim; ++l) xjm(p, l) += deriv(l, k) * xcure(p, k);

        // rhs of (linearized equation)
        for (int p = 0; p < ndim; ++p) rhs[p] = 0.0;

        for (int p = 0; p < ndim; ++p) rhs[p] = -Xsource[p];

        for (int k = 0; k < numnod; ++k)
          for (int p = 0; p < ndim; ++p) rhs[p] += funct(k) * xcure(p, k);

        double norm = 0.0;
        for (int p = 0; p < ndim; ++p) norm += rhs[p] * rhs[p];

        if (sqrt(norm) < WEARCONV) converged = true;

        // solve equation
        if (abs(xjm.Determinant()) < WEARSING) FOUR_C_THROW("*** WARNING: jacobi singular ***");

        double xjm_invert = xjm.Invert();
        if (abs(xjm_invert) < WEARSING) FOUR_C_THROW("Singular Jacobian for advection map");

        double deltae[3];
        for (int p = 0; p < ndim; ++p) deltae[p] = 0.0;

        for (int z = 0; z < ndim; ++z)
          for (int p = 0; p < ndim; ++p) deltae[z] -= xjm(z, p) * rhs[p];

        // incremental update
        for (int p = 0; p < ndim; ++p) e[p] += deltae[p];

        j = j + 1;
      }  // end loop
      //************************************************

      if (!converged) FOUR_C_THROW("Evaluation of element coordinates not converged!");

      // if material parameters are within the element, evaluate material
      // coordinates
      if (distype == Core::FE::CellType::hex8 or distype == Core::FE::CellType::hex20 or
          distype == Core::FE::CellType::hex27 or distype == Core::FE::CellType::quad4 or
          distype == Core::FE::CellType::quad8 or distype == Core::FE::CellType::quad9)
      {
        if (e[0] >= -1.0 - WEARADVMAP and e[0] <= 1.0 + WEARADVMAP and e[1] >= -1.0 - WEARADVMAP and
            e[1] <= 1.0 + WEARADVMAP and e[2] >= -1.0 - WEARADVMAP and e[2] <= 1.0 + WEARADVMAP)
          found = true;
      }
      else if (distype == Core::FE::CellType::tet4 or distype == Core::FE::CellType::tet10 or
               distype == Core::FE::CellType::tri3 or distype == Core::FE::CellType::tri6)
      {
        if (e[0] >= 0.0 - WEARADVMAP and e[0] <= 1.0 + WEARADVMAP and e[1] >= 0.0 - WEARADVMAP and
            e[1] <= 1.0 + WEARADVMAP and e[2] >= 0.0 - WEARADVMAP and e[2] <= 1.0 + WEARADVMAP)
          found = true;
      }
      else
        FOUR_C_THROW("Element type not supported!");

      double xmat[ndim];
      for (int p = 0; p < ndim; ++p) xmat[p] = 0.0;

      if (ndim == 2)
        Core::FE::shape_function_2D(funct, e[0], e[1], distype);
      else
        Core::FE::shape_function_3D(funct, e[0], e[1], e[2], distype);

      for (int k = 0; k < numnod; ++k)
        for (int p = 0; p < ndim; ++p)
          xmat[p] += funct(k) * (ele->Nodes()[k]->X()[p] + mydisp_target[k * ndim + p]);

      for (int p = 0; p < ndim; ++p) Xtarget[p] = xmat[p];

      return;
    };

  }  // namespace UTILS

}  // namespace Wear

FOUR_C_NAMESPACE_CLOSE

#endif
