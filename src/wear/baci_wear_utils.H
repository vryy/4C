/*----------------------------------------------------------------------*/
/*! \file

\brief  utils for wear algorithm

\level 2

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef BACI_WEAR_UTILS_H
#define BACI_WEAR_UTILS_H

/*----------------------------------------------------------------------*
 | headers                                                  farah 12/13 |
 *----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_lib_element.H"
#include "baci_lib_utils.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_wear_defines.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | forward declarations                                     farah 12/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class LocationArray;
}
/*----------------------------------------------------------------------*
 |                                                          farah 12/13 |
 *----------------------------------------------------------------------*/
namespace WEAR
{
  namespace UTILS
  {
    // advection map for elements
    template <CORE::FE::CellType distype>
    void av(DRT::Element* ele,                          // in
        double* Xtarget,                                // out
        double* Xsource,                                // in
        Teuchos::RCP<const Epetra_Vector> disp_source,  // in
        Teuchos::RCP<const Epetra_Vector> disp_target,  // in
        const std::vector<int>& lm,                     // in
        bool& found, double* e)
    {
      static const int numnod = CORE::FE::num_nodes<distype>;
      static const int ndim = CORE::FE::dim<distype>;

      CORE::LINALG::Matrix<numnod, 1> funct;
      CORE::LINALG::Matrix<ndim, numnod> xcure;
      CORE::LINALG::Matrix<ndim, ndim> xjm;
      CORE::LINALG::Matrix<ndim, numnod> deriv;

      // spatial displacements
      std::vector<double> mydisp_source(lm.size());
      DRT::UTILS::ExtractMyValues(*disp_source, mydisp_source, lm);

      // material displacements
      std::vector<double> mydisp_target(lm.size());
      DRT::UTILS::ExtractMyValues(*disp_target, mydisp_target, lm);

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
          CORE::FE::shape_function_2D(funct, e[0], e[1], distype);
          CORE::FE::shape_function_2D_deriv1(deriv, e[0], e[1], distype);
        }
        else if (ndim == 3)
        {
          CORE::FE::shape_function_3D(funct, e[0], e[1], e[2], distype);
          CORE::FE::shape_function_3D_deriv1(deriv, e[0], e[1], e[2], distype);
        }
        else
          dserror("Wrong dimension!");

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
        if (abs(xjm.Determinant()) < WEARSING) dserror("*** WARNING: jacobi singular ***");

        double xjm_invert = xjm.Invert();
        if (abs(xjm_invert) < WEARSING) dserror("Singular Jacobian for advection map");

        double deltae[3];
        for (int p = 0; p < ndim; ++p) deltae[p] = 0.0;

        for (int z = 0; z < ndim; ++z)
          for (int p = 0; p < ndim; ++p) deltae[z] -= xjm(z, p) * rhs[p];

        // incremental update
        for (int p = 0; p < ndim; ++p) e[p] += deltae[p];

        j = j + 1;
      }  // end loop
      //************************************************

      if (!converged) dserror("Evaluation of element coordinates not converged!");

      // if material parameters are within the element, evaluate material
      // coordinates
      if (distype == CORE::FE::CellType::hex8 or distype == CORE::FE::CellType::hex20 or
          distype == CORE::FE::CellType::hex27 or distype == CORE::FE::CellType::quad4 or
          distype == CORE::FE::CellType::quad8 or distype == CORE::FE::CellType::quad9)
      {
        if (e[0] >= -1.0 - WEARADVMAP and e[0] <= 1.0 + WEARADVMAP and e[1] >= -1.0 - WEARADVMAP and
            e[1] <= 1.0 + WEARADVMAP and e[2] >= -1.0 - WEARADVMAP and e[2] <= 1.0 + WEARADVMAP)
          found = true;
      }
      else if (distype == CORE::FE::CellType::tet4 or distype == CORE::FE::CellType::tet10 or
               distype == CORE::FE::CellType::tri3 or distype == CORE::FE::CellType::tri6)
      {
        if (e[0] >= 0.0 - WEARADVMAP and e[0] <= 1.0 + WEARADVMAP and e[1] >= 0.0 - WEARADVMAP and
            e[1] <= 1.0 + WEARADVMAP and e[2] >= 0.0 - WEARADVMAP and e[2] <= 1.0 + WEARADVMAP)
          found = true;
      }
      else
        dserror("Element type not supported!");

      double xmat[ndim];
      for (int p = 0; p < ndim; ++p) xmat[p] = 0.0;

      if (ndim == 2)
        CORE::FE::shape_function_2D(funct, e[0], e[1], distype);
      else
        CORE::FE::shape_function_3D(funct, e[0], e[1], e[2], distype);

      for (int k = 0; k < numnod; ++k)
        for (int p = 0; p < ndim; ++p)
          xmat[p] += funct(k) * (ele->Nodes()[k]->X()[p] + mydisp_target[k * ndim + p]);

      for (int p = 0; p < ndim; ++p) Xtarget[p] = xmat[p];

      return;
    };

  }  // namespace UTILS

}  // namespace WEAR

BACI_NAMESPACE_CLOSE

#endif  // WEAR_UTILS_H
