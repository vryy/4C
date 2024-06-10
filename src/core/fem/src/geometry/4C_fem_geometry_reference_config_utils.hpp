/*---------------------------------------------------------------------*/
/*! \file

 \brief Utility functions to calculate properties in the reference configuration

 \level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GEOMETRY_REFERENCE_CONFIG_UTILS_HPP
#define FOUR_C_FEM_GEOMETRY_REFERENCE_CONFIG_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /**
   * \brief Calculate the global position in the reference configuration for a given element at a
   * given position in parameter space \f$ \vec{\xi} \f$
   *
   *  \param [in]   element: element for which position shall be calculated
   *  \param [in]        xi: position in parameter space
   *  \param [in,out] coord: global position in reference configuration for element at parameter
   * space position xi
   *
   *  \author cschmidt \date 11/18 */
  template <int probdim, Core::FE::CellType distype>
  static void LocalToGlobalPositionAtXiRefConfig(const Core::Elements::Element* element,
      const Core::LinAlg::Matrix<Core::FE::dim<distype>, 1>& xi,
      Core::LinAlg::Matrix<probdim, 1>& coord)
  {
    static Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1> funct(true);
    static Core::LinAlg::Matrix<probdim, Core::FE::num_nodes<distype>> nodecoords(true);

    Core::FE::shape_function<distype>(xi, funct);

    const Core::Nodes::Node* const* nodes = element->Nodes();
    const int nodedim = nodes[0]->Dim();

    if (!nodes)
    {
      FOUR_C_THROW("ERROR: Did not get nodes of element!");
    }
    if (probdim != nodedim)
    {
      FOUR_C_THROW(
          "Problem dimension: %i and dimension of nodes: %i does not match!", probdim, nodedim);
    }

    for (int i = 0; i < Core::FE::num_nodes<distype>; ++i)
    {
      for (int j = 0; j < nodedim; ++j)
      {
        nodecoords(j, i) = nodes[i]->X()[j];
      }
    }

    coord.Multiply(1.0, nodecoords, funct, 0.0);
  }

  /**
   * \brief Calculate the normal in the reference configuration for a given element at a given
   * position in parameter space \f$ \vec{\xi} \f$
   *
   *  \param [in]    element: element for which position shall be calculated
   *  \param [in]         xi: position in parameter space
   *  \param [in,out] normal: normal in reference configuration for element at parameter space
   * position xi
   *
   *  \author cschmidt \date 11/18 */
  template <Core::FE::CellType distype>
  static void ComputeUnitNormalAtXiRefConfig(const Core::Elements::Element* element,
      const Core::LinAlg::Matrix<Core::FE::dim<distype>, 1>& xi, Core::LinAlg::Matrix<3, 1>& normal)
  {
    static Core::LinAlg::Matrix<3, Core::FE::dim<distype>> gxieta(true);
    static Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>> deriv(true);
    static Core::LinAlg::Matrix<3, Core::FE::num_nodes<distype>> nodecoords(true);

    Core::FE::shape_function_deriv1<distype>(xi, deriv);

    const Core::Nodes::Node* const* nodes = element->Nodes();
    const int nodedim = nodes[0]->Dim();

    if (!nodes)
    {
      FOUR_C_THROW("ERROR: Did not get nodes of element!");
    }
    if (nodedim != 3)
    {
      FOUR_C_THROW("ERROR: Only implemented for 3D cases so far!");
    }

    for (int i = 0; i < Core::FE::num_nodes<distype>; ++i)
    {
      for (int j = 0; j < nodedim; ++j)
      {
        nodecoords(j, i) = nodes[i]->X()[j];
      }
    }

    gxieta.MultiplyNT(1.0, nodecoords, deriv, 0.0);
    static Core::LinAlg::Matrix<3, 1> gxi(true);
    static Core::LinAlg::Matrix<3, 1> geta(true);
    static Core::LinAlg::Matrix<2, 1> first(true);
    static Core::LinAlg::Matrix<2, 1> second(true);
    first(0, 0) = 1.0;
    second(1, 0) = 1.0;
    gxi.Multiply(1.0, gxieta, first, 0.0);
    geta.Multiply(1.0, gxieta, second, 0.0);

    // clear, calculate and scale normal
    normal.Clear();
    normal.CrossProduct(gxi, geta);
    const double normnormal = normal.Norm2();
    normal.Scale(1.0 / normnormal);
  }
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
