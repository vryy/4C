/*! \file

\brief Utility functions for FBAR element technology

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_CALC_LIB_FBAR_HPP
#define BACI_SOLID_3D_ELE_CALC_LIB_FBAR_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type.hpp"
#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /*!
   * @brief A small container storing information needed to compute the linearization of an element
   * with FBAR
   */
  template <CORE::FE::CellType celltype>
  struct FBarLinearizationContainer
  {
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        Bop{};

    CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1> Hop{};

    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchygreen{};

    double fbar_factor = 1.0;
  };

  /*!
   * @brief Evaluate the fbar factor \f[ \frac{\mathbf{F}_{\mathrm{centroid}}}{\mathbf{F}}^{1/3}
   * \f]
   *
   * @param defgrd_centroid (in) : Deformation gradient evaluated at the element centroid
   * @param defgrd_gp (in) : Deformation gradient evaluated at the Gauss point
   * @return double : Fbar factor
   */
  inline double EvaluateFbarFactor(const double& defgrd_centroid, const double& defgrd_gp)
  {
    const double fbar_factor = std::pow(defgrd_centroid / defgrd_gp, 1.0 / 3.0);
    return fbar_factor;
  }

  /*!
   * @brief Evaluates the H-Operator used in F-bar of the specified element
   *
   * @tparam celltype : Cell type
   * @param shape_function_derivs (in) : Derivative of the shape functions w.r.t. XYZ at the Gauss
   * point
   * @param shape_function_derivs_centroid (in) : Derivative of the shape functions w.r.t. XYZ at
   * the element center
   * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient) evaluated at the Gauss point
   * @param spatial_material_mapping_centroid (in) : An object holding quantities of the spatial
   * material mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient) evaluated at the element centroid
   * @return CORE::LINALG::Matrix<num_dof_per_ele, 1> : H-Operator
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::dim<celltype> == 3, int> = 0>
  inline CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1>
  EvaluateFbarHOperator(
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::num_nodes<celltype>>&
          shape_function_derivs,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::num_nodes<celltype>>&
          shape_function_derivs_centroid,
      const DRT::ELEMENTS::SpatialMaterialMapping<celltype>& spatial_material_mapping,
      const DRT::ELEMENTS::SpatialMaterialMapping<celltype>& spatial_material_mapping_centroid)
  {
    // inverse deformation gradient at centroid
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invdefgrd_centroid;
    invdefgrd_centroid.Invert(spatial_material_mapping_centroid.deformation_gradient_);

    // inverse deformation gradient at gp
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invdefgrd;
    invdefgrd.Invert(spatial_material_mapping.deformation_gradient_);

    CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1> Hop(true);
    for (int idof = 0; idof < CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>; idof++)
    {
      for (int idim = 0; idim < CORE::FE::dim<celltype>; idim++)
      {
        Hop(idof) += invdefgrd_centroid(idim, idof % CORE::FE::dim<celltype>) *
                     shape_function_derivs_centroid(idim, idof / CORE::FE::dim<celltype>);
        Hop(idof) -= invdefgrd(idim, idof % CORE::FE::dim<celltype>) *
                     shape_function_derivs(idim, idof / CORE::FE::dim<celltype>);
      }
    }

    return Hop;
  }

  /*!
   * @brief Add fbar stiffness matrix contribution of one Gauss point
   *
   * @tparam celltype : Cell type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param Hop (in) : H-Operator
   * @param f_bar_factor (in) : f_bar_factor
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param cauchyGreen (in) : An object holding the right Cauchy-Green deformation tensor and
   * its inverse
   * @param stress_bar (in) : Deviatoric part of stress measures
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType celltype>
  inline void AddFbarStiffnessMatrix(
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2,
          CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>>& Bop,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1>& Hop,
      const double f_bar_factor, const double integration_fac,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchyGreen,
      const DRT::ELEMENTS::Stress<celltype> stress_bar,
      CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>,
          CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>>& stiffness_matrix)
  {
    constexpr int num_dof_per_ele = CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>;

    CORE::LINALG::Matrix<CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2, 1>
        rcg_bar_voigt;
    CORE::LINALG::VOIGT::Strains::MatrixToVector(cauchyGreen, rcg_bar_voigt);

    CORE::LINALG::Matrix<CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2, 1> ccg;
    ccg.MultiplyNN(stress_bar.cmat_, rcg_bar_voigt);

    // auxiliary integrated stress_bar
    CORE::LINALG::Matrix<num_dof_per_ele, 1> bopccg(false);
    bopccg.MultiplyTN(integration_fac * f_bar_factor / 3.0, Bop, ccg);

    CORE::LINALG::Matrix<num_dof_per_ele, 1> bops(false);
    bops.MultiplyTN(-integration_fac / f_bar_factor / 3.0, Bop, stress_bar.pk2_);

    for (int idof = 0; idof < num_dof_per_ele; idof++)
    {
      for (int jdof = 0; jdof < num_dof_per_ele; jdof++)
      {
        stiffness_matrix(idof, jdof) += Hop(jdof) * (bops(idof, 0) + bopccg(idof, 0));
      }
    }
  }
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE
#endif  // BACI_SOLID_3D_ELE_CALC_LIB_FBAR_HPP