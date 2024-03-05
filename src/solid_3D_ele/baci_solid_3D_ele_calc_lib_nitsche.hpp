/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_CALC_LIB_NITSCHE_HPP
#define BACI_SOLID_3D_ELE_CALC_LIB_NITSCHE_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_lib_element.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_utils_demangle.hpp"


BACI_NAMESPACE_OPEN

namespace MAT
{
  class So3Material;
}

namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    template <unsigned a, unsigned b, unsigned c>
    static inline void MultiplyTN(const CORE::LINALG::Matrix<a, b>& mat1,
        const CORE::LINALG::Matrix<a, c>& mat2, CORE::LINALG::SerialDenseMatrix& out)
    {
      out.reshape(b, c);
      CORE::LINALG::Matrix<b, c> out_mat(out.values(), true);

      out_mat.MultiplyTN(1.0, mat1, mat2, 0.0);
    }
  }  // namespace DETAILS

  /*!
   * @brief Evaluate the derivative of the Cauchy stress w.r.t. the nodal displacements
   *
   * @tparam celltype
   * @param d_F_dd (in) : derivative of the deformation gradient w.r.t. nodal displacements
   * @param d_cauchyndir_dF (in) : derivative of the cauchy stress w.r.t. deformation gradient
   * @param d_cauchyndir_dd (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateDCauchyNDirDDisplacements(
      const CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>&
          d_F_dd,
      const CORE::LINALG::Matrix<9, 1>& d_cauchyndir_dF,
      CORE::LINALG::SerialDenseMatrix& d_cauchyndir_dd)
  {
    DETAILS::MultiplyTN(d_F_dd, d_cauchyndir_dF, d_cauchyndir_dd);
  }

  /*!
   * @brief Evaluate the 2. derivative of the Cauchy stress w.r.t. the nodal displacements and the
   * normal
   *
   * @tparam celltype
   * @param d_F_dd (in) : derivative of the deformation gradient w.r.t. nodal displacements
   * @param d2_cauchyndir_dF_dn (in) : second derivative of the Cauchy stress w.r.t. the deformation
   * gradient and the normal
   * @param d2_cauchyndir_dd_dn (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateD2CauchyNDirDDisplacementsDNormal(
      const CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>&
          d_F_dd,
      const CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>& d2_cauchyndir_dF_dn,
      CORE::LINALG::SerialDenseMatrix& d2_cauchyndir_dd_dn)
  {
    DETAILS::MultiplyTN(d_F_dd, d2_cauchyndir_dF_dn, d2_cauchyndir_dd_dn);
  }

  /*!
   * @brief Evaluate the 2. derivative of the Cauchy stress w.r.t. the nodal displacements and the
   * direction
   *
   * @tparam celltype
   * @param d_F_dd (in) : derivative of the deformation gradient w.r.t. nodal displacements
   * @param d2_cauchyndir_dF_ddir (in) : second derivative of the Cauchy stress w.r.t. the
   * deformation gradient and the direction
   * @param d2_cauchyndir_dd_ddir (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateD2CauchyNDirDDisplacementsDDir(
      const CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>&
          d_F_dd,
      const CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>& d2_cauchyndir_dF_ddir,
      CORE::LINALG::SerialDenseMatrix& d2_cauchyndir_dd_ddir)
  {
    DETAILS::MultiplyTN(d_F_dd, d2_cauchyndir_dF_ddir, d2_cauchyndir_dd_ddir);
  }

  /*!
   * @brief Evaluate the 2. derivative of the Cauchy stress w.r.t. the nodal displacements
   *
   * @tparam celltype
   * @param d_F_dd (in) : derivative of the deformation gradient w.r.t. nodal displacements
   * @param d2_cauchyndir_dF2 (in) : second derivative of the Cauchy stress w.r.t. the
   * deformation gradient
   * @param d2_cauchyndir_dd_dd (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateD2CauchyNDirDDisplacements2(
      const CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>&
          d_F_dd,
      const CORE::LINALG::Matrix<9, 9>& d2_cauchyndir_dF2,
      CORE::LINALG::SerialDenseMatrix& d2_cauchyndir_dd_dd)
  {
    CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        d2_cauchyndir_dF2_d_F_dd(false);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    DETAILS::MultiplyTN(d_F_dd, d2_cauchyndir_dF2_d_F_dd, d2_cauchyndir_dd_dd);
  }

  /*!
   * @brief Evaluate the derivateive of the Cauchy stress w.r.t. xi
   *
   * @tparam celltype
   * @param d_F_dxi  (in) : derivative of the deformation gradient w.r.t. xi
   * @param d_cauchyndir_dF (in) : derivative of the Cauchy stress w.r.t. deformation gradient
   * @param d_cauchyndir_dxi (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateDCauchyNDirDXi(const CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>& d_F_dxi,
      const CORE::LINALG::Matrix<9, 1>& d_cauchyndir_dF,
      CORE::LINALG::Matrix<3, 1>& d_cauchyndir_dxi)
  {
    d_cauchyndir_dxi.MultiplyTN(1.0, d_F_dxi, d_cauchyndir_dF, 0.0);
  }

  /*!
   * @brief Evaluate the second derivative of the Cauchy stress w.r.t. the nodal displacements and
   * xi
   *
   * @tparam celltype
   * @param d2_F_dxi_dd (in) : second derivative of the deformation gradient w.r.t. xi and the nodal
   * displacements
   * @param d_cauchyndir_dF (in) : derivative of the Cauchy stress w.r.t. deformation gradient
   * @param d2_cauchyndir_dd_dxi (out) : result
   */
  template <CORE::FE::CellType celltype>
  void EvaluateD2CauchyNDirDDisplacementsDXi(
      const CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype> *
                                        CORE::FE::dim<celltype>>& d2_F_dxi_dd,
      const CORE::LINALG::Matrix<9, 1>& d_cauchyndir_dF,
      CORE::LINALG::SerialDenseMatrix& d2_cauchyndir_dd_dxi)
  {
    constexpr int num_dof = CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>;

    d2_cauchyndir_dd_dxi.reshape(num_dof, CORE::FE::dim<celltype>);

    for (int i = 0; i < CORE::FE::dim<celltype>; ++i)
    {
      for (int j = 0; j < CORE::FE::dim<celltype>; ++j)
      {
        for (int k = 0; k < CORE::FE::num_nodes<celltype>; ++k)
        {
          for (int l = 0; l < CORE::FE::dim<celltype>; ++l)
          {
            using VoigtMapping = CORE::LINALG::VOIGT::IndexMappings;
            d2_cauchyndir_dd_dxi(k * 3 + i, l) +=
                d_cauchyndir_dF(VoigtMapping::NonSymToVoigt9(i, j), 0) *
                d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                    CORE::FE::dim<celltype> * (CORE::FE::dim<celltype> * k + i) + l);
          }
        }
      }
    }
  }

  /*!
   * @brief A struct holding the Cauchy stress at a point xi in a specific direction and the
   * derivatives w.r.t. all input quantities
   *
   * @tparam dim
   */
  template <int dim>
  struct CauchyNDirAndLinearization
  {
    /// Cauchy stress in a specific direction
    double cauchy_n_dir = 0.0;

    /// @brief first derivatives
    ///@{

    /// first derivative w.r.t. normal
    CORE::LINALG::Matrix<dim, 1> d_cauchyndir_dn{};

    /// first derivative w.r.t. given direction
    CORE::LINALG::Matrix<dim, 1> d_cauchyndir_ddir{};

    /// first derivative w.r.t. xi
    CORE::LINALG::Matrix<dim, 1> d_cauchyndir_dxi{};

    /// first derivative w.r.t. displacements
    CORE::LINALG::SerialDenseMatrix d_cauchyndir_dd;
    ///@}

    /// @brief second derivatives
    ///@{

    /// second derivative w.r.t. displacements
    CORE::LINALG::SerialDenseMatrix d2_cauchyndir_dd2;

    /// second derivative w.r.t. displacements, normal
    CORE::LINALG::SerialDenseMatrix d2_cauchyndir_dd_dn;

    /// second derivative w.r.t. displacements, given direction
    CORE::LINALG::SerialDenseMatrix d2_cauchyndir_dd_ddir;

    /// second derivative w.r.t. displacements, xi
    CORE::LINALG::SerialDenseMatrix d2_cauchyndir_dd_dxi;
    ///@}
  };

  /// Check wheter a solid variant can evaluate the Cauchy stress at xi in a specifc direction
  /// including the derivatives
  template <typename T, int dim, typename AlwaysVoid = void>
  constexpr bool CanEvaluateEvaluateCauchyNDir = false;

  template <typename T, int dim>
  constexpr bool CanEvaluateEvaluateCauchyNDir<T, dim,
      std::void_t<decltype(std::declval<T>()->GetCauchyNDirAndDerivativesAtXi(
          std::declval<const DRT::Element&>(), std::declval<MAT::So3Material&>(),
          std::declval<const std::vector<double>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>()))>> = true;

  namespace DETAILS
  {
    template <int dim>
    struct EvaluateCauchyNDirAction
    {
      EvaluateCauchyNDirAction(const DRT::Element& e, MAT::So3Material& m,
          const std::vector<double>& d, const CORE::LINALG::Matrix<dim, 1>& x,
          const CORE::LINALG::Matrix<dim, 1>& normal, const CORE::LINALG::Matrix<dim, 1>& direction)
          : element(e), mat(m), disp(d), xi(x), n(normal), dir(direction)
      {
      }

      template <typename T, std::enable_if_t<CanEvaluateEvaluateCauchyNDir<T&, dim>, bool> = true>
      CauchyNDirAndLinearization<dim> operator()(T& cauchy_n_dir_evaluatable)
      {
        return cauchy_n_dir_evaluatable->GetCauchyNDirAndDerivativesAtXi(
            element, mat, disp, xi, n, dir);
      }

      template <typename T, std::enable_if_t<!CanEvaluateEvaluateCauchyNDir<T&, dim>, bool> = true>
      CauchyNDirAndLinearization<dim> operator()(T& other)
      {
        dserror(
            "Your element evaluation %s does not allow to evaluate the Cauchy stress at a specific "
            "point in a specific direction in the dimension dim=%d.",
            CORE::UTILS::TryDemangle(typeid(T).name()).c_str(), dim);
      }

      const DRT::Element& element;
      MAT::So3Material& mat;
      const std::vector<double>& disp;
      const CORE::LINALG::Matrix<dim, 1>& xi;
      const CORE::LINALG::Matrix<dim, 1>& n;
      const CORE::LINALG::Matrix<dim, 1>& dir;
    };
  }  // namespace DETAILS

  /*!
   * @brief Evaluate optional function @p GetCauchyNDirAndDerivativesAtXi(....) on SolidVariant @p
   * variant. If the function does not exist, an error will be thrown.
   *
   * @return CauchyNDirAndLinearization<dim>
   */
  template <int dim, typename VariantType>
  CauchyNDirAndLinearization<dim> GetCauchyNDirAndDerivativesAtXi(VariantType& variant,
      const DRT::Element& element, MAT::So3Material& mat, const std::vector<double>& disp,
      const CORE::LINALG::Matrix<dim, 1>& xi, const CORE::LINALG::Matrix<dim, 1>& n,
      const CORE::LINALG::Matrix<dim, 1>& dir)
  {
    return std::visit(
        DETAILS::EvaluateCauchyNDirAction<dim>(element, mat, disp, xi, n, dir), variant);
  }
}  // namespace DRT::ELEMENTS
BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_3D_ELE_CALC_LIB_NITSCHE_HPP