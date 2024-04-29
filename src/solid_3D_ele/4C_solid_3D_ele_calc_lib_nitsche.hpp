/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_NITSCHE_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_NITSCHE_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_lib_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_utils_demangle.hpp"


FOUR_C_NAMESPACE_OPEN

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

      out_mat.MultiplyTN(mat1, mat2);
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
    d2_cauchyndir_dF2_d_F_dd.Multiply(d2_cauchyndir_dF2, d_F_dd);
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
    d_cauchyndir_dxi.MultiplyTN(d_F_dxi, d_cauchyndir_dF);
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
  struct CauchyNDirLinearizations
  {
    /// @brief first derivatives
    ///@{

    /// first derivative w.r.t. normal
    CORE::LINALG::Matrix<dim, 1>* d_cauchyndir_dn = nullptr;

    /// first derivative w.r.t. given direction
    CORE::LINALG::Matrix<dim, 1>* d_cauchyndir_ddir = nullptr;

    /// first derivative w.r.t. xi
    CORE::LINALG::Matrix<dim, 1>* d_cauchyndir_dxi = nullptr;

    /// first derivative w.r.t. displacements
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd = nullptr;
    ///@}

    /// @brief second derivatives
    ///@{

    /// second derivative w.r.t. displacements
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd2 = nullptr;

    /// second derivative w.r.t. displacements, normal
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dn = nullptr;

    /// second derivative w.r.t. displacements, given direction
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_ddir = nullptr;

    /// second derivative w.r.t. displacements, xi
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dxi = nullptr;
    ///@}
  };

  /// Check wheter a solid variant can evaluate the Cauchy stress at xi in a specifc direction
  /// including the derivatives
  template <typename T, int dim, typename AlwaysVoid = void>
  constexpr bool CanEvaluateEvaluateCauchyNDir = false;

  template <typename T, int dim>
  constexpr bool CanEvaluateEvaluateCauchyNDir<T, dim,
      std::void_t<decltype(std::declval<T>()->GetCauchyNDirAtXi(std::declval<const DRT::Element&>(),
          std::declval<MAT::So3Material&>(), std::declval<const std::vector<double>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<CauchyNDirLinearizations<dim>&>()))>> = true;

  namespace DETAILS
  {
    template <int dim>
    struct EvaluateCauchyNDirAction
    {
      EvaluateCauchyNDirAction(const DRT::Element& e, MAT::So3Material& m,
          const std::vector<double>& d, const CORE::LINALG::Matrix<dim, 1>& x,
          const CORE::LINALG::Matrix<dim, 1>& normal, const CORE::LINALG::Matrix<dim, 1>& direction,
          CauchyNDirLinearizations<dim>& lins)
          : element(e), mat(m), disp(d), xi(x), n(normal), dir(direction), linearizations(lins)
      {
      }

      template <typename T, std::enable_if_t<CanEvaluateEvaluateCauchyNDir<T&, dim>, bool> = true>
      double operator()(T& cauchy_n_dir_evaluatable)
      {
        return cauchy_n_dir_evaluatable->GetCauchyNDirAtXi(
            element, mat, disp, xi, n, dir, linearizations);
      }

      template <typename T, std::enable_if_t<!CanEvaluateEvaluateCauchyNDir<T&, dim>, bool> = true>
      double operator()(T& other)
      {
        FOUR_C_THROW(
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
      CauchyNDirLinearizations<dim>& linearizations;
    };
  }  // namespace DETAILS

  /*!
   * @brief Evaluate optional function @p GetCauchyNDirAtXi(....) on SolidVariant @p
   * variant. If the function does not exist, an error will be thrown.
   *
   * @return double
   */
  template <int dim, typename VariantType>
  double GetCauchyNDirAtXi(VariantType& variant, const DRT::Element& element, MAT::So3Material& mat,
      const std::vector<double>& disp, const CORE::LINALG::Matrix<dim, 1>& xi,
      const CORE::LINALG::Matrix<dim, 1>& n, const CORE::LINALG::Matrix<dim, 1>& dir,
      CauchyNDirLinearizations<dim>& linearizations)
  {
    return std::visit(
        DETAILS::EvaluateCauchyNDirAction<dim>(element, mat, disp, xi, n, dir, linearizations),
        variant);
  }
}  // namespace DRT::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif