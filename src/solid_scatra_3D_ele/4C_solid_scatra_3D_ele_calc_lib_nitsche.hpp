/*! \file

\brief A library of free functions for a solid-scatra element with Nitsche contact

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_CALC_LIB_NITSCHE_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_CALC_LIB_NITSCHE_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele_calc_lib_nitsche.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{

  template <int dim>
  struct SolidScatraCauchyNDirLinearizations
  {
    /// all pure solid linearizations
    CauchyNDirLinearizations<dim> solid{};

    /// first derivative w.r.t. scalars
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_ds = nullptr;
  };

  template <CORE::FE::CellType celltype, typename SolidFormulation>
  CauchyNDirLinearizationDependencies<celltype>
  get_initialized_cauchy_n_dir_linearization_dependencies(
      const DRT::ELEMENTS::ElementFormulationDerivativeEvaluator<celltype, SolidFormulation>&
          evaluator,
      DRT::ELEMENTS::SolidScatraCauchyNDirLinearizations<3>& linearizations)
  {
    // Get pure solid dependencies
    CauchyNDirLinearizationDependencies<celltype> linearization_dependencies =
        DRT::ELEMENTS::get_initialized_cauchy_n_dir_linearization_dependencies(
            evaluator, linearizations.solid);

    // initialize dependencies for solid-scatra
    if (linearizations.d_cauchyndir_ds && !linearization_dependencies.d_cauchyndir_dF.has_value())
    {
      linearization_dependencies.d_cauchyndir_dF.emplace();
    }

    return linearization_dependencies;
  }
  /// Check wheter a solid-scatra variant can evaluate the Cauchy stress at xi in a specifc
  /// direction including the derivatives
  template <typename T, int dim, typename AlwaysVoid = void>
  constexpr bool can_evaluate_solid_scatra_cauchy_n_dir_at_xi = false;

  template <typename T, int dim>
  constexpr bool can_evaluate_solid_scatra_cauchy_n_dir_at_xi<T, dim,
      std::void_t<decltype(std::declval<T>()->GetCauchyNDirAtXi(std::declval<const DRT::Element&>(),
          std::declval<MAT::So3Material&>(), std::declval<const std::vector<double>&>(),
          std::declval<const std::optional<std::vector<double>>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<const CORE::LINALG::Matrix<dim, 1>&>(),
          std::declval<SolidScatraCauchyNDirLinearizations<dim>&>()))>> = true;


  namespace DETAILS
  {
    template <int dim>
    struct EvaluateSolidScatraCauchyNDirAction
    {
      EvaluateSolidScatraCauchyNDirAction(const DRT::Element& e, MAT::So3Material& m,
          const std::vector<double>& d, const std::optional<std::vector<double>>& s,
          const CORE::LINALG::Matrix<dim, 1>& x, const CORE::LINALG::Matrix<dim, 1>& normal,
          const CORE::LINALG::Matrix<dim, 1>& direction,
          SolidScatraCauchyNDirLinearizations<dim>& lins)
          : element(e),
            mat(m),
            disp(d),
            scalars(s),
            xi(x),
            n(normal),
            dir(direction),
            linearizations(lins)
      {
      }

      template <typename T,
          std::enable_if_t<can_evaluate_solid_scatra_cauchy_n_dir_at_xi<T&, dim>, bool> = true>
      double operator()(T& cauchy_n_dir_evaluatable)
      {
        return cauchy_n_dir_evaluatable->GetCauchyNDirAtXi(
            element, mat, disp, scalars, xi, n, dir, linearizations);
      }

      template <typename T,
          std::enable_if_t<!can_evaluate_solid_scatra_cauchy_n_dir_at_xi<T&, dim>, bool> = true>
      double operator()(T& other)
      {
        FOUR_C_THROW(
            "Your element evaluation %s does not allow to evaluate the Cauchy stress at a "
            "specific "
            "point in a specific direction in the dimension dim=%d.",
            CORE::UTILS::TryDemangle(typeid(T).name()).c_str(), dim);
      }

      const DRT::Element& element;
      MAT::So3Material& mat;
      const std::vector<double>& disp;
      const std::optional<std::vector<double>>& scalars;
      const CORE::LINALG::Matrix<dim, 1>& xi;
      const CORE::LINALG::Matrix<dim, 1>& n;
      const CORE::LINALG::Matrix<dim, 1>& dir;
      SolidScatraCauchyNDirLinearizations<dim>& linearizations;
    };
  }  // namespace DETAILS

  template <typename VariantType>
  double GetCauchyNDirAtXi(VariantType& variant, const DRT::Element& element, MAT::So3Material& mat,
      const std::vector<double>& disp, const std::optional<std::vector<double>>& scalars,
      const CORE::LINALG::Matrix<3, 1>& xi, const CORE::LINALG::Matrix<3, 1>& n,
      const CORE::LINALG::Matrix<3, 1>& dir, SolidScatraCauchyNDirLinearizations<3>& linearizations)
  {
    return std::visit(DETAILS::EvaluateSolidScatraCauchyNDirAction<3>(
                          element, mat, disp, scalars, xi, n, dir, linearizations),
        variant);
  }
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif