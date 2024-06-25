/*! \file

\brief A displacement based solid element formulation with MULF prestressing

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_MULF_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_MULF_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_calc_lib_mulf.H"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype>
  struct MulfLinearizationContainer
  {
    Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
        Bop{};
  };

  /*!
   * @brief A displacement based solid element formulation with MULF prestressing
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct MulfFormulation
  {
    static constexpr bool has_gauss_point_history = true;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;
    static constexpr bool is_prestress_updatable = true;

    using LinearizationContainer = MulfLinearizationContainer<celltype>;
    using GaussPointHistory = MulfHistoryData<celltype>;

    template <typename Evaluator>
    static auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping, MulfHistoryData<celltype>& history_data,
        Evaluator evaluator)
    {
      if (!history_data.is_setup)
      {
        history_data.inverse_jacobian = jacobian_mapping.inverse_jacobian_;
        history_data.is_setup = true;
      }

      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          evaluate_mulf_spatial_material_mapping(
              jacobian_mapping, shape_functions, element_nodes.displacements_, history_data);

      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> cauchygreen =
          evaluate_cauchy_green<celltype>(spatial_material_mapping);

      const Core::LinAlg::Matrix<Details::num_str<celltype>, 1> gl_strain =
          evaluate_green_lagrange_strain(cauchygreen);

      Core::LinAlg::Matrix<Details::num_str<celltype>,
          Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
          Bop = evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

      const MulfLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            return MulfLinearizationContainer<celltype>{
                evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping)};
          });

      return evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    get_linear_b_operator(const MulfLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void add_internal_force_vector(const MulfLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        MulfHistoryData<celltype>& history_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::ELEMENTS::add_internal_force_vector(
          linearization.Bop, stress, integration_factor, force_vector);
    }

    static void add_stiffness_matrix(const MulfLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, MulfHistoryData<celltype>& history_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
            Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::ELEMENTS::add_elastic_stiffness_matrix(
          linearization.Bop, stress, integration_factor, stiffness_matrix);
      Discret::ELEMENTS::add_geometric_stiffness_matrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }

    static void pack(
        const MulfHistoryData<celltype>& history_data, Core::Communication::PackBuffer& data)
    {
      Core::Communication::ParObject::add_to_pack(data, history_data.inverse_jacobian);
      Core::Communication::ParObject::add_to_pack(data, history_data.deformation_gradient);
      Core::Communication::ParObject::add_to_pack(data, history_data.is_setup);
    }

    static void unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        MulfHistoryData<celltype>& history_data)
    {
      Core::Communication::ParObject::extract_from_pack(
          position, data, history_data.inverse_jacobian);
      Core::Communication::ParObject::extract_from_pack(
          position, data, history_data.deformation_gradient);
      int is_setup_int;
      Core::Communication::ParObject::extract_from_pack(position, data, is_setup_int);
      history_data.is_setup = static_cast<bool>(is_setup_int);
    }

    static inline void update_prestress(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        MulfHistoryData<celltype>& history_data)
    {
      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> delta_defgrd =
          evaluate_mulf_deformation_gradient_update(
              shape_functions, element_nodes.displacements_, history_data);

      // update mulf history data only if prestress is active
      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> inv_delta_defgrd(
          delta_defgrd);
      inv_delta_defgrd.invert();


      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> invJ_new;

      invJ_new.multiply_tn(inv_delta_defgrd, history_data.inverse_jacobian);

      history_data.deformation_gradient = deformation_gradient;
      history_data.inverse_jacobian = std::move(invJ_new);
    }
  };

  template <Core::FE::CellType celltype>
  using MulfSolidIntegrator = SolidEleCalc<celltype, MulfFormulation<celltype>>;

  template <typename T, typename AlwaysVoid = void>
  constexpr bool IsPrestressUpdateable = false;

  template <typename T>
  constexpr bool IsPrestressUpdateable<T,
      std::void_t<decltype(std::declval<T>()->update_prestress(
          std::declval<const Core::Elements::Element&>(), std::declval<Mat::So3Material&>(),
          std::declval<const Core::FE::Discretization&>(), std::declval<const std::vector<int>&>(),
          std::declval<Teuchos::ParameterList&>()))>> = true;

  namespace Details
  {
    struct UpdatePrestressAction
    {
      UpdatePrestressAction(const Core::Elements::Element& e, Mat::So3Material& m,
          const Core::FE::Discretization& d, const std::vector<int>& lmvec,
          Teuchos::ParameterList& p)
          : element(e), mat(m), discretization(d), lm(lmvec), params(p)
      {
      }

      template <typename T, std::enable_if_t<IsPrestressUpdateable<T&>, bool> = true>
      void operator()(T& updateable)
      {
        updateable->update_prestress(element, mat, discretization, lm, params);
      }

      template <typename T, std::enable_if_t<!IsPrestressUpdateable<T&>, bool> = true>
      void operator()(T& other)
      {
        FOUR_C_THROW(
            "Your element evaluation %s does not allow to update prestress. You may need to add "
            "MULF to your element line definitions.",
            Core::UTILS::TryDemangle(typeid(T).name()).c_str());
      }

      const Core::Elements::Element& element;
      Mat::So3Material& mat;
      const Core::FE::Discretization& discretization;
      const std::vector<int>& lm;
      Teuchos::ParameterList& params;
    };
  }  // namespace Details

  template <typename VariantType>
  void update_prestress(VariantType& variant, const Core::Elements::Element& element,
      Mat::So3Material& mat, const Core::FE::Discretization& discretization,
      const std::vector<int>& lm, Teuchos::ParameterList& params)
  {
    std::visit(Details::UpdatePrestressAction(element, mat, discretization, lm, params), variant);
  }
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
