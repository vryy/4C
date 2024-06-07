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
    static auto Evaluate(const Core::Elements::Element& ele,
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
          EvaluateMulfSpatialMaterialMapping(
              jacobian_mapping, shape_functions, element_nodes.displacements_, history_data);

      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> cauchygreen =
          EvaluateCauchyGreen<celltype>(spatial_material_mapping);

      const Core::LinAlg::Matrix<Details::num_str<celltype>, 1> gl_strain =
          EvaluateGreenLagrangeStrain(cauchygreen);

      Core::LinAlg::Matrix<Details::num_str<celltype>,
          Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
          Bop = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

      const MulfLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            return MulfLinearizationContainer<celltype>{
                EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping)};
          });

      return evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    GetLinearBOperator(const MulfLinearizationContainer<celltype>& linearization)
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

    static void AddStiffnessMatrix(const MulfLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, MulfHistoryData<celltype>& history_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
            Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::ELEMENTS::AddElasticStiffnessMatrix(
          linearization.Bop, stress, integration_factor, stiffness_matrix);
      Discret::ELEMENTS::AddGeometricStiffnessMatrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }

    static void Pack(
        const MulfHistoryData<celltype>& history_data, Core::Communication::PackBuffer& data)
    {
      Core::Communication::ParObject::AddtoPack(data, history_data.inverse_jacobian);
      Core::Communication::ParObject::AddtoPack(data, history_data.deformation_gradient);
      Core::Communication::ParObject::AddtoPack(data, history_data.is_setup);
    }

    static void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        MulfHistoryData<celltype>& history_data)
    {
      Core::Communication::ParObject::ExtractfromPack(
          position, data, history_data.inverse_jacobian);
      Core::Communication::ParObject::ExtractfromPack(
          position, data, history_data.deformation_gradient);
      int is_setup_int;
      Core::Communication::ParObject::ExtractfromPack(position, data, is_setup_int);
      history_data.is_setup = static_cast<bool>(is_setup_int);
    }

    static inline void UpdatePrestress(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        MulfHistoryData<celltype>& history_data)
    {
      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> delta_defgrd =
          EvaluateMulfDeformationGradientUpdate(
              shape_functions, element_nodes.displacements_, history_data);

      // update mulf history data only if prestress is active
      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> inv_delta_defgrd(
          delta_defgrd);
      inv_delta_defgrd.Invert();


      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> invJ_new;

      invJ_new.MultiplyTN(inv_delta_defgrd, history_data.inverse_jacobian);

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
      std::void_t<decltype(std::declval<T>()->UpdatePrestress(
          std::declval<const Core::Elements::Element&>(), std::declval<Mat::So3Material&>(),
          std::declval<const Discret::Discretization&>(), std::declval<const std::vector<int>&>(),
          std::declval<Teuchos::ParameterList&>()))>> = true;

  namespace Details
  {
    struct UpdatePrestressAction
    {
      UpdatePrestressAction(const Core::Elements::Element& e, Mat::So3Material& m,
          const Discret::Discretization& d, const std::vector<int>& lmvec,
          Teuchos::ParameterList& p)
          : element(e), mat(m), discretization(d), lm(lmvec), params(p)
      {
      }

      template <typename T, std::enable_if_t<IsPrestressUpdateable<T&>, bool> = true>
      void operator()(T& updateable)
      {
        updateable->UpdatePrestress(element, mat, discretization, lm, params);
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
      const Discret::Discretization& discretization;
      const std::vector<int>& lm;
      Teuchos::ParameterList& params;
    };
  }  // namespace Details

  template <typename VariantType>
  void UpdatePrestress(VariantType& variant, const Core::Elements::Element& element,
      Mat::So3Material& mat, const Discret::Discretization& discretization,
      const std::vector<int>& lm, Teuchos::ParameterList& params)
  {
    std::visit(Details::UpdatePrestressAction(element, mat, discretization, lm, params), variant);
  }
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
