/*! \file

\brief A displacement based solid element formulation with MULF prestressing

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_CALC_MULF_HPP
#define BACI_SOLID_3D_ELE_CALC_MULF_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_lib_element.hpp"
#include "baci_solid_3D_ele_calc.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"
#include "baci_solid_3D_ele_calc_lib_io.hpp"
#include "baci_solid_3D_ele_calc_lib_mulf.H"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype>
  struct MulfLinearizationContainer
  {
    CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        Bop{};
  };

  /*!
   * @brief A displacement based solid element formulation with MULF prestressing
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct MulfFormulation
  {
    static constexpr bool has_gauss_point_history = true;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;
    static constexpr bool is_prestress_updatable = true;

    using LinearizationContainer = MulfLinearizationContainer<celltype>;
    using GaussPointHistory = MulfHistoryData<celltype>;

    template <typename Evaluator>
    static auto Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
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

      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchygreen =
          EvaluateCauchyGreen<celltype>(spatial_material_mapping);

      const CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1> gl_strain =
          EvaluateGreenLagrangeStrain(cauchygreen);

      CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
          CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
          Bop = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

      const MulfLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            return MulfLinearizationContainer<celltype>{
                EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping)};
          });

      return evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static inline SolidFormulationLinearization<celltype> EvaluateFullLinearization(
        const DRT::Element& ele, const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        MulfHistoryData<celltype>& history_data)
    {
      dserror(
          "The full linearization is not yet implemented for the displacement based formulation "
          "with MULF prestressing.");
    }

    static CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    GetLinearBOperator(const MulfLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void AddInternalForceVector(const MulfLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        MulfHistoryData<celltype>& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::AddInternalForceVector(
          linearization.Bop, stress, integration_factor, force_vector);
    }

    static void AddStiffnessMatrix(const MulfLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, MulfHistoryData<celltype>& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(
          linearization.Bop, stress, integration_factor, stiffness_matrix);
      DRT::ELEMENTS::AddGeometricStiffnessMatrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }

    static void Pack(const MulfHistoryData<celltype>& history_data, CORE::COMM::PackBuffer& data)
    {
      CORE::COMM::ParObject::AddtoPack(data, history_data.inverse_jacobian);
      CORE::COMM::ParObject::AddtoPack(data, history_data.deformation_gradient);
      CORE::COMM::ParObject::AddtoPack(data, history_data.is_setup);
    }

    static void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        MulfHistoryData<celltype>& history_data)
    {
      CORE::COMM::ParObject::ExtractfromPack(position, data, history_data.inverse_jacobian);
      CORE::COMM::ParObject::ExtractfromPack(position, data, history_data.deformation_gradient);
      int is_setup_int;
      CORE::COMM::ParObject::ExtractfromPack(position, data, is_setup_int);
      history_data.is_setup = static_cast<bool>(is_setup_int);
    }

    static inline void UpdatePrestress(const DRT::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        MulfHistoryData<celltype>& history_data)
    {
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> delta_defgrd =
          EvaluateMulfDeformationGradientUpdate(
              shape_functions, element_nodes.displacements_, history_data);

      // update mulf history data only if prestress is active
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> inv_delta_defgrd(
          delta_defgrd);
      inv_delta_defgrd.Invert();


      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invJ_new;

      invJ_new.MultiplyTN(inv_delta_defgrd, history_data.inverse_jacobian);

      history_data.deformation_gradient = deformation_gradient;
      history_data.inverse_jacobian = std::move(invJ_new);
    }
  };

  template <CORE::FE::CellType celltype>
  using MulfSolidIntegrator = SolidEleCalc<celltype, MulfFormulation<celltype>>;

  template <typename T, typename AlwaysVoid = void>
  constexpr bool IsPrestressUpdateable = false;

  template <typename T>
  constexpr bool IsPrestressUpdateable<T,
      std::void_t<decltype(std::declval<T>()->UpdatePrestress(std::declval<const DRT::Element&>(),
          std::declval<MAT::So3Material&>(), std::declval<const DRT::Discretization&>(),
          std::declval<const std::vector<int>&>(), std::declval<Teuchos::ParameterList&>()))>> =
      true;

  namespace DETAILS
  {
    struct UpdatePrestressAction
    {
      UpdatePrestressAction(const DRT::Element& e, MAT::So3Material& m,
          const DRT::Discretization& d, const std::vector<int>& lmvec, Teuchos::ParameterList& p)
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
        dserror(
            "Your element evaluation %s does not allow to update prestress. You may need to add "
            "MULF to your element line definitions.",
            CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      }

      const DRT::Element& element;
      MAT::So3Material& mat;
      const DRT::Discretization& discretization;
      const std::vector<int>& lm;
      Teuchos::ParameterList& params;
    };
  }  // namespace DETAILS

  template <typename VariantType>
  void UpdatePrestress(VariantType& variant, const DRT::Element& element, MAT::So3Material& mat,
      const DRT::Discretization& discretization, const std::vector<int>& lm,
      Teuchos::ParameterList& params)
  {
    std::visit(DETAILS::UpdatePrestressAction(element, mat, discretization, lm, params), variant);
  }
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_3D_ELE_CALC_MULF_HPP
