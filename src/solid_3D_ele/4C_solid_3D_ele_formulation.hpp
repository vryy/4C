/*! \file

\brief This file contains helper functions and type traits for the solid formulation

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_FORMULATION_HPP
#define FOUR_C_SOLID_3D_ELE_FORMULATION_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_utils_exceptions.hpp"

#include <tuple>
#include <type_traits>
#include <utility>


FOUR_C_NAMESPACE_OPEN


namespace Discret::ELEMENTS
{
  /*!
   * @brief a trait for the solid formulation determining whether the formulation has a gauss point
   * history
   *
   * @tparam SolidFormulation
   */
  template <typename SolidFormulation>
  constexpr bool has_gauss_point_history = SolidFormulation::has_gauss_point_history;

  /*!
   * @brief a trait for the solid formulation determining whether the formulation has a global
   * history
   *
   * @tparam SolidFormulation
   */
  template <typename SolidFormulation>
  constexpr bool has_global_history = SolidFormulation::has_global_history;

  /*!
   * @brief a trait for the solid formulation determining whether the formulation has preparation
   * data (i.e. data executed once per element)
   *
   * @tparam SolidFormulation
   */
  template <typename SolidFormulation>
  constexpr bool has_preparation_data = SolidFormulation::has_preparation_data;

  namespace Details
  {
    template <typename SolidFormulation, typename = void>
    struct IsPrestressUpdatable : std::false_type
    {
    };

    template <typename SolidFormulation>
    struct IsPrestressUpdatable<SolidFormulation,
        std::enable_if_t<SolidFormulation::is_prestress_updatable, void>> : std::true_type
    {
    };
  }  // namespace Details

  template <typename SolidFormulation>
  constexpr bool is_prestress_updatable = Details::IsPrestressUpdatable<SolidFormulation>::value;

  namespace Details
  {
    /*!
     * @brief A dummy type that is used if a solid formulation does not need preparation data
     */
    struct NoneType
    {
    };

    template <typename SolidFormulation, typename T = void>
    struct PreparationTypeTrait;

    template <typename SolidFormulation>
    struct PreparationTypeTrait<SolidFormulation,
        std::enable_if_t<has_preparation_data<SolidFormulation>, void>>
    {
      using type = typename SolidFormulation::PreparationData;
    };

    template <typename SolidFormulation>
    struct PreparationTypeTrait<SolidFormulation,
        std::enable_if_t<!has_preparation_data<SolidFormulation>, void>>
    {
      using type = NoneType;
    };
  }  // namespace Details

  /*!
   * @brief A type trait that is the type of the Preparation data if SolidFormulation needs to
   * prepare data, otherwise, it is @p Details::NoneType
   *
   * @tparam SolidFormulation
   */
  template <typename SolidFormulation>
  using PreparationData = typename Details::PreparationTypeTrait<SolidFormulation>::type;

  /*!
   * @brief An object holding the history data of the solid formulation.
   *
   * Can hold none, gauss point history and global history.
   *
   * @note The data is only stored if needed by the solid formulation
   *
   * @tparam SolidFormulation
   * @tparam T
   */
  template <typename SolidFormulation, typename T = void>
  struct SolidFormulationHistory;

  template <typename SolidFormulation>
  struct SolidFormulationHistory<SolidFormulation,
      std::enable_if_t<
          has_global_history<SolidFormulation> && has_gauss_point_history<SolidFormulation>, void>>
  {
    typename SolidFormulation::GlobalHistory global_history;
    std::vector<typename SolidFormulation::GaussPointHistory> gp_history;
  };

  template <typename SolidFormulation>
  struct SolidFormulationHistory<SolidFormulation,
      std::enable_if_t<
          !has_global_history<SolidFormulation> && has_gauss_point_history<SolidFormulation>, void>>
  {
    std::vector<typename SolidFormulation::GaussPointHistory> gp_history;
  };

  template <typename SolidFormulation>
  struct SolidFormulationHistory<SolidFormulation,
      std::enable_if_t<
          has_global_history<SolidFormulation> && !has_gauss_point_history<SolidFormulation>, void>>
  {
    typename SolidFormulation::GlobalHistory global_history;
  };

  template <typename SolidFormulation>
  struct SolidFormulationHistory<SolidFormulation,
      std::enable_if_t<!has_global_history<SolidFormulation> &&
                           !has_gauss_point_history<SolidFormulation>,
          void>>
  {
  };

  template <typename SolidFormulation>
  void ResizeGPHistory(SolidFormulationHistory<SolidFormulation>& history_data,
      [[maybe_unused]] const std::size_t num_gps)
  {
    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      history_data.gp_history.resize(num_gps);
    }
  }

  namespace Details
  {
    template <typename SolidFormulation>
    auto GetAdditionalPreparationTuple(const PreparationData<SolidFormulation>& preparation_data)
    {
      if constexpr (has_preparation_data<SolidFormulation>)
        return std::tie(preparation_data);
      else
        return std::tie();
    }

    template <typename SolidFormulation>
    auto GetAdditionalGlobalHistoryTuple(SolidFormulationHistory<SolidFormulation>& history_data)
    {
      if constexpr (has_global_history<SolidFormulation>)
        return std::tie(history_data.global_history);
      else
        return std::tie();
    }

    template <typename SolidFormulation>
    auto GetAdditionalGaussPointHistoryTuple(
        SolidFormulationHistory<SolidFormulation>& history_data, [[maybe_unused]] const int gp)
    {
      if constexpr (has_gauss_point_history<SolidFormulation>)
        return std::tie(history_data.gp_history[gp]);
      else
        return std::tie();
    }

    template <typename SolidFormulation>
    auto GetAdditionalGaussPointHistoryTuple(
        SolidFormulationHistory<SolidFormulation>& history_data)
    {
      static_assert(!has_gauss_point_history<SolidFormulation>,
          "The solid formulation has a Gauss point history and, therefore, needs the Gauss point "
          "id!");

      return std::tie();
    }


    template <typename SolidFormulation>
    auto GetAdditionalTuple(const PreparationData<SolidFormulation>& preparation_data,
        SolidFormulationHistory<SolidFormulation>& history_data, const int gp)
    {
      return std::tuple_cat(GetAdditionalPreparationTuple<SolidFormulation>(preparation_data),
          GetAdditionalGlobalHistoryTuple<SolidFormulation>(history_data),
          GetAdditionalGaussPointHistoryTuple<SolidFormulation>(history_data, gp));
    }

    template <typename SolidFormulation>
    auto GetAdditionalTuple(const PreparationData<SolidFormulation>& preparation_data,
        SolidFormulationHistory<SolidFormulation>& history_data)
    {
      return std::tuple_cat(GetAdditionalPreparationTuple<SolidFormulation>(preparation_data),
          GetAdditionalGlobalHistoryTuple<SolidFormulation>(history_data),
          GetAdditionalGaussPointHistoryTuple<SolidFormulation>(history_data));
    }
  }  // namespace Details

  /*!
   * @brief Pack the solid formulation history data
   *
   * @note Calls the respective @p Pack(...) method for the SolidFormulation if needed by the solid
   * formulation
   *
   * @tparam SolidFormulation
   * @param data (out) : Buffer where the data is packed to
   * @param solid_formulation_history (in) : History data to be packed
   */
  template <typename SolidFormulation>
  void Pack(Core::Communication::PackBuffer& data,
      const SolidFormulationHistory<SolidFormulation>& solid_formulation_history)
  {
    if constexpr (has_global_history<SolidFormulation>)
    {
      SolidFormulation::Pack(solid_formulation_history.global_history, data);
    }

    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      data.add_to_pack(solid_formulation_history.gp_history.size());
      for (const auto& item : solid_formulation_history.gp_history)
      {
        SolidFormulation::Pack(item, data);
      }
    }
  }

  /*!
   * @brief Unpack the solid formulation history data
   *
   * @note Calls the respective @p Unpack(...) method for the SolidFormulation if needed by the
   * solid formulation
   *
   * @tparam SolidFormulation
   * @param position (in/out) : Position to start to unpack (will be incremented)
   * @param data (in) : Data to unpack from
   * @param solid_formulation_history (out) : History data to be unpacked
   */
  template <typename SolidFormulation>
  void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
      SolidFormulationHistory<SolidFormulation>& solid_formulation_history)
  {
    if constexpr (has_global_history<SolidFormulation>)
    {
      SolidFormulation::Unpack(position, data, solid_formulation_history.global_history);
    }

    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      std::size_t num_gps;
      Core::Communication::ParObject::extract_from_pack(position, data, num_gps);
      solid_formulation_history.gp_history.resize(num_gps);
      for (auto& item : solid_formulation_history.gp_history)
      {
        SolidFormulation::Unpack(position, data, item);
      }
    }
  }

  /*!
   * @brief Calls the @p Prepare(...) on the solid formulation of the solid formulation needs to
   * prepare stuff
   *
   * @tparam SolidFormulation
   * @tparam celltype
   * @param ele (in) : Solid element
   * @param nodal_coordinates (in) : element coordinates
   * @param history_data (in/out) : history data
   * @return PreparationData<SolidFormulation>
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  PreparationData<SolidFormulation> Prepare(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& nodal_coordinates,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (has_preparation_data<SolidFormulation>)
    {
      if constexpr (has_global_history<SolidFormulation>)
        return SolidFormulation::Prepare(ele, nodal_coordinates, history_data.global_history);
      else
      {
        return SolidFormulation::Prepare(ele, nodal_coordinates);
      }
    }
    else
    {
      // nothing to prepare
      return {};
    }
  }

  /*!
   * @brief Evaluate the deformation gradient and Green-Lagrange strain tensor for the solid element
   * formulation.
   *
   * @note: This method does not support solid formulation with Gauss point history since the method
   * does not necessarily be called on a Gauss point. If called on a Gauss point, prefer the other
   * @p evaluate() call with the Gauss point id as parameter.
   */
  template <typename SolidFormulation, Core::FE::CellType celltype, typename Evaluator>
  inline auto evaluate(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, Evaluator&& evaluator)
  {
    return std::apply([](auto&&... args)
        { return SolidFormulation::evaluate(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(ele, element_nodes, xi, shape_functions, jacobian_mapping),
            Details::GetAdditionalTuple<SolidFormulation>(preparation_data, history_data),
            std::forward_as_tuple(evaluator)));
  }


  /*!
   * @brief Evaluate the deformation gradient and Green-Lagrange strain tensor for the solid element
   * formulation.
   *
   * @note: This method should be preferred if called at a Gauss point since it also supports
   * element formulations with gauss point history.
   */
  template <typename SolidFormulation, Core::FE::CellType celltype, typename Evaluator>
  inline auto evaluate(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp, Evaluator&& evaluator)
  {
    return std::apply([](auto&&... args)
        { return SolidFormulation::evaluate(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(ele, element_nodes, xi, shape_functions, jacobian_mapping),
            Details::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(evaluator)));
  }

  /*!
   * @brief Evaluate the derivative of the deformation gradient w.r.t. the displacements
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  inline Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
  evaluate_d_deformation_gradient_d_displacements(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
          deformation_gradient,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      FOUR_C_THROW("The Solid element formulation can only be evaluated at the Gauss points.");
    }
    else
    {
      return std::apply(
          [](auto&&... args)
          {
            return SolidFormulation::evaluate_d_deformation_gradient_d_displacements(
                std::forward<decltype(args)>(args)...);
          },
          std::tuple_cat(std::forward_as_tuple(ele, element_nodes, xi, shape_functions,
                             jacobian_mapping, deformation_gradient),
              Details::GetAdditionalTuple(preparation_data, history_data)));
    }
  }

  /*!
   * @brief Evaluate the derivative of the deformation gradient w.r.t. the xi
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  inline Core::LinAlg::Matrix<9, Core::FE::dim<celltype>> evaluate_d_deformation_gradient_d_xi(
      const Core::Elements::Element& ele, const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
          deformation_gradient,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      FOUR_C_THROW("The Solid element formulation can only be evaluated at the Gauss points.");
    }
    else
    {
      return std::apply(
          [](auto&&... args)
          {
            return SolidFormulation::evaluate_d_deformation_gradient_d_xi(
                std::forward<decltype(args)>(args)...);
          },
          std::tuple_cat(std::forward_as_tuple(ele, element_nodes, xi, shape_functions,
                             jacobian_mapping, deformation_gradient),
              Details::GetAdditionalTuple(preparation_data, history_data)));
    }
  }

  /*!
   * @brief Evaluate the second derivative of the deformation gradient w.r.t. the xi and the
   * displacements
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  inline Core::LinAlg::Matrix<9,
      Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
  evaluate_d_deformation_gradient_d_displacements_d_xi(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
          deformation_gradient,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      FOUR_C_THROW("The Solid element formulation can only be evaluated at the Gauss points.");
    }
    else
    {
      return std::apply(
          [](auto&&... args)
          {
            return SolidFormulation::evaluate_d_deformation_gradient_d_displacements_d_xi(
                std::forward<decltype(args)>(args)...);
          },
          std::tuple_cat(std::forward_as_tuple(ele, element_nodes, xi, shape_functions,
                             jacobian_mapping, deformation_gradient),
              Details::GetAdditionalTuple(preparation_data, history_data)));
    }
  }

  /*!
   * @brief Evaluate the linear B-Operator
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  static Core::LinAlg::Matrix<Core::FE::dim<celltype>*(Core::FE::dim<celltype> + 1) / 2,
      Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
  get_linear_b_operator(const typename SolidFormulation::LinearizationContainer& linearization)
  {
    return SolidFormulation::get_linear_b_operator(linearization);
  }

  /*!
   * @brief Add the internal force vector contribution of the Gauss point to @p force_vector
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  static inline void add_internal_force_vector(
      const typename SolidFormulation::LinearizationContainer& linearization,
      const Stress<celltype>& stress, const double integration_factor,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp,
      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1>&
          force_vector)
  {
    std::apply([](auto&&... args)
        { SolidFormulation::add_internal_force_vector(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(std::forward_as_tuple(linearization, stress, integration_factor),
            Details::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(force_vector)));
  }

  /*!
   * @brief Add stiffness matrix contribution of the Gauss point to @p stiffness_matrix
   */
  template <typename SolidFormulation, Core::FE::CellType celltype>
  static inline void add_stiffness_matrix(
      const typename SolidFormulation::LinearizationContainer& linearization,
      const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
      const double integration_factor, const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp,
      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
          Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
  {
    std::apply([](auto&&... args)
        { SolidFormulation::add_stiffness_matrix(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(linearization, jacobian_mapping, stress, integration_factor),
            Details::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(stiffness_matrix)));
  }

  template <typename SolidFormulation, Core::FE::CellType celltype>
  static inline void update_prestress(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (is_prestress_updatable<SolidFormulation>)
    {
      if constexpr (has_global_history<SolidFormulation>)
      {
        std::apply([](auto&&... args)
            { SolidFormulation::update_prestress(std::forward<decltype(args)>(args)...); },
            std::tuple_cat(std::forward_as_tuple(ele, element_nodes),
                Details::GetAdditionalPreparationTuple<SolidFormulation>(preparation_data),
                Details::GetAdditionalGlobalHistoryTuple<SolidFormulation>(history_data)));
      }
    }
    else
    {
      FOUR_C_THROW("The solid formulation does not support to update the prestress!");
    }
  }

  template <typename SolidFormulation, Core::FE::CellType celltype>
  static inline void update_prestress(const Core::Elements::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
          deformation_gradient,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, [[maybe_unused]] const int gp)
  {
    if constexpr (is_prestress_updatable<SolidFormulation>)
    {
      std::apply([](auto&&... args)
          { SolidFormulation::update_prestress(std::forward<decltype(args)>(args)...); },
          std::tuple_cat(std::forward_as_tuple(ele, element_nodes, xi, shape_functions,
                             jacobian_mapping, deformation_gradient),
              Details::GetAdditionalTuple(preparation_data, history_data, gp)));
    }
    else
    {
      FOUR_C_THROW("The solid formulation does not support to update the prestress!");
    }
  }


}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
