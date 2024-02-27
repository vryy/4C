/*! \file

\brief This file contains helper functions and type traits for the solid formulation

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_FORMULATION_HPP
#define BACI_SOLID_3D_ELE_FORMULATION_HPP

#include "baci_config.hpp"

#include "baci_comm_pack_buffer.hpp"
#include "baci_discretization_fem_general_cell_type.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"
#include "baci_utils_exceptions.hpp"

#include <tuple>
#include <type_traits>
#include <utility>


BACI_NAMESPACE_OPEN


namespace DRT::ELEMENTS
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


  namespace DETAILS
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
  }  // namespace DETAILS

  /*!
   * @brief A type trait that is the type of the Preparation data if SolidFormulation needs to
   * prepare data, otherwise, it is @p DETAILS::NoneType
   *
   * @tparam SolidFormulation
   */
  template <typename SolidFormulation>
  using PreparationData = typename DETAILS::PreparationTypeTrait<SolidFormulation>::type;

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

  namespace DETAILS
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
  }  // namespace DETAILS

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
  void Pack(CORE::COMM::PackBuffer& data,
      const SolidFormulationHistory<SolidFormulation>& solid_formulation_history)
  {
    if constexpr (has_global_history<SolidFormulation>)
    {
      SolidFormulation::Pack(solid_formulation_history.global_history, data);
    }

    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      data.AddtoPack(solid_formulation_history.gp_history.size());
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
      CORE::COMM::ParObject::ExtractfromPack(position, data, num_gps);
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
  template <typename SolidFormulation, CORE::FE::CellType celltype>
  PreparationData<SolidFormulation> Prepare(const DRT::Element& ele,
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
   * @p Evaluate() call with the Gauss point id as parameter.
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype, typename Evaluator>
  inline auto Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& element_nodes,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, Evaluator&& evaluator)
  {
    return std::apply([](auto&&... args)
        { return SolidFormulation::Evaluate(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(ele, element_nodes, xi, shape_functions, jacobian_mapping),
            DETAILS::GetAdditionalTuple<SolidFormulation>(preparation_data, history_data),
            std::forward_as_tuple(evaluator)));
  }


  /*!
   * @brief Evaluate the deformation gradient and Green-Lagrange strain tensor for the solid element
   * formulation.
   *
   * @note: This method should be preferred if called at a Gauss point since it also supports
   * element formulations with gauss point history.
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype, typename Evaluator>
  inline auto Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& element_nodes,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp, Evaluator&& evaluator)
  {
    return std::apply([](auto&&... args)
        { return SolidFormulation::Evaluate(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(ele, element_nodes, xi, shape_functions, jacobian_mapping),
            DETAILS::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(evaluator)));
  }

  /*!
   * @brief Evaluate the full linearization of the solid formulation
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype>
  inline SolidFormulationLinearization<celltype> EvaluateFullLinearization(const DRT::Element& ele,
      const ElementNodes<celltype>& element_nodes,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
          deformation_gradient,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data)
  {
    if constexpr (has_gauss_point_history<SolidFormulation>)
    {
      dserror("The Solid element formulation can only be evaluated at the Gauss points.");
    }
    else
    {
      return std::apply(
          [](auto&&... args) {
            return SolidFormulation::EvaluateFullLinearization(
                std::forward<decltype(args)>(args)...);
          },
          std::tuple_cat(std::forward_as_tuple(ele, element_nodes, xi, shape_functions,
                             jacobian_mapping, deformation_gradient),
              DETAILS::GetAdditionalTuple(preparation_data, history_data)));
    }
  }

  /*!
   * @brief Evaluate the linear B-Operator
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype>
  static CORE::LINALG::Matrix<CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2,
      CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
  GetLinearBOperator(const typename SolidFormulation::LinearizationContainer& linearization)
  {
    return SolidFormulation::GetLinearBOperator(linearization);
  }

  /*!
   * @brief Add the internal force vector contribution of the Gauss point to @p force_vector
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype>
  static inline void AddInternalForceVector(
      const typename SolidFormulation::LinearizationContainer& linearization,
      const Stress<celltype>& stress, const double integration_factor,
      const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp,
      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
          force_vector)
  {
    std::apply([](auto&&... args)
        { SolidFormulation::AddInternalForceVector(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(std::forward_as_tuple(linearization, stress, integration_factor),
            DETAILS::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(force_vector)));
  }

  /*!
   * @brief Add stiffness matrix contribution of the Gauss point to @p stiffness_matrix
   */
  template <typename SolidFormulation, CORE::FE::CellType celltype>
  static inline void AddStiffnessMatrix(
      const typename SolidFormulation::LinearizationContainer& linearization,
      const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
      const double integration_factor, const PreparationData<SolidFormulation>& preparation_data,
      SolidFormulationHistory<SolidFormulation>& history_data, const int gp,
      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
          CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
  {
    std::apply([](auto&&... args)
        { SolidFormulation::AddStiffnessMatrix(std::forward<decltype(args)>(args)...); },
        std::tuple_cat(
            std::forward_as_tuple(linearization, jacobian_mapping, stress, integration_factor),
            DETAILS::GetAdditionalTuple(preparation_data, history_data, gp),
            std::forward_as_tuple(stiffness_matrix)));
  }

}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_3D_ELE_FORMULATION_HPP
