/*! \file

\brief Factory of solid elements

\level 1
*/

#include "4C_solid_3D_ele_factory.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_eas.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_calc_mulf.hpp"
#include "4C_solid_3D_ele_calc_mulf_fbar.hpp"
#include "4C_solid_3D_ele_properties.hpp"
#include "4C_utils_exceptions.hpp"

#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief A template class that is taking different element formulation switches as template
   * parameter. If implemented, the struct defines the type of the solid evaluation.
   *
   * @tparam celltype
   * @tparam kinem : Kinematic type (linear, nonliner)
   * @tparam ele_tech : Element technology (none, FBar, EAS mild and full)
   * @tparam prestress_technology : Prestress technology (none or mulf)
   * @tparam Enable : A dummy parameter for enabling a subset of switches.
   */
  template <CORE::FE::CellType celltype, INPAR::STR::KinemType kinem,
      DRT::ELEMENTS::ElementTechnology ele_tech,
      DRT::ELEMENTS::PrestressTechnology prestress_technology, typename Enable = void>
  struct SolidCalculationFormulation
  {
  };

  /*!
   * @brief Standard nonlinear displacement based total lagrangian formulation valid for all
   * celltypes
   */
  template <CORE::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::none, DRT::ELEMENTS::PrestressTechnology::none>
  {
    using type = DRT::ELEMENTS::DisplacementBasedSolidIntegrator<celltype>;
  };

  /*!
   * @brief Displacement based formulation valid for all celltypes with linear kinematics (i.e.
   * small displacements)
   */
  template <CORE::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, INPAR::STR::KinemType::linear,
      DRT::ELEMENTS::ElementTechnology::none, DRT::ELEMENTS::PrestressTechnology::none>
  {
    using type = DRT::ELEMENTS::DisplacementBasedLinearKinematicsSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with mild EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<CORE::FE::CellType::hex8,
      INPAR::STR::KinemType::nonlinearTotLag, DRT::ELEMENTS::ElementTechnology::eas_mild,
      DRT::ELEMENTS::PrestressTechnology::none>
  {
    using type = DRT::ELEMENTS::SolidEleCalcEas<CORE::FE::CellType::hex8,
        STR::ELEMENTS::EasType::eastype_h8_9>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with full EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<CORE::FE::CellType::hex8,
      INPAR::STR::KinemType::nonlinearTotLag, DRT::ELEMENTS::ElementTechnology::eas_full,
      DRT::ELEMENTS::PrestressTechnology::none>
  {
    using type = DRT::ELEMENTS::SolidEleCalcEas<CORE::FE::CellType::hex8,
        STR::ELEMENTS::EasType::eastype_h8_21>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with F-Bar for hex8 and pyramid 5
   */
  template <CORE::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::fbar, DRT::ELEMENTS::PrestressTechnology::none,
      std::enable_if_t<celltype == CORE::FE::CellType::hex8 ||
                       celltype == CORE::FE::CellType::pyramid5>>
  {
    using type = DRT::ELEMENTS::FBarSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear formulation with F-Bar and MULF prestressing for hex8 and pyramid5
   */
  template <CORE::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::fbar, DRT::ELEMENTS::PrestressTechnology::mulf,
      std::enable_if_t<celltype == CORE::FE::CellType::hex8 ||
                       celltype == CORE::FE::CellType::pyramid5>>
  {
    using type = DRT::ELEMENTS::MulfFBarSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear modified updated lagrangian prestressing for all celltypes
   */
  template <CORE::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::none, DRT::ELEMENTS::PrestressTechnology::mulf>
  {
    using type = DRT::ELEMENTS::SolidEleCalc<celltype, DRT::ELEMENTS::MulfFormulation<celltype>>;
  };
}  // namespace

DRT::ELEMENTS::SolidCalcVariant DRT::ELEMENTS::CreateSolidCalculationInterface(
    CORE::FE::CellType celltype, const DRT::ELEMENTS::SolidElementProperties& element_properties)
{
  // We have 4 different element properties and each combination results in a different element
  // formulation.
  return CORE::FE::CellTypeSwitch<DETAILS::ImplementedSolidCellTypes>(celltype,
      [&](auto celltype_t)
      {
        return switch_kinematic_type(element_properties.kintype,
            [&](auto kinemtype_t)
            {
              return ElementTechnologySwitch(element_properties.element_technology,
                  [&](auto eletech_t)
                  {
                    return PrestressTechnologySwitch(element_properties.prestress_technology,
                        [&](auto prestress_tech_t) -> SolidCalcVariant
                        {
                          constexpr CORE::FE::CellType celltype_c = celltype_t();
                          constexpr INPAR::STR::KinemType kinemtype_c = kinemtype_t();
                          constexpr ElementTechnology eletech_c = eletech_t();
                          constexpr PrestressTechnology prestress_tech_c = prestress_tech_t();
                          if constexpr (is_valid_type<SolidCalculationFormulation<celltype_c,
                                            kinemtype_c, eletech_c, prestress_tech_c>>)
                          {
                            return typename SolidCalculationFormulation<celltype_c, kinemtype_c,
                                eletech_c, prestress_tech_c>::type();
                          }

                          FOUR_C_THROW(
                              "Your element formulation with cell type %s, kinematic type %s,"
                              " elememt technology %s and prestress type %s oes not exist ",
                              CORE::FE::celltype_string<celltype_t()>,
                              INPAR::STR::KinemTypeString(element_properties.kintype).c_str(),
                              ElementTechnologyString(element_properties.element_technology)
                                  .c_str(),
                              PrestressTechnologyString(element_properties.prestress_technology)
                                  .c_str());
                        });
                  });
            });
      });
}

FOUR_C_NAMESPACE_CLOSE
