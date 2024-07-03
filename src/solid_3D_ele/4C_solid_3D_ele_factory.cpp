/*! \file

\brief Factory of solid elements

\level 1
*/

#include "4C_solid_3D_ele_factory.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
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
  template <Core::FE::CellType celltype, Inpar::Solid::KinemType kinem,
      Discret::ELEMENTS::ElementTechnology ele_tech,
      Discret::ELEMENTS::PrestressTechnology prestress_technology, typename Enable = void>
  struct SolidCalculationFormulation
  {
  };

  /*!
   * @brief Standard nonlinear displacement based total lagrangian formulation valid for all
   * celltypes
   */
  template <Core::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, Inpar::Solid::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::none, Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::DisplacementBasedSolidIntegrator<celltype>;
  };

  /*!
   * @brief Displacement based formulation valid for all celltypes with linear kinematics (i.e.
   * small displacements)
   */
  template <Core::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, Inpar::Solid::KinemType::linear,
      Discret::ELEMENTS::ElementTechnology::none, Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::DisplacementBasedLinearKinematicsSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with mild EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<Core::FE::CellType::hex8,
      Inpar::Solid::KinemType::nonlinearTotLag, Discret::ELEMENTS::ElementTechnology::eas_mild,
      Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::ELEMENTS::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>;
  };

  /*!
   * @brief Small displacements formulation with mild EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<Core::FE::CellType::hex8, Inpar::Solid::KinemType::linear,
      Discret::ELEMENTS::ElementTechnology::eas_mild, Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::ELEMENTS::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with full EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<Core::FE::CellType::hex8,
      Inpar::Solid::KinemType::nonlinearTotLag, Discret::ELEMENTS::ElementTechnology::eas_full,
      Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::ELEMENTS::EasType::eastype_h8_21, Inpar::Solid::KinemType::nonlinearTotLag>;
  };

  /*!
   * @brief Small displacements formulation with full EAS for hex8
   */
  template <>
  struct SolidCalculationFormulation<Core::FE::CellType::hex8, Inpar::Solid::KinemType::linear,
      Discret::ELEMENTS::ElementTechnology::eas_full, Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::SolidEleCalcEas<Core::FE::CellType::hex8,
        Solid::ELEMENTS::EasType::eastype_h8_21, Inpar::Solid::KinemType::linear>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with F-Bar for hex8 and pyramid 5
   */
  template <Core::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, Inpar::Solid::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::fbar, Discret::ELEMENTS::PrestressTechnology::none,
      std::enable_if_t<celltype == Core::FE::CellType::hex8 ||
                       celltype == Core::FE::CellType::pyramid5>>
  {
    using type = Discret::ELEMENTS::FBarSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear formulation with F-Bar and MULF prestressing for hex8 and pyramid5
   */
  template <Core::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, Inpar::Solid::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::fbar, Discret::ELEMENTS::PrestressTechnology::mulf,
      std::enable_if_t<celltype == Core::FE::CellType::hex8 ||
                       celltype == Core::FE::CellType::pyramid5>>
  {
    using type = Discret::ELEMENTS::MulfFBarSolidIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear modified updated lagrangian prestressing for all celltypes
   */
  template <Core::FE::CellType celltype>
  struct SolidCalculationFormulation<celltype, Inpar::Solid::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::none, Discret::ELEMENTS::PrestressTechnology::mulf>
  {
    using type =
        Discret::ELEMENTS::SolidEleCalc<celltype, Discret::ELEMENTS::MulfFormulation<celltype>>;
  };
}  // namespace

Discret::ELEMENTS::SolidCalcVariant Discret::ELEMENTS::create_solid_calculation_interface(
    Core::FE::CellType celltype,
    const Discret::ELEMENTS::SolidElementProperties& element_properties)
{
  // We have 4 different element properties and each combination results in a different element
  // formulation.
  return Core::FE::CellTypeSwitch<Details::ImplementedSolidCellTypes>(celltype,
      [&](auto celltype_t)
      {
        return switch_kinematic_type(element_properties.kintype,
            [&](auto kinemtype_t)
            {
              return element_technology_switch(element_properties.element_technology,
                  [&](auto eletech_t)
                  {
                    return prestress_technology_switch(element_properties.prestress_technology,
                        [&](auto prestress_tech_t) -> SolidCalcVariant
                        {
                          constexpr Core::FE::CellType celltype_c = celltype_t();
                          constexpr Inpar::Solid::KinemType kinemtype_c = kinemtype_t();
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
                              Core::FE::celltype_string<celltype_t()>,
                              Inpar::Solid::KinemTypeString(element_properties.kintype).c_str(),
                              element_technology_string(element_properties.element_technology)
                                  .c_str(),
                              prestress_technology_string(element_properties.prestress_technology)
                                  .c_str());
                        });
                  });
            });
      });
}

FOUR_C_NAMESPACE_CLOSE
