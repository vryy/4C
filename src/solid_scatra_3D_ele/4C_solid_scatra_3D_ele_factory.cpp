/*! \file

\brief Factory of solid-scatra elements

\level 1
*/


#include "4C_solid_scatra_3D_ele_factory.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief A template class that is taking different element formulation switches as template
   * parameter. If implemented, the struct defines the type of the solid formulation in the
   * solid-scatra element evaluation.
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
  struct SolidScatraCalculationFormulation
  {
  };

  /*!
   * @brief Standard nonlinear displacement based total lagrangian formulation valid for all
   * celltypes
   */
  template <CORE::FE::CellType celltype>
  struct SolidScatraCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::none, DRT::ELEMENTS::PrestressTechnology::none>
  {
    using type = DRT::ELEMENTS::DETAILS::DisplacementBasedSolidScatraIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with F-Bar for hex8 and pyramid 5
   */
  template <CORE::FE::CellType celltype>
  struct SolidScatraCalculationFormulation<celltype, INPAR::STR::KinemType::nonlinearTotLag,
      DRT::ELEMENTS::ElementTechnology::fbar, DRT::ELEMENTS::PrestressTechnology::none,
      std::enable_if_t<celltype == CORE::FE::CellType::hex8>>
  {
    using type = DRT::ELEMENTS::DETAILS::FBarSolidScatraIntegrator<celltype>;
  };
}  // namespace

void DRT::ELEMENTS::AddToPack(
    CORE::COMM::PackBuffer& data, const DRT::ELEMENTS::SolidScatraElementProperties& properties)
{
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(properties.impltype));

  AddToPack(data, properties.solid);
}

void DRT::ELEMENTS::ExtractFromPack(std::size_t& position, const std::vector<char>& data,
    DRT::ELEMENTS::SolidScatraElementProperties& properties)
{
  properties.impltype =
      static_cast<INPAR::SCATRA::ImplType>(CORE::COMM::ParObject::ExtractInt(position, data));

  ExtractFromPack(position, data, properties.solid);
}

DRT::ELEMENTS::SolidScatraCalcVariant DRT::ELEMENTS::CreateSolidScatraCalculationInterface(
    CORE::FE::CellType celltype, const DRT::ELEMENTS::SolidElementProperties& element_properties)
{
  // We have 4 different element properties and each combination results in a different element
  // formulation.
  return CORE::FE::CellTypeSwitch<DETAILS::ImplementedSolidScatraCellTypes>(celltype,
      [&](auto celltype_t)
      {
        return switch_kinematic_type(element_properties.kintype,
            [&](auto kinemtype_t)
            {
              return ElementTechnologySwitch(element_properties.element_technology,
                  [&](auto eletech_t)
                  {
                    return PrestressTechnologySwitch(element_properties.prestress_technology,
                        [&](auto prestress_tech_t) -> SolidScatraCalcVariant
                        {
                          constexpr CORE::FE::CellType celltype_c = celltype_t();
                          constexpr INPAR::STR::KinemType kinemtype_c = kinemtype_t();
                          constexpr ElementTechnology eletech_c = eletech_t();
                          constexpr PrestressTechnology prestress_tech_c = prestress_tech_t();
                          if constexpr (is_valid_type<SolidScatraCalculationFormulation<celltype_c,
                                            kinemtype_c, eletech_c, prestress_tech_c>>)
                          {
                            return typename SolidScatraCalculationFormulation<celltype_c,
                                kinemtype_c, eletech_c, prestress_tech_c>::type();
                          }

                          FOUR_C_THROW(
                              "Your element formulation with cell type %s, kinematic type %s,"
                              " elememt technology %s and prestress type %s oes not exist in the "
                              "solid-scatra context.",
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