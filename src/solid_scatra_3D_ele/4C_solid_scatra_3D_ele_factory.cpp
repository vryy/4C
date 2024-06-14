/*! \file

\brief Factory of solid-scatra elements

\level 1
*/


#include "4C_solid_scatra_3D_ele_factory.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
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
  template <Core::FE::CellType celltype, Inpar::STR::KinemType kinem,
      Discret::ELEMENTS::ElementTechnology ele_tech,
      Discret::ELEMENTS::PrestressTechnology prestress_technology, typename Enable = void>
  struct SolidScatraCalculationFormulation
  {
  };

  /*!
   * @brief Standard nonlinear displacement based total lagrangian formulation valid for all
   * celltypes
   */
  template <Core::FE::CellType celltype>
  struct SolidScatraCalculationFormulation<celltype, Inpar::STR::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::none, Discret::ELEMENTS::PrestressTechnology::none>
  {
    using type = Discret::ELEMENTS::Details::DisplacementBasedSolidScatraIntegrator<celltype>;
  };

  /*!
   * @brief Nonlinear total lagrangian formulation with F-Bar for hex8 and pyramid 5
   */
  template <Core::FE::CellType celltype>
  struct SolidScatraCalculationFormulation<celltype, Inpar::STR::KinemType::nonlinearTotLag,
      Discret::ELEMENTS::ElementTechnology::fbar, Discret::ELEMENTS::PrestressTechnology::none,
      std::enable_if_t<celltype == Core::FE::CellType::hex8>>
  {
    using type = Discret::ELEMENTS::Details::FBarSolidScatraIntegrator<celltype>;
  };
}  // namespace

void Discret::ELEMENTS::add_to_pack(Core::Communication::PackBuffer& data,
    const Discret::ELEMENTS::SolidScatraElementProperties& properties)
{
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(properties.impltype));

  Discret::ELEMENTS::add_to_pack(data, properties.solid);
}

void Discret::ELEMENTS::extract_from_pack(std::size_t& position, const std::vector<char>& data,
    Discret::ELEMENTS::SolidScatraElementProperties& properties)
{
  properties.impltype = static_cast<Inpar::ScaTra::ImplType>(
      Core::Communication::ParObject::extract_int(position, data));

  Discret::ELEMENTS::ExtractFromPack(position, data, properties.solid);
}

Discret::ELEMENTS::SolidScatraCalcVariant
Discret::ELEMENTS::create_solid_scatra_calculation_interface(Core::FE::CellType celltype,
    const Discret::ELEMENTS::SolidElementProperties& element_properties)
{
  // We have 4 different element properties and each combination results in a different element
  // formulation.
  return Core::FE::CellTypeSwitch<Details::ImplementedSolidScatraCellTypes>(celltype,
      [&](auto celltype_t)
      {
        return switch_kinematic_type(element_properties.kintype,
            [&](auto kinemtype_t)
            {
              return element_technology_switch(element_properties.element_technology,
                  [&](auto eletech_t)
                  {
                    return prestress_technology_switch(element_properties.prestress_technology,
                        [&](auto prestress_tech_t) -> SolidScatraCalcVariant
                        {
                          constexpr Core::FE::CellType celltype_c = celltype_t();
                          constexpr Inpar::STR::KinemType kinemtype_c = kinemtype_t();
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
                              Core::FE::celltype_string<celltype_t()>,
                              Inpar::STR::KinemTypeString(element_properties.kintype).c_str(),
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