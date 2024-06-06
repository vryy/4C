/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_poro.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_poro_3D_ele_calc_pressure_based.hpp"

#include <memory>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  /*!
   *  @brief struct for managing solidporo element properties
   */
  struct SolidPoroElementProperties
  {
    //! porosity implementation type (physics)
    Inpar::Poro::PoroType porotype{Inpar::Poro::PoroType::undefined};

    //! scalar transport implementation type (physics)
    Inpar::ScaTra::ImplType impltype{Inpar::ScaTra::ImplType::impltype_undefined};
  };

  namespace Details
  {
    using ImplementedSolidPoroCellTypes = Core::FE::CelltypeSequence<Core::FE::CellType::hex8,
        Core::FE::CellType::hex27, Core::FE::CellType::tet4, Core::FE::CellType::tet10>;
    using PoroPressureBasedEvaluators =
        Core::FE::apply_celltype_sequence<SolidPoroPressureBasedEleCalc,
            ImplementedSolidPoroCellTypes>;

    using SolidPoroEvaluators = Core::FE::Join<PoroPressureBasedEvaluators>;
  }  // namespace Details
  using SolidPoroCalcVariant = CreateVariantType<Details::SolidPoroEvaluators>;

  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoro;


  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(
      Core::Elements::Element& ele, Inpar::Poro::PoroType porotype);

  template <Core::FE::CellType celltype>
  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(Inpar::Poro::PoroType porotype);
}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif