/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_inpar_poro.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_lib_element.hpp"
#include "baci_solid_3D_ele_factory_lib.hpp"
#include "baci_solid_poro_3D_ele_calc_pressure_based.hpp"

#include <memory>
#include <variant>

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /*!
   *  @brief struct for managing solidporo element properties
   */
  struct SolidPoroElementProperties
  {
    //! porosity implementation type (physics)
    INPAR::PORO::PoroType porotype{INPAR::PORO::PoroType::undefined};

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype{INPAR::SCATRA::ImplType::impltype_undefined};
  };

  namespace DETAILS
  {
    using ImplementedSolidPoroCellTypes = CORE::FE::celltype_sequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10>;
    using PoroPressureBasedEvaluators =
        CORE::FE::apply_celltype_sequence<SolidPoroPressureBasedEleCalc,
            ImplementedSolidPoroCellTypes>;

    using SolidPoroEvaluators = CORE::FE::Join<PoroPressureBasedEvaluators>;
  }  // namespace DETAILS
  using SolidPoroCalcVariant = CreateVariantType<DETAILS::SolidPoroEvaluators>;

  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoro;


  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(
      DRT::Element& ele, INPAR::PORO::PoroType porotype);

  template <CORE::FE::CellType celltype>
  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(INPAR::PORO::PoroType porotype);
}  // namespace DRT::ELEMENTS


BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_PORO_3D_ELE_FACTORY_H