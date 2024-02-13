/*! \file

\brief Factory of solid-scatra elements

\level 1
*/

#ifndef BACI_SOLID_SCATRA_ELE_FACTORY_HPP
#define BACI_SOLID_SCATRA_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_solid_ele_calc_displacement_based.hpp"
#include "baci_solid_ele_factory_lib.hpp"
#include "baci_solid_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /*!
   *  @brief struct for managing solidscatra element properties
   */
  struct SolidScatraElementProperties
  {
    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype{INPAR::SCATRA::ImplType::impltype_undefined};
  };

  namespace DETAILS
  {

    using ImplementedSolidScatraCellTypes = CORE::FE::celltype_sequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10>;

    template <CORE::FE::CellType celltype>
    using DisplacementBasedSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedFormulation<celltype>,
            DisplacementBasedPreparationData, DisplacementBasedHistoryData>;

    using DisplacementBasedSolidScatraEvaluator =
        CORE::FE::apply_celltype_sequence<DisplacementBasedSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes>;

    using SolidScatraEvaluators = CORE::FE::Join<DisplacementBasedSolidScatraEvaluator>;
  }  // namespace DETAILS

  /// Variant holding the different implementations for the solid-scatra element
  using SolidScatraCalcVariant = CreateVariantType<DETAILS::SolidScatraEvaluators>;

  SolidScatraCalcVariant CreateSolidScatraCalculationInterface(CORE::FE::CellType celltype);

  template <CORE::FE::CellType celltype>
  SolidScatraCalcVariant CreateSolidScatraCalculationInterface();
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_SCATRA_ELE_FACTORY_HPP