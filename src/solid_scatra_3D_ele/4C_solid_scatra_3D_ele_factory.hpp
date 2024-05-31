/*! \file

\brief Factory of solid-scatra elements

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"
#include "4C_solid_scatra_3D_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{

  /*!
   *  @brief struct for managing solidscatra element properties
   */
  struct SolidScatraElementProperties
  {
    SolidElementProperties solid{};

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype{INPAR::SCATRA::ImplType::impltype_undefined};
  };


  void AddToPack(
      CORE::COMM::PackBuffer& data, const DRT::ELEMENTS::SolidScatraElementProperties& properties);

  void ExtractFromPack(std::size_t& position, const std::vector<char>& data,
      DRT::ELEMENTS::SolidScatraElementProperties& properties);

  namespace DETAILS
  {

    using ImplementedSolidScatraCellTypes = CORE::FE::CelltypeSequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10>;

    // Displacement based integrators
    template <CORE::FE::CellType celltype>
    using DisplacementBasedSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedFormulation<celltype>>;
    using DisplacementBasedSolidScatraEvaluator =
        CORE::FE::apply_celltype_sequence<DisplacementBasedSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes>;

    // FBar evaluators
    template <CORE::FE::CellType celltype>
    using FBarSolidScatraIntegrator = SolidScatraEleCalc<celltype, FBarFormulation<celltype>>;
    using FbarScatraEvaluators = CORE::FE::apply_celltype_sequence<FBarSolidScatraIntegrator,
        CORE::FE::CelltypeSequence<CORE::FE::CellType::hex8>>;

    using SolidScatraEvaluators =
        CORE::FE::Join<DisplacementBasedSolidScatraEvaluator, FbarScatraEvaluators>;
  }  // namespace DETAILS

  /// Variant holding the different implementations for the solid-scatra element
  using SolidScatraCalcVariant = CreateVariantType<DETAILS::SolidScatraEvaluators>;

  SolidScatraCalcVariant CreateSolidScatraCalculationInterface(
      CORE::FE::CellType celltype, const DRT::ELEMENTS::SolidElementProperties& element_properties);
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif