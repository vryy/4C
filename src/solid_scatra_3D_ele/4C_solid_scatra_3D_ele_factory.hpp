/*! \file

\brief Factory of solid-scatra elements

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"
#include "4C_solid_scatra_3D_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{

  /*!
   *  @brief struct for managing solidscatra element properties
   */
  struct SolidScatraElementProperties
  {
    SolidElementProperties solid{};

    //! scalar transport implementation type (physics)
    Inpar::ScaTra::ImplType impltype{Inpar::ScaTra::ImplType::impltype_undefined};
  };


  void add_to_pack(Core::Communication::PackBuffer& data,
      const Discret::ELEMENTS::SolidScatraElementProperties& properties);

  void extract_from_pack(std::size_t& position, const std::vector<char>& data,
      Discret::ELEMENTS::SolidScatraElementProperties& properties);

  namespace Details
  {

    using ImplementedSolidScatraCellTypes = Core::FE::CelltypeSequence<Core::FE::CellType::hex8,
        Core::FE::CellType::hex27, Core::FE::CellType::tet4, Core::FE::CellType::tet10>;

    // Displacement based integrators
    template <Core::FE::CellType celltype>
    using DisplacementBasedSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedFormulation<celltype>>;
    using DisplacementBasedSolidScatraEvaluator =
        Core::FE::apply_celltype_sequence<DisplacementBasedSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes>;

    // FBar evaluators
    template <Core::FE::CellType celltype>
    using FBarSolidScatraIntegrator = SolidScatraEleCalc<celltype, FBarFormulation<celltype>>;
    using FbarScatraEvaluators = Core::FE::apply_celltype_sequence<FBarSolidScatraIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8>>;

    using SolidScatraEvaluators =
        Core::FE::Join<DisplacementBasedSolidScatraEvaluator, FbarScatraEvaluators>;
  }  // namespace Details

  /// Variant holding the different implementations for the solid-scatra element
  using SolidScatraCalcVariant = CreateVariantType<Details::SolidScatraEvaluators>;

  SolidScatraCalcVariant create_solid_scatra_calculation_interface(Core::FE::CellType celltype,
      const Discret::ELEMENTS::SolidElementProperties& element_properties);
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif