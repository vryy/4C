// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"
#include "4C_solid_scatra_3D_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
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
      const Discret::Elements::SolidScatraElementProperties& properties);

  void extract_from_pack(Core::Communication::UnpackBuffer& buffer,
      Discret::Elements::SolidScatraElementProperties& properties);

  namespace Internal
  {

    using ImplementedSolidScatraCellTypes =
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::hex27,
            Core::FE::CellType::tet4, Core::FE::CellType::tet10, Core::FE::CellType::nurbs27>;

    // Displacement based integrators
    template <Core::FE::CellType celltype>
    using DisplacementBasedSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedFormulation<celltype>>;
    using DisplacementBasedSolidScatraEvaluator =
        Core::FE::apply_celltype_sequence<DisplacementBasedSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes>;

    // Displacement based integrators with linear kinematics
    template <Core::FE::CellType celltype>
    using DisplacementBasedLinearKinematicsSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedLinearKinematicsFormulation<celltype>>;
    using DisplacementBasedLinearKinematicsSolidScatraEvaluator =
        Core::FE::apply_celltype_sequence<DisplacementBasedLinearKinematicsSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes>;


    // FBar evaluators
    template <Core::FE::CellType celltype>
    using FBarSolidScatraIntegrator = SolidScatraEleCalc<celltype, FBarFormulation<celltype>>;
    using FbarScatraEvaluators = Core::FE::apply_celltype_sequence<FBarSolidScatraIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8>>;

    using SolidScatraEvaluators = Core::FE::Join<DisplacementBasedSolidScatraEvaluator,
        DisplacementBasedLinearKinematicsSolidScatraEvaluator, FbarScatraEvaluators>;
  }  // namespace Internal

  /// Variant holding the different implementations for the solid-scatra element
  using SolidScatraCalcVariant = CreateVariantType<Internal::SolidScatraEvaluators>;

  SolidScatraCalcVariant create_solid_scatra_calculation_interface(Core::FE::CellType celltype,
      const Discret::Elements::SolidElementProperties& element_properties);
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif