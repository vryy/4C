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
#include "4C_solid_3D_ele_calc_eas.hpp"
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
  template <unsigned dim>
  struct SolidScatraElementProperties
  {
    SolidElementProperties<dim> solid{};

    //! scalar transport implementation type (physics)
    Inpar::ScaTra::ImplType impltype{Inpar::ScaTra::ImplType::impltype_undefined};
  };

  namespace Internal
  {

    template <unsigned dim>
    using ImplementedSolidScatraCellTypes = std::conditional_t<dim == 3,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::hex27,
            Core::FE::CellType::tet4, Core::FE::CellType::tet10, Core::FE::CellType::nurbs27>,
        Core::FE::CelltypeSequence<Core::FE::CellType::quad4, Core::FE::CellType::quad9,
            Core::FE::CellType::tri3, Core::FE::CellType::tri6>>;

    // Displacement based integrators
    template <Core::FE::CellType celltype>
    using DisplacementBasedSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedFormulation<celltype>>;

    template <unsigned dim>
    using DisplacementBasedSolidScatraEvaluator =
        Core::FE::apply_celltype_sequence<DisplacementBasedSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes<dim>>;

    // Displacement based integrators with linear kinematics
    template <Core::FE::CellType celltype>
    using DisplacementBasedLinearKinematicsSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, DisplacementBasedLinearKinematicsFormulation<celltype>>;
    template <unsigned dim>
    using DisplacementBasedLinearKinematicsSolidScatraEvaluator =
        Core::FE::apply_celltype_sequence<DisplacementBasedLinearKinematicsSolidScatraIntegrator,
            ImplementedSolidScatraCellTypes<dim>>;


    // FBar evaluators
    template <Core::FE::CellType celltype>
    using FBarSolidScatraIntegrator = SolidScatraEleCalc<celltype, FBarFormulation<celltype>>;
    using FbarScatraEvaluators = Core::FE::apply_celltype_sequence<FBarSolidScatraIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8>>;

    // Eas evaluators
    template <Core::FE::CellType celltype, Discret::Elements::EasType eas_type,
        Inpar::Solid::KinemType kinem_type>
    using EASSolidScatraIntegrator =
        SolidScatraEleCalc<celltype, EASFormulation<celltype, eas_type, kinem_type>>;
    template <unsigned dim>
    using EASScatraEvaluators = std::conditional_t<dim == 3,
        Core::FE::BaseTypeList<
            EASSolidScatraIntegrator<Core::FE::CellType::hex8,
                Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>,
            EASSolidScatraIntegrator<Core::FE::CellType::hex8,
                Discret::Elements::EasType::eastype_h8_21,
                Inpar::Solid::KinemType::nonlinearTotLag>>,
        Core::FE::BaseTypeList<EASSolidScatraIntegrator<Core::FE::CellType::quad4,
            Discret::Elements::EasType::eastype_q4_4, Inpar::Solid::KinemType::nonlinearTotLag>>>;

    template <unsigned dim>
    using SolidScatraEvaluators = std::conditional_t<dim == 3,
        Core::FE::Join<DisplacementBasedSolidScatraEvaluator<dim>,
            DisplacementBasedLinearKinematicsSolidScatraEvaluator<dim>, FbarScatraEvaluators,
            EASScatraEvaluators<dim>>,
        Core::FE::Join<DisplacementBasedSolidScatraEvaluator<dim>,
            DisplacementBasedLinearKinematicsSolidScatraEvaluator<dim>, EASScatraEvaluators<dim>>>;
  }  // namespace Internal

  /// Variant holding the different implementations for the solid-scatra element
  template <unsigned dim>
  using SolidScatraCalcVariant = CreateVariantType<Internal::SolidScatraEvaluators<dim>>;

  template <unsigned dim>
  SolidScatraCalcVariant<dim> create_solid_scatra_calculation_interface(Core::FE::CellType celltype,
      const Discret::Elements::SolidElementProperties<dim>& element_properties);
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif