// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_3D_ele_calc_eas.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_3D_ele_calc_mulf.hpp"
#include "4C_solid_3D_ele_calc_mulf_fbar.hpp"
#include "4C_solid_3D_ele_calc_shell_ans.hpp"
#include "4C_solid_3D_ele_calc_shell_eas_ans.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid::Elements
{
  enum class EasType;
}
namespace Discret::Elements
{
  using ImplementedSolidCellTypes3D = Core::FE::CelltypeSequence<Core::FE::CellType::hex8,
      Core::FE::CellType::hex18, Core::FE::CellType::hex20, Core::FE::CellType::hex27,
      Core::FE::CellType::nurbs27, Core::FE::CellType::tet4, Core::FE::CellType::tet10,
      Core::FE::CellType::wedge6, Core::FE::CellType::pyramid5>;

  using ImplementedSolidCellTypes2D = Core::FE::CelltypeSequence<Core::FE::CellType::quad4,
      Core::FE::CellType::quad8, Core::FE::CellType::quad9, Core::FE::CellType::tri3,
      Core::FE::CellType::tri6, Core::FE::CellType::nurbs9>;

  template <unsigned dim>
  using ImplementedSolidCellTypes =
      std::conditional_t<dim == 3, ImplementedSolidCellTypes3D, ImplementedSolidCellTypes2D>;
  namespace Internal
  {
    template <unsigned dim>
    using DisplacementBasedEvaluators =
        Core::FE::apply_celltype_sequence<DisplacementBasedSolidIntegrator,
            ImplementedSolidCellTypes<dim>>;

    template <unsigned dim>
    using DisplacementBasedLinearKinematicsEvaluators =
        Core::FE::apply_celltype_sequence<DisplacementBasedLinearKinematicsSolidIntegrator,
            ImplementedSolidCellTypes<dim>>;

    using FbarEvaluators = Core::FE::apply_celltype_sequence<FBarSolidIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>>;

    template <unsigned dim>
    using EASEvaluators = std::conditional_t<dim == 3,
        Core::FE::BaseTypeList<
            EASSolidIntegrator<Core::FE::CellType::hex8, Discret::Elements::EasType::eastype_h8_9,
                Inpar::Solid::KinemType::nonlinearTotLag>,
            EASSolidIntegrator<Core::FE::CellType::hex8, Discret::Elements::EasType::eastype_h8_21,
                Inpar::Solid::KinemType::nonlinearTotLag>,
            EASSolidIntegrator<Core::FE::CellType::hex8, Discret::Elements::EasType::eastype_sh8_7,
                Inpar::Solid::KinemType::nonlinearTotLag>,
            EASSolidIntegrator<Core::FE::CellType::hex8, Discret::Elements::EasType::eastype_h8_9,
                Inpar::Solid::KinemType::linear>,
            EASSolidIntegrator<Core::FE::CellType::hex8, Discret::Elements::EasType::eastype_h8_21,
                Inpar::Solid::KinemType::linear>>,
        Core::FE::BaseTypeList<EASSolidIntegrator<Core::FE::CellType::quad4,
            Discret::Elements::EasType::eastype_q4_4, Inpar::Solid::KinemType::nonlinearTotLag>>>;

    using MulfEvaluators =
        Core::FE::apply_celltype_sequence<MulfSolidIntegrator, ImplementedSolidCellTypes<3>>;
    using FBarMulfEvaluators = Core::FE::apply_celltype_sequence<MulfFBarSolidIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>>;

    using SolidShellEvaluators = Core::FE::apply_celltype_sequence<ANSSolidShellIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::wedge6>>;

    using SolidShellEasEvaluators =
        Core::FE::BaseTypeList<EasAnsSolidShellIntegrator<Core::FE::CellType::hex8,
                                   Discret::Elements::EasType::eastype_sh8_7>,
            EasAnsSolidShellIntegrator<Core::FE::CellType::wedge6,
                Discret::Elements::EasType::eastype_sw6_1>>;

    template <unsigned dim>
    using SolidEvaluators = std::conditional_t<dim == 3,
        Core::FE::Join<DisplacementBasedEvaluators<dim>,
            DisplacementBasedLinearKinematicsEvaluators<dim>, FbarEvaluators, EASEvaluators<dim>,
            MulfEvaluators, FBarMulfEvaluators, SolidShellEvaluators, SolidShellEasEvaluators>,
        Core::FE::Join<DisplacementBasedEvaluators<dim>,
            DisplacementBasedLinearKinematicsEvaluators<dim>, EASEvaluators<dim>>>;
  }  // namespace Internal

  template <unsigned dim>
  using SolidCalcVariant = CreateVariantType<Internal::SolidEvaluators<dim>>;

  template <unsigned dim>
  SolidCalcVariant<dim> create_solid_calculation_interface(Core::FE::CellType celltype,
      const Discret::Elements::SolidElementProperties<dim>& element_properties,
      const SolidIntegrationRules<dim>& integration_rules);

}  // namespace Discret::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
