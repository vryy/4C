/*! \file

\brief Factory of solid-scatra elements

\level 1
*/


#include "4C_solid_scatra_3D_ele_factory.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SolidScatraCalcVariant DRT::ELEMENTS::CreateSolidScatraCalculationInterface(
    CORE::FE::CellType celltype)
{
  return CORE::FE::CellTypeSwitch<DETAILS::ImplementedSolidScatraCellTypes>(celltype,
      [&](auto celltype_t) { return CreateSolidScatraCalculationInterface<celltype_t()>(); });
}

template <CORE::FE::CellType celltype>
DRT::ELEMENTS::SolidScatraCalcVariant DRT::ELEMENTS::CreateSolidScatraCalculationInterface()
{
  return DETAILS::DisplacementBasedSolidScatraIntegrator<celltype>();
}

FOUR_C_NAMESPACE_CLOSE