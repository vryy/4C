/*! \file

\brief Factory of solid-scatra elements

\level 1
*/


#include "baci_solid_scatra_ele_factory.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"

BACI_NAMESPACE_OPEN

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

BACI_NAMESPACE_CLOSE