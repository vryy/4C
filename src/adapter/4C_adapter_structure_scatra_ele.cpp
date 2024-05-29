/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_adapter_structure_scatra_ele.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_membrane_scatra.hpp"
#include "4C_shell7p_ele_scatra.hpp"
#include "4C_so3_scatra.hpp"
#include "4C_solid_scatra_3D_ele.hpp"
#include "4C_truss3_scatra.hpp"
#include "4C_w1_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::SCATRA::ImplType ADAPTER::GetScaTraImplType(CORE::Elements::Element* ele)
{
  INPAR::SCATRA::ImplType impltype(INPAR::SCATRA::impltype_undefined);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = ele->ElementType().Name();

  // tet 4 solid scatra
  if (eletypename == "So_tet4ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(
             ele))
            ->ImplType();
  }
  // tet10 solid scatra
  else if (eletypename == "So_tet10ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>*>(
             ele))
            ->ImplType();
  }
  // HEX 8 Elements
  // hex8 solid scatra
  else if (eletypename == "So_hex8ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
             ele))
            ->ImplType();
  }
  // hex8fbar solid scatra
  else if (eletypename == "So_hex8fbarScatraType")
  {
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>*>(ele))
            ->ImplType();
  }
  // hex27 solid scatra
  else if (eletypename == "So_hex27ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>*>(
             ele))
            ->ImplType();
  }
  // wedge6
  else if (eletypename == "So_weg6ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoWeg6, CORE::FE::CellType::wedge6>*>(
             ele))
            ->ImplType();
  }
  // wall scatra elements
  else if (eletypename == "Wall1ScatraType")
  {
    impltype = (dynamic_cast<DRT::ELEMENTS::Wall1Scatra*>(ele))->ImplType();
  }
  // shell scatra elements
  else if (eletypename == "Shell7pScatraType")
  {
    impltype = (dynamic_cast<DRT::ELEMENTS::Shell7pScatra*>(ele))->ImplType();
  }
  // membrane3 scatra element
  else if (eletypename == "MembraneScatra_tri3Type")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>*>(ele))->ImplType();
  }
  // membrane6 scatra element
  else if (eletypename == "MembraneScatra_tri6Type")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>*>(ele))->ImplType();
  }
  // membrane4 scatra element
  else if (eletypename == "MembraneScatra_quad4Type")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>*>(ele))->ImplType();
  }
  // membrane9 scatra element
  else if (eletypename == "MembraneScatra_quad9Type")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>*>(ele))->ImplType();
  }
  // truss3 scatra element
  else if (eletypename == "Truss3ScatraType")
  {
    impltype = (dynamic_cast<DRT::ELEMENTS::Truss3Scatra*>(ele))->ImplType();
  }
  // SolidScatra element
  else if (eletypename == "SolidScatraType")
  {
    impltype = (dynamic_cast<DRT::ELEMENTS::SolidScatra*>(ele))->ImplType();
  }
  else
  {
    if (!(eletypename == "Bele3Type")) return impltype;
    impltype = INPAR::SCATRA::impltype_no_physics;
  }

  return impltype;
}
FOUR_C_NAMESPACE_CLOSE
