/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid for scalar transport coupling problems

\level 2

 *----------------------------------------------------------------------*/

#include "4C_poroelast_scatra_utils_clonestrategy.hpp"

#include "4C_adapter_structure_scatra_ele.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_so3_nurbs27.hpp"
#include "4C_so3_poro_p1_scatra.hpp"
#include "4C_so3_poro_scatra.hpp"
#include "4C_so3_scatra.hpp"
#include "4C_solid_poro_3D_ele.hpp"
#include "4C_w1_poro_p1_scatra.hpp"
#include "4C_w1_poro_scatra.hpp"

FOUR_C_NAMESPACE_OPEN


INPAR::SCATRA::ImplType POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::GetImplType(
    CORE::Elements::Element* ele  //! element whose SCATRA::ImplType shall be determined
)
{
  // the element type name, needed to cast correctly in the following
  const std::string eletypename = ele->ElementType().Name();

  // TET 4 Elements
  // tet 4 solid poro scatra
  if (eletypename == "So_tet4PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(ele))
        ->ImplType();
  }
  // tet4 solid porop1 scatra
  else if (eletypename == "So_tet4PoroP1ScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(ele))
        ->ImplType();
  }
  // tet 10 solid poro scatra
  else if (eletypename == "So_tet10PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>*>(ele))
        ->ImplType();
  }
  // HEX 8 Elements
  // hex8 solid poro scatra
  else if (eletypename == "So_hex8PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(ele))
        ->ImplType();
  }
  // hex8 solid porop1 scatra
  else if (eletypename == "So_hex8PoroP1ScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(ele))
        ->ImplType();
  }
  // hex27 solid poro scatra
  else if (eletypename == "So_hex27PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>*>(ele))
        ->ImplType();
  }
  // nurbs 27
  else if (eletypename == "So_nurbs27PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27,
                CORE::FE::CellType::nurbs27>*>(ele))
        ->ImplType();
  }
  // Solidporo
  else if (eletypename == "SolidPoroType")
  {
    return (dynamic_cast<DRT::ELEMENTS::SolidPoro*>(ele))->GetImplType();
  }
  // wall poro scatra elements
  // quad 4
  else if (eletypename == "WallQuad4PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::quad4>*>(ele))
        ->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::quad9>*>(ele))
        ->ImplType();
  }
  // nurbs 4
  else if (eletypename == "WallNurbs4PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::nurbs4>*>(ele))
        ->ImplType();
  }
  // nurbs 9
  else if (eletypename == "WallNurbs9PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::nurbs9>*>(ele))
        ->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::tri3>*>(ele))
        ->ImplType();
  }
  // wall poro p1 elements
  // quad 4
  else if (eletypename == "WallQuad4PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroP1Scatra<CORE::FE::CellType::quad4>*>(ele))
        ->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroP1Scatra<CORE::FE::CellType::quad9>*>(ele))
        ->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1PoroP1Scatra<CORE::FE::CellType::tri3>*>(ele))
        ->ImplType();
  }
  // call base class routine
  else
    return ADAPTER::GetScaTraImplType(ele);
}

bool POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone elements based on poro or scatra elements
  if (POROELASTSCATRA::UTILS::IsPoroScatraElement(actele) or
      POROELAST::UTILS::IsPoroElement(actele))
  {
    // we only support fluid elements here
    eletype.emplace_back("TRANSP");
    return true;
  }

  return false;
}


void POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::set_element_data(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele, const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans == nullptr) FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());


  // set material
  trans->SetMaterial(matid, oldele);
  // set distype as well!
  trans->SetDisType(oldele->Shape());

  // now check whether ImplType is reasonable and if set the ImplType
  INPAR::SCATRA::ImplType impltype =
      POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::GetImplType(oldele);

  if (impltype == INPAR::SCATRA::impltype_undefined)
    FOUR_C_THROW(
        "PoroScatraCloneStrategy copies scatra discretization from structure discretization, but "
        "the "
        "STRUCTURE elements that are defined in the .dat file are either not meant to be copied "
        "to scatra elements "
        "or the ImplType is set 'Undefined' which is not meaningful for the created scatra "
        "discretization! "
        "Use SOLIDSCATRA, WALLSCATRA, SHELLSCATRA, SOLIDPOROSCATRA, SOLIDPOROP1SCATRA, SOLIDPORO"
        "WALLPOROSCATRA or "
        "WALLPOROP1SCATRA Elements with meaningful ImplType instead!");

  trans->SetImplType(impltype);
}

std::map<std::string, std::string>
POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::conditions_to_copy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"TransportLineNeumann", "LineNeumann"}, {"TransportSurfaceNeumann", "SurfaceNeumann"},
      {"TransportVolumeNeumann", "VolumeNeumann"},
      {"TransportNeumannInflow", "TransportNeumannInflow"}, {"FSICoupling", "FSICoupling"},
      {"ScaTraCoupling", "ScaTraCoupling"}, {"ScaTraFluxCalc", "ScaTraFluxCalc"},
      {"Initfield", "Initfield"}, {"PoroCoupling", "PoroCoupling"}, {"Initfield", "Initfield"},
      {"ArtScatraCouplConNodebased", "ArtScatraCouplConNodebased"},
      {"PoroMultiphaseScatraOxyPartPressCalcCond", "PoroMultiphaseScatraOxyPartPressCalcCond"},
      {"TransportRobin", "TransportRobin"},
      {"ArtScatraCouplConNodeToPoint", "ArtScatraCouplConNodeToPoint"}};
}

void POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  CORE::Materials::MaterialType mtype =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if ((mtype != CORE::Materials::m_scatra) && (mtype != CORE::Materials::m_elchmat) &&
      (mtype != CORE::Materials::m_electrode) && (mtype != CORE::Materials::m_matlist) &&
      (mtype != CORE::Materials::m_matlist_reactions) && (mtype != CORE::Materials::m_myocard) &&
      (mtype != CORE::Materials::m_thermostvenant))
    FOUR_C_THROW("Material with ID %d is not admissible for scalar transport elements", matid);
}

bool POROELASTSCATRA::UTILS::PoroelastCloneStrategyforScatraElements::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro or scatra element (we support submeshes here)
  if (POROELASTSCATRA::UTILS::IsPoroScatraElement(actele) or
      POROELAST::UTILS::IsPoroElement(actele))
  {
    // we only support fluid elements here
    eletype.emplace_back("FLUIDPORO");
    return true;
  }

  return false;
}
FOUR_C_NAMESPACE_CLOSE
