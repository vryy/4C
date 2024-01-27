/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid for scalar transport coupling problems

\level 2

 *----------------------------------------------------------------------*/

#include "baci_poroelast_scatra_utils_clonestrategy.H"

#include "baci_adapter_structure_scatra_ele.H"
#include "baci_global_data.H"
#include "baci_mat_fluidporo.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_structporo.H"
#include "baci_poroelast_scatra_utils.H"
#include "baci_poroelast_utils_clonestrategy.H"
#include "baci_scatra_ele.H"
#include "baci_so3_nurbs27.H"
#include "baci_so3_poro_p1_scatra.H"
#include "baci_so3_poro_scatra.H"
#include "baci_so3_scatra.H"
#include "baci_solid_poro_ele.H"
#include "baci_w1_poro_p1_scatra.H"
#include "baci_w1_poro_scatra.H"

BACI_NAMESPACE_OPEN


INPAR::SCATRA::ImplType POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::GetImplType(
    DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
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
            DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>*>(ele))
        ->ImplType();
  }
  // tet4 solid porop1 scatra
  else if (eletypename == "So_tet4PoroP1ScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>*>(
            ele))
        ->ImplType();
  }
  // tet 10 solid poro scatra
  else if (eletypename == "So_tet10PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>*>(
            ele))
        ->ImplType();
  }
  // HEX 8 Elements
  // hex8 solid poro scatra
  else if (eletypename == "So_hex8PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>*>(ele))
        ->ImplType();
  }
  // hex8 solid porop1 scatra
  else if (eletypename == "So_hex8PoroP1ScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>*>(
            ele))
        ->ImplType();
  }
  // hex27 solid poro scatra
  else if (eletypename == "So_hex27PoroScatraType")
  {
    return (
        dynamic_cast<
            DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>*>(
            ele))
        ->ImplType();
  }
  // nurbs 27
  else if (eletypename == "So_nurbs27PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
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
    return (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad4>*>(ele))
        ->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad9>*>(ele))
        ->ImplType();
  }
  // nurbs 4
  else if (eletypename == "WallNurbs4PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs4>*>(ele))
        ->ImplType();
  }
  // nurbs 9
  else if (eletypename == "WallNurbs9PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs9>*>(ele))
        ->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::tri3>*>(ele))
        ->ImplType();
  }
  // wall poro p1 elements
  // quad 4
  else if (eletypename == "WallQuad4PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<CORE::FE::CellType::quad4>*>(ele))
        ->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<CORE::FE::CellType::quad9>*>(ele))
        ->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroP1ScatraType")
  {
    return (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<CORE::FE::CellType::tri3>*>(ele))
        ->ImplType();
  }
  // call base class routine
  else
    return ADAPTER::GetScaTraImplType(ele);
}

bool POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
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


void POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans == nullptr) dserror("unsupported element type '%s'", typeid(*newele).name());


  // set material
  trans->SetMaterial(matid, oldele);
  // set distype as well!
  trans->SetDisType(oldele->Shape());

  // now check whether ImplType is reasonable and if set the ImplType
  INPAR::SCATRA::ImplType impltype =
      POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::GetImplType(oldele);

  if (impltype == INPAR::SCATRA::impltype_undefined)
    dserror(
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
POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::ConditionsToCopy() const
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

void POROELASTSCATRA::UTILS::PoroScatraCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = GLOBAL::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_elchmat) &&
      (mtype != INPAR::MAT::m_electrode) && (mtype != INPAR::MAT::m_matlist) &&
      (mtype != INPAR::MAT::m_matlist_reactions) && (mtype != INPAR::MAT::m_myocard) &&
      (mtype != INPAR::MAT::m_thermostvenant))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}

bool POROELASTSCATRA::UTILS::PoroelastCloneStrategyforScatraElements::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
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
BACI_NAMESPACE_CLOSE
