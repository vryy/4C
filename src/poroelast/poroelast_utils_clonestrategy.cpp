/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

 *----------------------------------------------------------------------*/

#include "poroelast_utils_clonestrategy.H"

#include "poroelast_utils.H"

#include "adapter_structure_scatra_ele.H"

#include "lib_globalproblem.H"

#include "fluid_ele_poro.H"

#include "mat_par_bundle.H"
#include "mat_fluidporo.H"
#include "mat_structporo.H"

#include "so3_nurbs27.H"
#include "so3_scatra.H"
#include "so3_poro_scatra.H"
#include "so3_poro_p1_scatra.H"

#include "w1_poro_scatra.H"
#include "w1_poro_p1_scatra.H"

#include "scatra_ele.H"

std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::ConditionsToCopy()
    const
{
  return {{"PoroDirichlet", "Dirichlet"}, {"PoroPointNeumann", "PointNeumann"},
      {"PoroLineNeumann", "LineNeumann"}, {"PoroSurfaceNeumann", "SurfaceNeumann"},
      {"PoroVolumeNeumann", "VolumeNeumann"}, {"NoPenetration", "NoPenetration"},
      {"PoroPartInt", "PoroPartInt"}, {"PoroCoupling", "PoroCoupling"},
      {"FSICoupling", "FSICoupling"}, {"FPSICoupling", "FPSICoupling"},
      {"PoroPresInt", "PoroPresInt"}, {"Mortar", "Mortar"}, {"SurfFlowRate", "SurfFlowRate"},
      {"LineFlowRate", "LineFlowRate"}, {"ImmersedSearchbox", "ImmersedSearchbox"},
      {"XFEMSurfFPIMono", "XFEMSurfFPIMono"}};
}

void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_fluidporo))
    dserror("Material with ID %d is not admissible for fluid poroelasticity elements", matid);
}

void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);
  if (fluid != Teuchos::null)
  {
    fluid->SetMaterial(matid);
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())
        ->SetInitialPorosity(
            Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->InitPorosity());
    fluid->SetDisType(oldele->Shape());  // set distype as well!
    fluid->SetIsAle(true);
    auto* so_base = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
    if (so_base)
      fluid->SetKinematicType(so_base->KinematicType());
    else
      dserror(" dynamic cast from DRT::Element* to DRT::ELEMENTS::So_base* failed ");

    SetAnisotropicPermeabilityDirectionsOntoFluid(newele, oldele);
    SetAnisotropicPermeabilityNodalCoeffsOntoFluid(newele, oldele);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
}

void POROELAST::UTILS::PoroelastCloneStrategy::SetAnisotropicPermeabilityDirectionsOntoFluid(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele)
{
  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = oldele->ElementType().Name();

  if (eletypename == "So_tet4PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_tet10PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_hex27PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad4>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallQuad9PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad9>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<DRT::Element::tri3>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }

  // Anisotropic permeability not yet supported for p1 or scatra type elements. Do nothing.
}

void POROELAST::UTILS::PoroelastCloneStrategy::SetAnisotropicPermeabilityNodalCoeffsOntoFluid(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele)
{
  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = oldele->ElementType().Name();

  if (eletypename == "So_tet4PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad4>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<DRT::Element::tri3>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }

  // Nodal anisotropic permeability not yet supported for higher order, p1, or scatra type elements.
  // Do nothing.
}

bool POROELAST::UTILS::PoroelastCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support fluid elements here
    eletype.emplace_back("FLUIDPORO");
    return true;
  }

  return false;
}

INPAR::SCATRA::ImplType POROELAST::UTILS::PoroScatraCloneStrategy::GetImplType(
    DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
)
{
  INPAR::SCATRA::ImplType impltype(INPAR::SCATRA::impltype_undefined);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = ele->ElementType().Name();

  // TET 4 Elements
  // tet 4 solid poro scatra
  if (eletypename == "So_tet4PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(
             ele))
            ->ImplType();
  }
  // tet4 solid porop1 scatra
  else if (eletypename == "So_tet4PoroP1ScatraType")
  {
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(ele))
            ->ImplType();
  }
  // tet 10 solid poro scatra
  else if (eletypename == "So_tet10PoroScatraType")
  {
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>*>(ele))
            ->ImplType();
  }
  // HEX 8 Elements
  // hex8 solid poro scatra
  else if (eletypename == "So_hex8PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
             ele))
            ->ImplType();
  }
  // hex8 solid porop1 scatra
  else if (eletypename == "So_hex8PoroP1ScatraType")
  {
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(ele))
            ->ImplType();
  }
  // hex27 solid poro scatra
  else if (eletypename == "So_hex27PoroScatraType")
  {
    impltype =
        (dynamic_cast<
             DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(ele))
            ->ImplType();
  }
  // nurbs 27
  else if (eletypename == "So_nurbs27PoroScatraType")
  {
    impltype = (dynamic_cast<DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27,
                    DRT::Element::nurbs27>*>(ele))
                   ->ImplType();
  }
  // wall poro scatra elements
  // quad 4
  else if (eletypename == "WallQuad4PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::quad4>*>(ele))->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::quad9>*>(ele))->ImplType();
  }
  // nurbs 4
  else if (eletypename == "WallNurbs4PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::nurbs4>*>(ele))->ImplType();
  }
  // nurbs 9
  else if (eletypename == "WallNurbs9PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::nurbs9>*>(ele))->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro_Scatra<DRT::Element::tri3>*>(ele))->ImplType();
  }
  // wall poro p1 elements
  // quad 4
  else if (eletypename == "WallQuad4PoroP1ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>*>(ele))->ImplType();
  }
  // quad 9
  else if (eletypename == "WallQuad9PoroP1ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>*>(ele))->ImplType();
  }
  // tri 3
  else if (eletypename == "WallTri3PoroP1ScatraType")
  {
    impltype =
        (dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>*>(ele))->ImplType();
  }
  // call base class routine
  else
    impltype = ADAPTER::GetScaTraImplType(ele);

  return impltype;
}

bool POROELAST::UTILS::PoroScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support transport elements here
    eletype.emplace_back("TRANSP");
    return true;
  }

  return false;
}

void POROELAST::UTILS::PoroScatraCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set material
    trans->SetMaterial(matid, oldele);
    // set distype as well!
    trans->SetDisType(oldele->Shape());

    // now check whether ImplType is reasonable and if set the ImplType
    INPAR::SCATRA::ImplType impltype =
        POROELAST::UTILS::PoroScatraCloneStrategy::GetImplType(oldele);
    if (impltype == INPAR::SCATRA::impltype_undefined)
    {
      dserror(
          "PoroScatraCloneStrategy copies scatra discretization from structure discretization, but "
          "the "
          "STRUCTURE elements that are defined in the .dat file are either not meant to be copied "
          "to scatra elements "
          "or the ImplType is set 'Undefined' which is not meaningful for the created scatra "
          "discretization! "
          "Use SOLIDSCATRA, WALLSCATRA, SHELLSCATRA, SOLIDPOROSCATRA, SOLIDPOROP1SCATRA, "
          "WALLPOROSCATRA or "
          "WALLPOROP1SCATRA Elements with meaningful ImplType instead!");
    }
    else
      trans->SetImplType(impltype);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
}

std::map<std::string, std::string> POROELAST::UTILS::PoroScatraCloneStrategy::ConditionsToCopy()
    const
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

void POROELAST::UTILS::PoroScatraCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_elchmat) &&
      (mtype != INPAR::MAT::m_electrode) && (mtype != INPAR::MAT::m_matlist) &&
      (mtype != INPAR::MAT::m_matlist_reactions) && (mtype != INPAR::MAT::m_myocard) &&
      (mtype != INPAR::MAT::m_thermostvenant))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}
