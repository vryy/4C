/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

 *----------------------------------------------------------------------*/

#include "baci_poroelast_utils_clonestrategy.H"

#include "baci_fluid_ele_poro.H"
#include "baci_global_data.H"
#include "baci_mat_fluidporo.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_structporo.H"
#include "baci_poroelast_utils.H"
#include "baci_so3_element_service.H"
#include "baci_so3_poro.H"
#include "baci_so3_poro_p1_eletypes.H"
#include "baci_solid_poro_ele.H"
#include "baci_w1_poro.H"

BACI_NAMESPACE_OPEN

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
      {"XFEMSurfFPIMono", "XFEMSurfFPIMono"}, {"FluidNeumannInflow", "FluidNeumannInflow"}};
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
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_tet10PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "So_hex27PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>*>(
             oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<CORE::FE::CellType::quad4>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallQuad9PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<CORE::FE::CellType::quad9>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<CORE::FE::CellType::tri3>*>(oldele))
            ->GetAnisotropicPermeabilityDirections());
  }
  else if (eletypename == "SolidPoroType")
  {
    fluid->SetAnisotropicPermeabilityDirections(
        (dynamic_cast<DRT::ELEMENTS::SolidPoro*>(oldele))->GetAnisotropicPermeabilityDirections());
  }

  // Anisotropic permeability not yet supported for p1 type elements. Do nothing.
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
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>*>(
             oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>*>(
             oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<CORE::FE::CellType::quad4>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1_Poro<CORE::FE::CellType::tri3>*>(oldele))
            ->GetAnisotropicPermeabilityNodalCoeffs());
  }
  else if (eletypename == "SolidPoroType")
  {
    fluid->SetAnisotropicPermeabilityNodalCoeffs(
        (dynamic_cast<DRT::ELEMENTS::SolidPoro*>(oldele))->GetAnisotropicPermeabilityNodalCoeffs());
  }

  // Nodal anisotropic permeability not yet supported for higher order or p1 elements.
  // Do nothing.
}

bool POROELAST::UTILS::PoroelastCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
  if (IsPoroElement(actele))
  {
    // we only support fluid elements here
    eletype.emplace_back("FLUIDPORO");
    return true;
  }

  return false;
}

BACI_NAMESPACE_CLOSE
