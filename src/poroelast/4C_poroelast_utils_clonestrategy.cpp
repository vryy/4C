/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

 *----------------------------------------------------------------------*/

#include "4C_poroelast_utils_clonestrategy.hpp"

#include "4C_fluid_ele_poro.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_so3_element_service.hpp"
#include "4C_so3_poro.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"
#include "4C_solid_poro_3D_ele.hpp"
#include "4C_w1_poro.hpp"

FOUR_C_NAMESPACE_OPEN

std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::conditions_to_copy()
    const
{
  return {{"PoroDirichlet", "Dirichlet"}, {"PoroPointNeumann", "PointNeumann"},
      {"PoroLineNeumann", "LineNeumann"}, {"PoroSurfaceNeumann", "SurfaceNeumann"},
      {"PoroVolumeNeumann", "VolumeNeumann"}, {"no_penetration", "no_penetration"},
      {"PoroPartInt", "PoroPartInt"}, {"PoroCoupling", "PoroCoupling"},
      {"FSICoupling", "FSICoupling"}, {"fpsi_coupling", "fpsi_coupling"},
      {"PoroPresInt", "PoroPresInt"}, {"Mortar", "Mortar"}, {"SurfFlowRate", "SurfFlowRate"},
      {"LineFlowRate", "LineFlowRate"}, {"ImmersedSearchbox", "ImmersedSearchbox"},
      {"XFEMSurfFPIMono", "XFEMSurfFPIMono"}, {"FluidNeumannInflow", "FluidNeumannInflow"}};
}

void POROELAST::UTILS::PoroelastCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  CORE::Materials::MaterialType mtype =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if ((mtype != CORE::Materials::m_fluidporo))
    FOUR_C_THROW("Material with ID %d is not admissible for fluid poroelasticity elements", matid);
}

void POROELAST::UTILS::PoroelastCloneStrategy::set_element_data(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele, const int matid,
    const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);
  if (fluid != Teuchos::null)
  {
    fluid->SetMaterial(0, MAT::Factory(matid));
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())
        ->SetInitialPorosity(
            Teuchos::rcp_static_cast<MAT::StructPoro>(oldele->Material())->InitPorosity());
    fluid->SetDisType(oldele->Shape());  // set distype as well!
    fluid->SetIsAle(true);
    auto* so_base = dynamic_cast<DRT::ELEMENTS::SoBase*>(oldele);
    if (so_base)
      fluid->SetKinematicType(so_base->KinematicType());
    else
      FOUR_C_THROW(
          " dynamic cast from CORE::Elements::Element* to DRT::ELEMENTS::So_base* failed ");

    set_anisotropic_permeability_directions_onto_fluid(newele, oldele);
    set_anisotropic_permeability_nodal_coeffs_onto_fluid(newele, oldele);
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
}

void POROELAST::UTILS::PoroelastCloneStrategy::set_anisotropic_permeability_directions_onto_fluid(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele)
{
  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = oldele->ElementType().Name();

  if (eletypename == "So_tet4PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(
             oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "So_tet10PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>*>(
             oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
             oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "So_hex27PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>*>(
             oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>*>(oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "WallQuad9PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad9>*>(oldele))
            ->get_anisotropic_permeability_directions());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->set_anisotropic_permeability_directions(
        (dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>*>(oldele))
            ->get_anisotropic_permeability_directions());
  }

  // Anisotropic permeability not yet supported for p1 type elements. Do nothing.
}

void POROELAST::UTILS::PoroelastCloneStrategy::set_anisotropic_permeability_nodal_coeffs_onto_fluid(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele)
{
  Teuchos::RCP<DRT::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidPoro>(newele);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = oldele->ElementType().Name();

  if (eletypename == "So_tet4PoroType")
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(
             oldele))
            ->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (eletypename == "So_hex8PoroType")
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        (dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
             oldele))
            ->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (eletypename == "WallQuad4PoroType")
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>*>(oldele))
            ->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (eletypename == "WallTri3PoroType")
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        (dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>*>(oldele))
            ->get_anisotropic_permeability_nodal_coeffs());
  }

  // Nodal anisotropic permeability not yet supported for higher order or p1 elements.
  // Do nothing.
}

bool POROELAST::UTILS::PoroelastCloneStrategy::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
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

FOUR_C_NAMESPACE_CLOSE
