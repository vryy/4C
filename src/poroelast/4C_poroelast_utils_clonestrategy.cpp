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
#include "4C_solid_poro_3D_ele_pressure_based.hpp"
#include "4C_solid_poro_3D_ele_pressure_velocity_based.hpp"
#include "4C_w1_poro.hpp"

FOUR_C_NAMESPACE_OPEN

std::map<std::string, std::string> PoroElast::Utils::PoroelastCloneStrategy::conditions_to_copy()
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

void PoroElast::Utils::PoroelastCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::instance()->materials()->parameter_by_id(matid)->type();
  if ((mtype != Core::Materials::m_fluidporo))
    FOUR_C_THROW("Material with ID %d is not admissible for fluid poroelasticity elements", matid);
}

void PoroElast::Utils::PoroelastCloneStrategy::set_element_data(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele, const int matid,
    const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  Teuchos::RCP<Discret::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::FluidPoro>(newele);
  if (fluid != Teuchos::null)
  {
    fluid->set_material(0, Mat::factory(matid));
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<Mat::PAR::FluidPoro*>(fluid->material()->parameter())
        ->set_initial_porosity(
            Teuchos::rcp_static_cast<Mat::StructPoro>(oldele->material())->init_porosity());
    fluid->set_dis_type(oldele->shape());  // set distype as well!
    fluid->set_is_ale(true);
    auto* solid_poro_pressure_velocity_based =
        dynamic_cast<Discret::ELEMENTS::SolidPoroPressureVelocityBased*>(oldele);
    auto* so_base = dynamic_cast<Discret::ELEMENTS::SoBase*>(oldele);
    if (solid_poro_pressure_velocity_based)
    {
      fluid->set_kinematic_type(solid_poro_pressure_velocity_based->kinematic_type());
    }
    else if (so_base)
      fluid->set_kinematic_type(so_base->kinematic_type());
    else
      FOUR_C_THROW(
          " dynamic cast from Core::Elements::Element* to Discret::ELEMENTS::So_base* or "
          "Discret::ELEMENTS::SolidPoroPressureVelocityBased failed ");

    set_anisotropic_permeability_directions_onto_fluid(newele, oldele);
    set_anisotropic_permeability_nodal_coeffs_onto_fluid(newele, oldele);
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
}

void PoroElast::Utils::PoroelastCloneStrategy::set_anisotropic_permeability_directions_onto_fluid(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele)
{
  Teuchos::RCP<Discret::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::FluidPoro>(newele);

  if (const auto* const so_tet4_poro_ele = dynamic_cast<
          Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>*>(oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        so_tet4_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const so_tet10_poro_ele = dynamic_cast<
               Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>*>(
               oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        so_tet10_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const so_hex8_poro_ele = dynamic_cast<
               Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>*>(
               oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        so_hex8_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const so_hex27_poro_ele = dynamic_cast<
               Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>*>(
               oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        so_hex27_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const wall1_quad4_poro_ele =
               dynamic_cast<Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::quad4>*>(oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        wall1_quad4_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const wall1_quad9_poro_ele =
               dynamic_cast<Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::quad9>*>(oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        wall1_quad9_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const wall1_tri3_poro_ele =
               dynamic_cast<Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::tri3>*>(oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        wall1_tri3_poro_ele->get_anisotropic_permeability_directions());
  }
  else if (const auto* const solid_poro_ele =
               dynamic_cast<const Discret::ELEMENTS::SolidPoroPressureVelocityBased* const>(oldele))
  {
    fluid->set_anisotropic_permeability_directions(
        solid_poro_ele->get_anisotropic_permeability_directions());
  }

  // Anisotropic permeability not yet supported for p1 type elements. Do nothing.
}

void PoroElast::Utils::PoroelastCloneStrategy::set_anisotropic_permeability_nodal_coeffs_onto_fluid(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele)
{
  Teuchos::RCP<Discret::ELEMENTS::FluidPoro> fluid =
      Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::FluidPoro>(newele);

  if (const auto* const so_tet4_poro_ele = dynamic_cast<
          Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>*>(oldele))
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        so_tet4_poro_ele->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (const auto* const so_hex8_poro_ele = dynamic_cast<
               Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>*>(
               oldele))
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        so_hex8_poro_ele->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (const auto* const wall1_hex8_poro_ele =
               dynamic_cast<const Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::quad4>* const>(
                   oldele))
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        wall1_hex8_poro_ele->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (const auto* const wall1_tri3_poro_ele =
               dynamic_cast<const Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::tri3>* const>(
                   oldele))
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        wall1_tri3_poro_ele->get_anisotropic_permeability_nodal_coeffs());
  }
  else if (const auto* const solid_poro_ele =
               dynamic_cast<const Discret::ELEMENTS::SolidPoroPressureVelocityBased* const>(oldele))
  {
    fluid->set_anisotropic_permeability_nodal_coeffs(
        solid_poro_ele->get_anisotropic_permeability_nodal_coeffs());
  }

  // Nodal anisotropic permeability not yet supported for higher order or p1 elements.
  // Do nothing.
}

bool PoroElast::Utils::PoroelastCloneStrategy::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element only if it is a poro element (we support submeshes here)
  if (is_poro_element(actele))
  {
    // we only support fluid elements here
    eletype.emplace_back("FLUIDPORO");
    return true;
  }

  return false;
}

FOUR_C_NAMESPACE_CLOSE
