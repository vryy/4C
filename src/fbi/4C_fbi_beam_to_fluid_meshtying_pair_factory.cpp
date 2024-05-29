/*----------------------------------------------------------------------*/
/*! \file

\brief Create the pairs for beam to fluid meshtying depending on the input parameters..

\level 2
*/


#include "4C_fbi_beam_to_fluid_meshtying_pair_factory.hpp"

#include "4C_fbi_beam_to_fluid_meshtying_pair_gauss_point.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_pair_mortar.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_inpar_fbi.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> FBI::PairFactory::CreatePair(
    std::vector<CORE::Elements::Element const*> const& ele_ptrs,
    const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr)
{
  // Cast the fluid element.
  DRT::ELEMENTS::Fluid const* fluidele = dynamic_cast<DRT::ELEMENTS::Fluid const*>(ele_ptrs[1]);
  CORE::FE::CellType shape = fluidele->Shape();

  // Get the meshtying discretization method.
  INPAR::FBI::BeamToFluidDiscretization meshtying_discretization =
      params_ptr->get_contact_discretization();

  // Check which contact discretization is wanted.
  if (meshtying_discretization == INPAR::FBI::BeamToFluidDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case CORE::FE::CellType::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
      case CORE::FE::CellType::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
      case CORE::FE::CellType::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
      case CORE::FE::CellType::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
      case CORE::FE::CellType::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
      default:
        FOUR_C_THROW("Wrong element type for fluid element.");
    }
  }
  else if (meshtying_discretization == INPAR::FBI::BeamToFluidDiscretization::mortar)
  {
    INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions mortar_shape_function =
        params_ptr->get_mortar_shape_function_type();

    switch (mortar_shape_function)
    {
      case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line2:
      {
        switch (shape)
        {
          case CORE::FE::CellType::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>());
          case CORE::FE::CellType::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>());
          case CORE::FE::CellType::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>());
          case CORE::FE::CellType::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>());
          case CORE::FE::CellType::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>());
          default:
            FOUR_C_THROW("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line3:
      {
        switch (shape)
        {
          case CORE::FE::CellType::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>());
          case CORE::FE::CellType::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>());
          case CORE::FE::CellType::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>());
          case CORE::FE::CellType::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>());
          case CORE::FE::CellType::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>());
          default:
            FOUR_C_THROW("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::line4:
      {
        switch (shape)
        {
          case CORE::FE::CellType::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line4>());
          case CORE::FE::CellType::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line4>());
          case CORE::FE::CellType::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line4>());
          case CORE::FE::CellType::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line4>());
          case CORE::FE::CellType::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToFluidMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line4>());
          default:
            FOUR_C_THROW("Wrong element type for solid element.");
        }
        break;
      }
      default:
        FOUR_C_THROW("Wrong mortar shape function.");
    }
  }
  else
    FOUR_C_THROW("discretization type not yet implemented!\n");

  // Default return value.
  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
