/*----------------------------------------------------------------------*/
/*! \file

\brief Create the pairs for beam to fluid meshtying depending on the input parameters..

\level 2
\maintainer Nora Hagmeyer
*/


#include "beam_to_fluid_meshtying_pair_factory.H"
#include "beam_to_fluid_meshtying_pair_gauss_point.H"
#include "beam_to_fluid_meshtying_params.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_inpar/inpar_fbi.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> FBI::PairFactory::CreatePair(
    std::vector<DRT::Element const*> const& ele_ptrs,
    const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr)
{
  // Cast the fluid element.
  DRT::ELEMENTS::Fluid const* fluidele = dynamic_cast<DRT::ELEMENTS::Fluid const*>(ele_ptrs[1]);
  DRT::Element::DiscretizationType shape = fluidele->Shape();

  // Get the meshtying discretization method.
  INPAR::FBI::BeamToFluidDiscretization meshtying_discretization =
      params_ptr->GetContactDiscretization();

  // Check which contact discretization is wanted.
  if (meshtying_discretization == INPAR::FBI::BeamToFluidDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
      default:
        dserror("Wrong element type for fluid element.");
    }
  }
  else
    dserror("Discretization type not yet implemented!\n");

  // Default return value.
  return Teuchos::null;
}
