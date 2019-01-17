/*!
\file beam_to_solid_volume_meshtying_pair_factory.cpp

\brief Create the pairs for beam to solid volume meshtying depending on the input parameters..

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair_factory.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point.H"
#include "beam_to_solid_volume_meshtying_params.H"

#include "../drt_so3/so_base.H"
#include "../drt_inpar/inpar_beaminteraction.H"


Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairFactory(
    std::vector<DRT::Element const*> const& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr)
{
  // Cast the solid element.
  DRT::ELEMENTS::So_base const* solidele = dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]);
  DRT::Element::DiscretizationType shape = solidele->Shape();

  // Get the contact discretization method.
  INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization contact_discretization =
      params_ptr->BeamToSolidVolumeMeshtyingParams()->GetContactDiscretization();

  // Check which contact discretization is wanted.
  if (contact_discretization ==
      INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
      default:
        dserror("Wrong element type for solid element.");
    }
  }
  else if (contact_discretization ==
           INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization::gauss_point_to_segment)
  {
    dserror("Mortar not yet implemented.");
  }
  else
  {
    dserror("Wrong contact discretization for beam to solid volume meshtying.");
  }

  // Default return value.
  return Teuchos::null;
}
