/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for rotational meshtying between a 3D beam and a 3D solid element.

\level 3
*/


#include "beam_to_solid_volume_meshtying_pair_mortar_rotation.H"

#include "../drt_geometry_pair/geometry_pair_element_functions.H"


/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::BeamToSolidVolumeMeshtyingPairMortarRotation()
    : BeamToSolidVolumeMeshtyingPairMortar<beam, solid, mortar>()
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_GB,
    LINALG::SparseMatrix& global_GS, LINALG::SparseMatrix& global_FB,
    LINALG::SparseMatrix& global_FS, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

#define initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(            \
    mortar, mortar_rot)                                                                     \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex8, mortar,    \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex20, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex27, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_tet4, mortar,    \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_tet10, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_nurbs27, mortar, \
      mortar_rot>;

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line4);

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line4);

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line4);
}  // namespace BEAMINTERACTION
