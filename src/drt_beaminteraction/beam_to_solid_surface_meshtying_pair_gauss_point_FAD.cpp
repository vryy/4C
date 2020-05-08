/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element. The
coupling terms are evaluated using FAD.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_gauss_point_FAD.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_calc_utils.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"

#include "Epetra_FEVector.h"


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<scalar_type, beam,
    surface>::BeamToSolidSurfaceMeshtyingPairGaussPointFAD()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<scalar_type, beam,
    surface>::EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix)
{
  dserror("Not yet implemented");
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
    Sacado::ELRFad::SLFad<Sacado::ELRFad::SLFad<double,
                              GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
        GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>;
