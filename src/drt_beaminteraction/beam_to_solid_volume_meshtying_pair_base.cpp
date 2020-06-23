/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair_base.H"

#include "beam_contact_pair.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam,
    solid>::BeamToSolidVolumeMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::Setup()
{
  // Call setup of base class first.
  base_class::Setup();

  // Set reference nodal positions for the solid element
  for (unsigned int n = 0; n < solid::n_nodes_; ++n)
  {
    const DRT::Node* node = this->Element2()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_(3 * n + d) = node->X()[d];
  }

  // Initialize current nodal positions for the solid element
  for (unsigned int i = 0; i < solid::n_dof_; i++) ele2pos_(i) = 0.0;
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::CreateGeometryPair(
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Call the method of the base class.
  BeamContactPair::CreateGeometryPair(geometry_evaluation_data_ptr);

  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, beam, solid>(
      geometry_evaluation_data_ptr);
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    LINALG::Matrix<solid::n_dof_, 1, double> solid_coupling_ref;
    this->GetCouplingReferencePosition(beam_coupling_ref, solid_coupling_ref);
    CastGeometryPair()->PreEvaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Call the method in the parent class.
  base_class::ResetState(beam_centerline_dofvec, solid_nodal_dofvec);

  // Solid element.
  for (unsigned int i = 0; i < solid::n_dof_; i++)
  {
    ele2pos_(i) = FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + solid::n_dof_, beam::n_dof_ + i, solid_nodal_dofvec[i]);
  }
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::SetRestartDisplacement(
    const std::vector<std::vector<double>>& centerline_restart_vec_)
{
  // Call the parent method.
  base_class::SetRestartDisplacement(centerline_restart_vec_);

  // We only set the restart displacement, if the current section has the restart coupling flag.
  if (this->Params()->BeamToSolidVolumeMeshtyingParams()->GetCoupleRestartState())
  {
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
      ele1posref_offset_(i_dof) = centerline_restart_vec_[0][i_dof];

    // Add the displacement at the restart step to the solid reference position.
    for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
      ele2posref_offset_(i_dof) = centerline_restart_vec_[1][i_dof];
  }
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::GetPairVisualization(
    Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::GetPairVisualization(visualization_writer, visualization_params);

  // If a writer exists for segmentation point data, add the segmentation point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_segmentation =
          visualization_writer->GetVisualizationWriter("btsvc-segmentation");
  if (visualization_segmentation != Teuchos::null)
  {
    // Setup variables.
    LINALG::Matrix<3, 1, scalar_type> X;
    LINALG::Matrix<3, 1, scalar_type> u;
    LINALG::Matrix<3, 1, scalar_type> r;

    // Get the visualization vectors.
    std::vector<double>& point_coordinates =
        visualization_segmentation->GetMutablePointCoordinateVector();
    std::vector<double>& displacement =
        visualization_segmentation->GetMutablePointDataVector("displacement");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the left and right boundary point of the segment.
      for (const auto& segmentation_point : {segment.GetEtaA(), segment.GetEtaB()})
      {
        GEOMETRYPAIR::EvaluatePosition<beam>(
            segmentation_point, this->ele1posref_, X, this->Element1());
        GEOMETRYPAIR::EvaluatePosition<beam>(
            segmentation_point, this->ele1pos_, r, this->Element1());
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(FADUTILS::CastToDouble(X(dim)));
          displacement.push_back(FADUTILS::CastToDouble(u(dim)));
        }
      }
    }
  }

  // If a writer exists for integration point data, add the integration point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_integration_points =
          visualization_writer->GetVisualizationWriter("btsvc-integration-points");
  if (visualization_integration_points != Teuchos::null)
  {
    // Setup variables.
    LINALG::Matrix<3, 1, scalar_type> X;
    LINALG::Matrix<3, 1, scalar_type> u;
    LINALG::Matrix<3, 1, scalar_type> r;
    LINALG::Matrix<3, 1, scalar_type> r_solid;
    LINALG::Matrix<3, 1, scalar_type> force_integration_point;

    // Get the visualization vectors.
    std::vector<double>& point_coordinates =
        visualization_integration_points->GetMutablePointCoordinateVector();
    std::vector<double>& displacement =
        visualization_integration_points->GetMutablePointDataVector("displacement");
    std::vector<double>& force =
        visualization_integration_points->GetMutablePointDataVector("force");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the integration points.
      for (const auto& projection_point : segment.GetProjectionPoints())
      {
        this->EvaluateBeamPosition(projection_point, X, true);
        this->EvaluateBeamPosition(projection_point, r, false);
        u = r;
        u -= X;
        GEOMETRYPAIR::EvaluatePosition<solid>(
            projection_point.GetXi(), this->ele2pos_, r_solid, this->Element2());
        EvaluatePenaltyForce(r, r_solid, force_integration_point);
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(FADUTILS::CastToDouble(X(dim)));
          displacement.push_back(FADUTILS::CastToDouble(u(dim)));
          force.push_back(FADUTILS::CastToDouble(force_integration_point(dim)));
        }
      }
    }
  }
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::EvaluatePenaltyForce(
    const LINALG::Matrix<3, 1, scalar_type>& r_beam,
    const LINALG::Matrix<3, 1, scalar_type>& r_solid,
    LINALG::Matrix<3, 1, scalar_type>& force) const
{
  // TODO: call this function also in the Evaluate methods.
  // The base implementation of the force is a simple linear penalty law.
  force = r_solid;
  force -= r_beam;
  force.Scale(this->Params()->BeamToSolidVolumeMeshtyingParams()->GetPenaltyParameter());
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<beam, solid>::GetCouplingReferencePosition(
    LINALG::Matrix<beam::n_dof_, 1, double>& beam_coupling_ref,
    LINALG::Matrix<solid::n_dof_, 1, double>& solid_coupling_ref) const
{
  // Add the offset to the reference position.
  beam_coupling_ref = this->ele1posref_;
  beam_coupling_ref += this->ele1posref_offset_;
  solid_coupling_ref = ele2posref_;
  solid_coupling_ref += ele2posref_offset_;
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_tet10>;
  template class BeamToSolidVolumeMeshtyingPairBase<t_hermite, t_nurbs27>;
}  // namespace BEAMINTERACTION
