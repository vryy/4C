/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a 3D fluid element.

\level 2

*/


#include "baci_fbi_beam_to_fluid_meshtying_pair_base.hpp"

#include "baci_beam3_euler_bernoulli.hpp"
#include "baci_beam3_kirchhoff.hpp"
#include "baci_beam3_reissner.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "baci_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "baci_beaminteraction_contact_pair.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_geometry_pair_element_functions.hpp"
#include "baci_geometry_pair_factory.hpp"
#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "baci_geometry_pair_line_to_volume.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

BACI_NAMESPACE_OPEN



template <typename beam, typename fluid>
BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::BeamToFluidMeshtyingPairBase()
    : BeamToSolidVolumeMeshtyingPairBase<beam, fluid>()
{
  // Empty constructor.
}
/*------------------------------------------------------------------------------------------------*/

template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::Setup()
{
  this->CheckInit();

  BeamToSolidVolumeMeshtyingPairBase<beam, fluid>::Setup();

  // Initialize current nodal velocities for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) this->ele1vel_(i) = 0.0;

  // Initialize current nodal velocities for fluid element
  for (unsigned int i = 0; i < fluid::n_dof_; i++) this->ele2vel_(i) = 0.0;

  this->issetup_ = true;
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam,
    fluid>::ResetState(  // todo somehow hand in nodal velocities
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& fluid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < beam::n_dof_; i++)
  {
    this->ele1pos_(i) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + fluid::n_dof_, i, beam_centerline_dofvec[i]);
    this->ele1poscur_(i) = beam_centerline_dofvec[i];
    this->ele1vel_(i) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + fluid::n_dof_, i, beam_centerline_dofvec[beam::n_dof_ + i]);
  }

  // Fluid element.
  for (unsigned int i = 0; i < fluid::n_dof_; i++)
  {
    this->ele2pos_(i) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + fluid::n_dof_, beam::n_dof_ + i, fluid_nodal_dofvec[i]);
    this->ele2poscur_(i) = fluid_nodal_dofvec[i];
    this->ele2vel_(i) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + fluid::n_dof_, beam::n_dof_ + i, fluid_nodal_dofvec[fluid::n_dof_ + i]);
  }
}

template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::CreateGeometryPair(
    const DRT::Element* element1, const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Set up the geometry pair
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, beam, fluid>(
      element1, element2, geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->PreEvaluate(ele1poscur_, ele2poscur_, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::Print(std::ostream& out) const
{
  this->CheckInitSetup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToFluidMeshtyingPair"
      << "\nBeam EleGID:  " << this->Element1()->Id()
      << "\nFluid EleGID: " << this->Element2()->Id();

  out << "\n\nele1 dofvec: " << this->ele1pos_;
  out << "\nele2 dofvec: " << this->ele2pos_;
  out << "\nn_segments: " << this->line_to_3D_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam,
    fluid>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  this->CheckInitSetup();

  // Only display information if a segment exists for this pair.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << this->Element1()->Id() << ", fluid ID " << this->Element2()->Id() << ":";
  out << " n_segments = " << this->line_to_3D_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < this->line_to_3D_segments_.size();
       index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << this->line_to_3D_segments_[index_segment].GetEtaA() << ", "
        << this->line_to_3D_segments_[index_segment].GetEtaB() << "]";
    out << ", Gauss points = "
        << this->line_to_3D_segments_[index_segment].GetNumberOfProjectionPoints();
    out << "\n";
  }
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::GetPairVisualization(
    Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  BeamContactPair::GetPairVisualization(visualization_writer, visualization_params);


  // If a writer exists for integration point data, add the integration point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization =
      visualization_writer->GetVisualizationWriter("integration-points");
  if (visualization != Teuchos::null)
  {
    // Setup variables.
    CORE::LINALG::Matrix<3, 1, scalar_type> X;
    CORE::LINALG::Matrix<3, 1, scalar_type> u;
    CORE::LINALG::Matrix<3, 1, scalar_type> r;
    CORE::LINALG::Matrix<3, 1, scalar_type> r_fluid;
    CORE::LINALG::Matrix<3, 1, scalar_type> v_beam;
    CORE::LINALG::Matrix<3, 1, scalar_type> force_integration_point;

    // Get the visualization vectors.
    auto& visualization_data = visualization->GetVisualizationData();
    std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates();
    std::vector<double>& displacement = visualization_data.GetPointData<double>("displacement");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the integration points.
      for (const auto& projection_point : segment.GetProjectionPoints())
      {
        EvaluateBeamPosition(projection_point, X, true);
        EvaluateBeamPosition(projection_point, r, false);
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(CORE::FADUTILS::CastToDouble(X(dim)));
          displacement.push_back(CORE::FADUTILS::CastToDouble(u(dim)));
        }
      }
    }
  }


  // If a writer exists for segmentation point data, add the segmentation point data.
  visualization = visualization_writer->GetVisualizationWriter("segmentation");
  if (visualization != Teuchos::null)
  {
    // Setup variables.
    CORE::LINALG::Matrix<3, 1, scalar_type> X;
    CORE::LINALG::Matrix<3, 1, scalar_type> u;
    CORE::LINALG::Matrix<3, 1, scalar_type> r;

    // Get the visualization vectors.
    auto& visualization_data = visualization->GetVisualizationData();
    std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates();
    std::vector<double>& displacement = visualization_data.GetPointData<double>("displacement");

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
          point_coordinates.push_back(CORE::FADUTILS::CastToDouble(X(dim)));
          displacement.push_back(CORE::FADUTILS::CastToDouble(u(dim)));
        }
      }
    }
  }
}

template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::EvaluateBeamPosition(
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
    CORE::LINALG::Matrix<3, 1, scalar_type>& r_beam, bool reference) const
{
  if (reference)
    GEOMETRYPAIR::EvaluatePosition<beam>(
        integration_point.GetEta(), this->ele1posref_, r_beam, this->Element1());
  else
    GEOMETRYPAIR::EvaluatePosition<beam>(
        integration_point.GetEta(), this->ele1pos_, r_beam, this->Element1());
}

/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;

BACI_NAMESPACE_CLOSE
