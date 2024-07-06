/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a 3D fluid element.

\level 2

*/


#include "4C_fbi_beam_to_fluid_meshtying_pair_base.hpp"

#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_volume.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN



template <typename Beam, typename Fluid>
BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::BeamToFluidMeshtyingPairBase()
    : BeamToSolidVolumeMeshtyingPairBase<Beam, Fluid>()
{
  // Empty constructor.
}
/*------------------------------------------------------------------------------------------------*/

template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::setup()
{
  this->check_init();

  BeamToSolidVolumeMeshtyingPairBase<Beam, Fluid>::setup();

  // Initialize the element data containers
  ele1vel_ = GEOMETRYPAIR::InitializeElementData<Beam, scalar_type>::initialize(this->element1());
  ele2vel_ = GEOMETRYPAIR::InitializeElementData<Fluid, scalar_type>::initialize(this->element2());
  ele1poscur_ = GEOMETRYPAIR::InitializeElementData<Beam, double>::initialize(this->element1());
  ele2poscur_ = GEOMETRYPAIR::InitializeElementData<Fluid, double>::initialize(this->element2());

  // Initialize current nodal velocities for beam element
  for (unsigned int i = 0; i < Beam::n_dof_; i++) this->ele1vel_.element_position_(i) = 0.0;

  // Initialize current nodal velocities for fluid element
  for (unsigned int i = 0; i < Fluid::n_dof_; i++) this->ele2vel_.element_position_(i) = 0.0;

  this->issetup_ = true;
}

/**
 *
 */
template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam,
    Fluid>::reset_state(  // todo somehow hand in nodal velocities
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& fluid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < Beam::n_dof_; i++)
  {
    this->ele1pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        Beam::n_dof_ + Fluid::n_dof_, i, beam_centerline_dofvec[i]);
    this->ele1poscur_.element_position_(i) = beam_centerline_dofvec[i];
    this->ele1vel_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        Beam::n_dof_ + Fluid::n_dof_, i, beam_centerline_dofvec[Beam::n_dof_ + i]);
  }

  // Fluid element.
  for (unsigned int i = 0; i < Fluid::n_dof_; i++)
  {
    this->ele2pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        Beam::n_dof_ + Fluid::n_dof_, Beam::n_dof_ + i, fluid_nodal_dofvec[i]);
    this->ele2poscur_.element_position_(i) = fluid_nodal_dofvec[i];
    this->ele2vel_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        Beam::n_dof_ + Fluid::n_dof_, Beam::n_dof_ + i, fluid_nodal_dofvec[Fluid::n_dof_ + i]);
  }
}

template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::create_geometry_pair(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Set up the geometry pair
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, Beam, Fluid>(
      element1, element2, geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::pre_evaluate()
{
  // Call pre_evaluate on the geometry Pair.
  if (!this->meshtying_is_evaluated_)
  {
    this->cast_geometry_pair()->pre_evaluate(ele1poscur_, ele2poscur_, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::print(std::ostream& out) const
{
  this->check_init_setup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToFluidMeshtyingPair"
      << "\nBeam EleGID:  " << this->element1()->id()
      << "\nFluid EleGID: " << this->element2()->id();

  out << "\n\nele1 dofvec: " << this->ele1pos_.element_position_;
  out << "\nele2 dofvec: " << this->ele2pos_.element_position_;
  out << "\nn_segments: " << this->line_to_3D_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}


/**
 *
 */
template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam,
    Fluid>::print_summary_one_line_per_active_segment_pair(std::ostream& out) const
{
  this->check_init_setup();

  // Only display information if a segment exists for this pair.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << this->element1()->id() << ", fluid ID " << this->element2()->id() << ":";
  out << " n_segments = " << this->line_to_3D_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < this->line_to_3D_segments_.size();
       index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << this->line_to_3D_segments_[index_segment].get_etadata() << ", "
        << this->line_to_3D_segments_[index_segment].get_eta_b() << "]";
    out << ", Gauss points = "
        << this->line_to_3D_segments_[index_segment].get_number_of_projection_points();
    out << "\n";
  }
}

/**
 *
 */
template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::get_pair_visualization(
    Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  BeamContactPair::get_pair_visualization(visualization_writer, visualization_params);


  // If a writer exists for integration point data, add the integration point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization =
      visualization_writer->get_visualization_writer("integration-points");
  if (visualization != Teuchos::null)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, scalar_type> X;
    Core::LinAlg::Matrix<3, 1, scalar_type> u;
    Core::LinAlg::Matrix<3, 1, scalar_type> r;
    Core::LinAlg::Matrix<3, 1, scalar_type> r_fluid;
    Core::LinAlg::Matrix<3, 1, scalar_type> v_beam;
    Core::LinAlg::Matrix<3, 1, scalar_type> force_integration_point;

    // Get the visualization vectors.
    auto& visualization_data = visualization->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the integration points.
      for (const auto& projection_point : segment.get_projection_points())
      {
        evaluate_beam_position(projection_point, X, true);
        evaluate_beam_position(projection_point, r, false);
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(Core::FADUtils::CastToDouble(X(dim)));
          displacement.push_back(Core::FADUtils::CastToDouble(u(dim)));
        }
      }
    }
  }


  // If a writer exists for segmentation point data, add the segmentation point data.
  visualization = visualization_writer->get_visualization_writer("segmentation");
  if (visualization != Teuchos::null)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, scalar_type> X;
    Core::LinAlg::Matrix<3, 1, scalar_type> u;
    Core::LinAlg::Matrix<3, 1, scalar_type> r;

    // Get the visualization vectors.
    auto& visualization_data = visualization->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the left and right boundary point of the segment.
      for (const auto& segmentation_point : {segment.get_etadata(), segment.get_eta_b()})
      {
        GEOMETRYPAIR::EvaluatePosition<Beam>(segmentation_point, this->ele1posref_, X);
        GEOMETRYPAIR::EvaluatePosition<Beam>(segmentation_point, this->ele1pos_, r);
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(Core::FADUtils::CastToDouble(X(dim)));
          displacement.push_back(Core::FADUtils::CastToDouble(u(dim)));
        }
      }
    }
  }
}

template <typename Beam, typename Fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<Beam, Fluid>::evaluate_beam_position(
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
    Core::LinAlg::Matrix<3, 1, scalar_type>& r_beam, bool reference) const
{
  if (reference)
    GEOMETRYPAIR::EvaluatePosition<Beam>(integration_point.get_eta(), this->ele1posref_, r_beam);
  else
    GEOMETRYPAIR::EvaluatePosition<Beam>(integration_point.get_eta(), this->ele1pos_, r_beam);
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

FOUR_C_NAMESPACE_CLOSE
