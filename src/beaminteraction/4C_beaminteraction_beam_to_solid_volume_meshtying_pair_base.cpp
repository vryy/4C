/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_base.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_volume.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename Beam, typename Solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam,
    Solid>::BeamToSolidVolumeMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::setup()
{
  // Call setup of base class first.
  base_class::setup();

  // Get the solid element data container
  ele2posref_ = GEOMETRYPAIR::InitializeElementData<Solid, double>::initialize(this->element2());
  ele2pos_ = GEOMETRYPAIR::InitializeElementData<Solid, scalar_type>::initialize(this->element2());

  // Set reference nodal positions for the solid element
  for (unsigned int n = 0; n < Solid::n_nodes_; ++n)
  {
    const Core::Nodes::Node* node = this->element2()->nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_.element_position_(3 * n + d) = node->x()[d];
  }

  // Initialize current nodal positions for the solid element
  for (unsigned int i = 0; i < Solid::n_dof_; i++) ele2pos_.element_position_(i) = 0.0;
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::create_geometry_pair(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  this->geometry_pair_ = GEOMETRYPAIR::geometry_pair_line_to_volume_factory<double, Beam, Solid>(
      element1, element2, geometry_evaluation_data_ptr);
}


/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::pre_evaluate()
{
  // Call pre_evaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    GEOMETRYPAIR::ElementData<Beam, double> beam_coupling_ref;
    GEOMETRYPAIR::ElementData<Solid, double> solid_coupling_ref;
    this->get_coupling_reference_position(beam_coupling_ref, solid_coupling_ref);
    cast_geometry_pair()->pre_evaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::reset_state(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Call the method in the parent class.
  base_class::reset_state(beam_centerline_dofvec, solid_nodal_dofvec);

  // Solid element.
  for (unsigned int i = 0; i < Solid::n_dof_; i++)
  {
    ele2pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<scalar_type>::apply(
        Beam::n_dof_ + Solid::n_dof_, Beam::n_dof_ + i, solid_nodal_dofvec[i]);
  }
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::set_restart_displacement(
    const std::vector<std::vector<double>>& centerline_restart_vec_)
{
  // Call the parent method.
  base_class::set_restart_displacement(centerline_restart_vec_);

  // We only set the restart displacement, if the current section has the restart coupling flag.
  if (this->params()->beam_to_solid_volume_meshtying_params()->get_couple_restart_state())
  {
    for (unsigned int i_dof = 0; i_dof < Beam::n_dof_; i_dof++)
      ele1posref_offset_(i_dof) = centerline_restart_vec_[0][i_dof];

    // Add the displacement at the restart step to the solid reference position.
    for (unsigned int i_dof = 0; i_dof < Solid::n_dof_; i_dof++)
      ele2posref_offset_(i_dof) = centerline_restart_vec_[1][i_dof];
  }
}


/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam, Solid>::get_pair_visualization(
    Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::get_pair_visualization(visualization_writer, visualization_params);

  // Get the writers.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_segmentation =
      visualization_writer->get_visualization_writer("btsvc-segmentation");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
      visualization_integration_points =
          visualization_writer->get_visualization_writer("btsvc-integration-points");
  if (visualization_segmentation.is_null() and visualization_integration_points.is_null()) return;

  const Teuchos::RCP<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>& output_params_ptr =
      visualization_params
          .get<Teuchos::RCP<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>>(
              "btsvc-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->get_write_unique_i_ds_flag();

  if (visualization_segmentation != Teuchos::null)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, scalar_type> X;
    Core::LinAlg::Matrix<3, 1, scalar_type> u;
    Core::LinAlg::Matrix<3, 1, scalar_type> r;

    // Get the visualization vectors.
    auto& visualization_data = visualization_segmentation->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");

    std::vector<double>* pair_beam_id = nullptr;
    std::vector<double>* pair_solid_id = nullptr;
    if (write_unique_ids)
    {
      pair_beam_id = &(visualization_data.get_point_data<double>("uid_0_pair_beam_id"));
      pair_solid_id = &(visualization_data.get_point_data<double>("uid_1_pair_solid_id"));
    }

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the left and right boundary point of the segment.
      for (const auto& segmentation_point : {segment.get_etadata(), segment.get_eta_b()})
      {
        GEOMETRYPAIR::evaluate_position<Beam>(segmentation_point, this->ele1posref_, X);
        GEOMETRYPAIR::evaluate_position<Beam>(segmentation_point, this->ele1pos_, r);
        u = r;
        u -= X;
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(Core::FADUtils::cast_to_double(X(dim)));
          displacement.push_back(Core::FADUtils::cast_to_double(u(dim)));
        }

        if (write_unique_ids)
        {
          pair_beam_id->push_back(this->element1()->id());
          pair_solid_id->push_back(this->element2()->id());
        }
      }
    }
  }

  // If a writer exists for integration point data, add the integration point data.
  if (visualization_integration_points != Teuchos::null)
  {
    // Setup variables.
    Core::LinAlg::Matrix<3, 1, double> X;
    Core::LinAlg::Matrix<3, 1, double> u;
    Core::LinAlg::Matrix<3, 1, double> r;
    Core::LinAlg::Matrix<3, 1, double> r_solid;
    Core::LinAlg::Matrix<3, 1, double> force_integration_point;

    // Get the visualization vectors.
    auto& visualization_data = visualization_integration_points->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");
    std::vector<double>& force = visualization_data.get_point_data<double>("force");

    std::vector<double>* pair_beam_id = nullptr;
    std::vector<double>* pair_solid_id = nullptr;
    if (write_unique_ids)
    {
      pair_beam_id = &(visualization_data.get_point_data<double>("uid_0_pair_beam_id"));
      pair_solid_id = &(visualization_data.get_point_data<double>("uid_1_pair_solid_id"));
    }

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the integration points.
      for (const auto& projection_point : segment.get_projection_points())
      {
        this->evaluate_beam_position_double(projection_point, X, true);
        this->evaluate_beam_position_double(projection_point, r, false);
        u = r;
        u -= X;
        GEOMETRYPAIR::evaluate_position<Solid>(projection_point.get_xi(),
            GEOMETRYPAIR::ElementDataToDouble<Solid>::to_double(this->ele2pos_), r_solid);
        evaluate_penalty_force_double(r, r_solid, force_integration_point);
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(X(dim));
          displacement.push_back(u(dim));
          force.push_back(force_integration_point(dim));
        }

        if (write_unique_ids)
        {
          pair_beam_id->push_back(this->element1()->id());
          pair_solid_id->push_back(this->element2()->id());
        }
      }
    }
  }
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam,
    Solid>::evaluate_penalty_force_double(const Core::LinAlg::Matrix<3, 1, double>& r_beam,
    const Core::LinAlg::Matrix<3, 1, double>& r_solid,
    Core::LinAlg::Matrix<3, 1, double>& force) const
{
  // The base implementation of the force is a simple linear penalty law.
  force = r_solid;
  force -= r_beam;
  force.scale(this->params()->beam_to_solid_volume_meshtying_params()->get_penalty_parameter());
}

/**
 *
 */
template <typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairBase<Beam,
    Solid>::get_coupling_reference_position(GEOMETRYPAIR::ElementData<Beam, double>&
                                                beam_coupling_ref,
    GEOMETRYPAIR::ElementData<Solid, double>& solid_coupling_ref) const
{
  // Add the offset to the reference position.
  beam_coupling_ref = this->ele1posref_;
  beam_coupling_ref.element_position_ += this->ele1posref_offset_;
  solid_coupling_ref = ele2posref_;
  solid_coupling_ref.element_position_ += ele2posref_offset_;
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

FOUR_C_NAMESPACE_CLOSE
