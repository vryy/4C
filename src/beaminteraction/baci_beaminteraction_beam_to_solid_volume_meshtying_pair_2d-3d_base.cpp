/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for 2D-3D beam-to-solid volume mesh tying.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_base.hpp"

#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "baci_geometry_pair_line_to_volume_gauss_point_projection_cross_section.hpp"
#include "baci_geometry_pair_utility_classes.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DBase<beam, solid>::CreateGeometryPair(
    const DRT::Element* element1, const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_3d_evaluation_data = Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineTo3DEvaluationData>(
      geometry_evaluation_data_ptr, true);

  // Explicitly create the cross section projection geometry pair here and check that the correct
  // parameter is set in the input file.
  INPAR::GEOMETRYPAIR::LineTo3DStrategy strategy = line_to_3d_evaluation_data->GetStrategy();
  if (strategy != INPAR::GEOMETRYPAIR::LineTo3DStrategy::gauss_point_projection_cross_section)
    dserror(
        "The 2D-3D beam-to-volume mesh tying pair only works with the cross section projection "
        "geometry pair. This has to be specified in the input file.");
  this->geometry_pair_ = Teuchos::rcp(
      new GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, beam,
          solid>(element1, element2, line_to_3d_evaluation_data));
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DBase<beam,
    solid>::EvaluateBeamPositionDouble(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>&
                                           integration_point,
    CORE::LINALG::Matrix<3, 1, double>& r_beam, bool reference) const
{
  auto evaluate_position = [&](const auto& q, auto& r_beam)
  {
    const auto eta = integration_point.GetEta();
    CORE::LINALG::Matrix<3, 3, double> triad;
    GetTriadAtXiDouble(eta, triad, reference);
    CORE::LINALG::Matrix<3, 1, double> r_cross_section_ref, r_cross_section_cur;
    r_cross_section_ref(0) = 0.0;
    r_cross_section_ref(1) = integration_point.GetEtaCrossSection()(0);
    r_cross_section_ref(2) = integration_point.GetEtaCrossSection()(1);
    r_cross_section_cur.Multiply(triad, r_cross_section_ref);
    GEOMETRYPAIR::EvaluatePosition<beam>(eta, q, r_beam);
    r_beam += r_cross_section_cur;
  };

  if (reference)
  {
    evaluate_position(this->ele1posref_, r_beam);
  }
  else
  {
    evaluate_position(GEOMETRYPAIR::ElementDataToDouble<beam>::ToDouble(this->ele1pos_), r_beam);
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPair2D3DBase<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<t_hermite, t_tet10>;
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
