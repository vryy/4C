/*!
\file beam_to_solid_volume_meshtying_pair.cpp

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair.H"

#include "beam_contact_pair.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "../linalg/linalg_utils.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::BeamToSolidVolumeMeshtyingPair()
    : BeamContactPair(), meshtying_is_evaluated_(false), line_to_volume_segments_()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
    std::vector<DRT::Element const*> elements)
{
  // Call Init of base class, the geometry pair will be created and initialized there.
  BeamContactPair::Init(params_ptr, geometry_evaluation_data_ptr, elements);
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::Setup()
{
  CheckInit();

  // Call setup of base class first.
  BeamContactPair::Setup();

  // Initialize reference nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1posref_(i) = 0.0;

  // Set reference nodal positions (and tangents) for beam element
  for (unsigned int n = 0; n < beam::n_nodes_; ++n)
  {
    const DRT::Node* node = Element1()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele1posref_(3 * beam::n_val_ * n + d) = node->X()[d];

    // tangents
    if (beam::n_val_ == 2)
    {
      LINALG::Matrix<3, 1> tan;
      const DRT::ElementType& eot = Element1()->ElementType();
      if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        dserror("ERROR: Beam3tosolidmeshtying: beam::n_val_=2 detected for beam3 element");
      }
      else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(Element1());
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror(
              "ERROR: Beam3tosolidmeshtying: beam::n_val_=2 detected for beam3r element w/o "
              "Hermite CL");
      }
      else if (eot == DRT::ELEMENTS::Beam3kType::Instance())
      {
        const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(Element1());
        tan = ele->Tref()[n];
      }
      else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
      {
        const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(Element1());
        tan = ele->Tref()[n];
      }
      else
      {
        dserror("ERROR: Beam3tosolidmeshtying: Invalid beam element type");
      }

      for (int d = 0; d < 3; ++d) ele1posref_(3 * beam::n_val_ * n + d + 3) = tan(d, 0);
    }
  }

  // Initialize reference nodal positions for solid surface element
  for (unsigned int i = 0; i < solid::n_dof_; i++) ele2posref_(i) = 0.0;

  // Set reference nodal positions for solid surface element
  for (unsigned int n = 0; n < solid::n_nodes_; ++n)
  {
    const DRT::Node* node = Element2()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_(3 * n + d) = node->X()[d];
  }

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1pos_(i) = 0.0;

  // Initialize current nodal positions for solid surface element
  for (unsigned int i = 0; i < solid::n_dof_; i++) ele2pos_(i) = 0.0;

  issetup_ = true;
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::CreateGeometryPair(
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
{
  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, beam, solid>(
      geometry_evaluation_data_ptr);
}


/**
 *
 */
template <typename beam, typename solid>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::CastGeometryPair() const
{
  return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>(
      geometry_pair_, true);
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    CastGeometryPair()->PreEvaluate(ele1posref_, ele2posref_, line_to_volume_segments_);
  }
}


/**
 *
 */
template <typename beam, typename solid>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    CastGeometryPair()->Evaluate(ele1posref_, ele2posref_, line_to_volume_segments_);
    meshtying_is_evaluated_ = true;
  }

  // Initialize variables for position and force vectors.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_beam;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_solid;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> force;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, beam::n_dof_, 1> force_element_1(true);
  LINALG::TMatrix<TYPE_BTS_VMT_AD, solid::n_dof_, 1> force_element_2(true);

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;
  double penalty_parameter = Params()->BeamToSolidVolumeMeshtyingParams()->GetPenaltyParameter();

  // Calculate the meshtying forces.
  // Loop over segments.
  for (unsigned int i_segment = 0; i_segment < line_to_volume_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * line_to_volume_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < line_to_volume_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& projected_gauss_point =
          line_to_volume_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the jacobian in the reference configuration.
      CastGeometryPair()->GetElement1PositionDerivative(
          projected_gauss_point.GetEta(), ele1posref_, dr_beam_ref);

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the current positions on beam and solid.
      CastGeometryPair()->GetElement1Position(projected_gauss_point.GetEta(), ele1pos_, r_beam);
      CastGeometryPair()->GetElement2Position(projected_gauss_point.GetXi(), ele2pos_, r_solid);

      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the beam.
      force = r_solid;
      force -= r_beam;
      force.Scale(penalty_parameter);

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
      // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_2(i_dof) -= force(i_dir) * r_solid(i_dir).dx(i_dof + beam::n_dof_) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }


  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != NULL) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != NULL) forcevec2->Size(solid::n_dof_);
    if (stiffmat11 != NULL) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != NULL) stiffmat12->Shape(beam::n_dof_, solid::n_dof_);
    if (stiffmat21 != NULL) stiffmat21->Shape(solid::n_dof_, beam::n_dof_);
    if (stiffmat22 != NULL) stiffmat22->Shape(solid::n_dof_, solid::n_dof_);

    if (forcevec1 != NULL && forcevec2 != NULL)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = force_element_1(i_dof).val();
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        (*forcevec2)(i_dof) = force_element_2(i_dof).val();
    }

    if (stiffmat11 != NULL && stiffmat12 != NULL && stiffmat21 != NULL && stiffmat22 != NULL)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < beam::n_dof_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(i_dof_2);

      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2);
          (*stiffmat21)(i_dof_2, i_dof_1) = -force_element_2(i_dof_2).dx(i_dof_1);
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < solid::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) = -force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2);
    }
  }

  // Check if there are meshtying contributions.
  if ((forcevec1 != NULL and forcevec1->NormInf() > 0) or
      (forcevec2 != NULL and forcevec2->NormInf() > 0) or
      (stiffmat11 != NULL and stiffmat11->NormInf() > 0) or
      (stiffmat12 != NULL and stiffmat12->NormInf() > 0) or
      (stiffmat21 != NULL and stiffmat21->NormInf() > 0) or
      (stiffmat22 != NULL and stiffmat22->NormInf() > 0))
    return true;
  else
    return false;
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < beam::n_dof_; i++)
  {
    ele1pos_(i) = TYPE_BTS_VMT_AD(beam::n_dof_ + solid::n_dof_, i, beam_centerline_dofvec[i]);
  }

  // Solid element.
  for (unsigned int i = 0; i < solid::n_dof_; i++)
  {
    ele2pos_(i) =
        TYPE_BTS_VMT_AD(beam::n_dof_ + solid::n_dof_, beam::n_dof_ + i, solid_nodal_dofvec[i]);
  }
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam, solid>::Print(std::ostream& out) const
{
  CheckInitSetup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToSolidVolumeMeshtyingPair"
      << "\nBeam EleGID:  " << Element1()->Id() << "\nSolid EleGID: " << Element2()->Id();

  out << "\n\nele1 dofvec: " << ele1pos_;
  out << "\nele2 dofvec: " << ele2pos_;
  out << "\nn_segments: " << line_to_volume_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<beam,
    solid>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Only display information if a segment exists for this pair.
  if (line_to_volume_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << Element1()->Id() << ", solid ID " << Element2()->Id() << ":";
  out << " n_segments = " << line_to_volume_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < line_to_volume_segments_.size();
       index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << line_to_volume_segments_[index_segment].GetEtaA() << ", "
        << line_to_volume_segments_[index_segment].GetEtaB() << "]";
    out << ", Gauss points = "
        << line_to_volume_segments_[index_segment].GetNumberOfProjectionPoints();
    out << "\n";
  }
}


/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
