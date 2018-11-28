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
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::BeamToSolidVolumeMeshtyingPair()
    : BeamContactPair(), line_to_volume_segments_()
{
  // Empty constructor.
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Init(
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
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Setup()
{
  CheckInit();

  // Call setup of base class first.
  BeamContactPair::Setup();

  // Initialize reference nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; i++) ele1posref_(i) = 0.0;

  // Set reference nodal positions (and tangents) for beam element
  for (unsigned int n = 0; n < numnodes; ++n)
  {
    const DRT::Node* node = Element1()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele1posref_(3 * numnodalvalues * n + d) = node->X()[d];

    // tangents
    if (numnodalvalues == 2)
    {
      LINALG::Matrix<3, 1> tan;
      const DRT::ElementType& eot = Element1()->ElementType();
      if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        dserror("ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3 element");
      }
      else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(Element1());
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror(
              "ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3r element w/o "
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

      for (int d = 0; d < 3; ++d) ele1posref_(3 * numnodalvalues * n + d + 3) = tan(d, 0);
    }
  }

  // Initialize reference nodal positions for solid surface element
  for (unsigned int i = 0; i < 3 * numnodessol; i++) ele2posref_(i) = 0.0;

  // Set reference nodal positions for solid surface element
  for (unsigned int n = 0; n < numnodessol; ++n)
  {
    const DRT::Node* node = Element2()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_(3 * n + d) = node->X()[d];
  }

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; i++) ele1pos_(i) = 0.0;

  // Initialize current nodal positions for solid surface element
  for (unsigned int i = 0; i < 3 * numnodessol; i++) ele2pos_(i) = 0.0;

  issetup_ = true;
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::
    CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
{
  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, numnodes, numnodalvalues,
      numnodessol, 1>(geometry_evaluation_data_ptr);
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
Teuchos::RCP<
    GEOMETRYPAIR::GeometryPairLineToVolume<double, numnodes, numnodalvalues, numnodessol, 1>>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::CastGeometryPair() const
{
  return Teuchos::rcp_dynamic_cast<
      GEOMETRYPAIR::GeometryPairLineToVolume<double, numnodes, numnodalvalues, numnodessol, 1>>(
      geometry_pair_, true);
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  CastGeometryPair()->PreEvaluate(ele1posref_, ele2posref_, line_to_volume_segments_);
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::Evaluate(LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair.
  CastGeometryPair()->Evaluate(ele1posref_, ele2posref_, line_to_volume_segments_);

  // Initialize variables for position and force vectors.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_beam;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_solid;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> force;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, n_dof_element_1_, 1> force_element_1(true);
  LINALG::TMatrix<TYPE_BTS_VMT_AD, n_dof_element_2_, 1> force_element_2(true);

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;
  double penalty_parameter = Params()->BeamToSolidVolumeMeshtyingParams()->PenaltyParameter();

  // Calculate the meshtying forces.
  // Loop over segments.
  for (unsigned int i_segment = 0; i_segment < line_to_volume_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * line_to_volume_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0; i_gp < line_to_volume_segments_[i_segment].GetGaussPoints().size();
         i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& projected_gauss_point =
          line_to_volume_segments_[i_segment].GetGaussPoints()[i_gp];

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
      for (unsigned int i_dof = 0; i_dof < n_dof_element_1_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < n_dof_element_2_; i_dof++)
      {
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        {
          force_element_2(i_dof) -= force(i_dir) * r_solid(i_dir).dx(i_dof + n_dof_element_1_) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
        }
      }
    }
  }


  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != NULL) forcevec1->Size(n_dof_element_1_);
    if (forcevec2 != NULL) forcevec2->Size(n_dof_element_2_);
    if (stiffmat11 != NULL) stiffmat11->Shape(n_dof_element_1_, n_dof_element_1_);
    if (stiffmat12 != NULL) stiffmat12->Shape(n_dof_element_1_, n_dof_element_2_);
    if (stiffmat21 != NULL) stiffmat21->Shape(n_dof_element_2_, n_dof_element_1_);
    if (stiffmat22 != NULL) stiffmat22->Shape(n_dof_element_2_, n_dof_element_2_);

    if (forcevec1 != NULL && forcevec2 != NULL)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < n_dof_element_1_; i_dof++)
        (*forcevec1)(i_dof) = force_element_1(i_dof).val();
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < n_dof_element_2_; i_dof++)
        (*forcevec2)(i_dof) = force_element_2(i_dof).val();
    }

    if (stiffmat11 != NULL && stiffmat12 != NULL && stiffmat21 != NULL && stiffmat22 != NULL)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < n_dof_element_1_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < n_dof_element_1_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(i_dof_2);

      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < n_dof_element_1_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < n_dof_element_2_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) =
              -force_element_1(i_dof_1).dx(n_dof_element_1_ + i_dof_2);
          (*stiffmat21)(i_dof_2, i_dof_1) = -force_element_2(i_dof_2).dx(i_dof_1);
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < n_dof_element_2_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < n_dof_element_2_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) =
              -force_element_2(i_dof_1).dx(n_dof_element_1_ + i_dof_2);
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
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::ResetState(const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < n_dof_element_1_; i++)
  {
    ele1pos_(i) =
        TYPE_BTS_VMT_AD(n_dof_element_1_ + n_dof_element_2_, i, beam_centerline_dofvec[i]);
  }

  // Solid element.
  for (unsigned int i = 0; i < n_dof_element_2_; i++)
  {
    ele2pos_(i) = TYPE_BTS_VMT_AD(
        n_dof_element_1_ + n_dof_element_2_, n_dof_element_1_ + i, solid_nodal_dofvec[i]);
  }
}


/**
 *
 */
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Print(
    std::ostream& out) const
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
template <unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes,
    numnodalvalues>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

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
    out << ", Gauss points = " << line_to_volume_segments_[index_segment].GetNumerOfGaussPoints();
    out << "\n";
  }
}


/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<8, 2, 2>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<20, 2, 2>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<27, 2, 2>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<4, 2, 2>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<10, 2, 2>;
