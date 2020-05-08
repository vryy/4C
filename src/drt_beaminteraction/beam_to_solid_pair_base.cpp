/*----------------------------------------------------------------------*/
/*! \file

\brief Base element for interactions between a beam and a solid.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_pair_base.H"

#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"

#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"

#include "Sacado.hpp"


/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::BeamToSolidPairBase()
    : BeamContactPair(), line_to_3D_segments_()
{
}


/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::Setup()
{
  CheckInit();

  // Call setup of base class first.
  BeamContactPair::Setup();

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
              "Hermite centerline");
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

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1pos_(i) = 0.0;

  issetup_ = true;
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < beam::n_dof_; i++)
  {
    ele1pos_(i) = scalar_type_fad(beam::n_dof_ + solid::n_dof_, i, beam_centerline_dofvec[i]);
  }
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::SetRestartDisplacement(
    const std::vector<std::vector<double>>& centerline_restart_vec_)
{
  // Call the parent method.
  BeamContactPair::SetRestartDisplacement(centerline_restart_vec_);
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::Print(
    std::ostream& out) const
{
  CheckInitSetup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToSolidPairBase"
      << "\nBeam EleGID:  " << Element1()->Id() << "\nSolid EleGID: " << Element2()->Id();

  out << "\n\nbeam dofvec: " << ele1pos_;
  out << "\nn_segments: " << line_to_3D_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam,
    solid>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Only display information if a segment exists for this pair.
  if (line_to_3D_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << Element1()->Id() << ", solid ID " << Element2()->Id() << ":";
  out << " n_segments = " << line_to_3D_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < line_to_3D_segments_.size(); index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << line_to_3D_segments_[index_segment].GetEtaA() << ", "
        << line_to_3D_segments_[index_segment].GetEtaB() << "]";
    out << ", Gauss points = " << line_to_3D_segments_[index_segment].GetNumberOfProjectionPoints();
    out << "\n";
  }
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type_fad, beam, solid>::EvaluateBeamPosition(
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
    LINALG::Matrix<3, 1, scalar_type_fad>& r_beam, bool reference) const
{
  if (reference)
    GEOMETRYPAIR::EvaluatePosition<beam>(
        integration_point.GetEta(), ele1posref_, r_beam, this->Element1());
  else
    GEOMETRYPAIR::EvaluatePosition<beam>(
        integration_point.GetEta(), ele1pos_, r_beam, this->Element1());
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_hex8::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_hex20::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_hex27::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tet4::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tet10::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double,
        GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs27::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>;

template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tri3::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tri6::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad4::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad8::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>;

template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::DFad<Sacado::ELRFad::DFad<double>>, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>;
template class BEAMINTERACTION::BeamToSolidPairBase<
    Sacado::ELRFad::SLFad<Sacado::ELRFad::SLFad<double,
                              GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
        GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>;
