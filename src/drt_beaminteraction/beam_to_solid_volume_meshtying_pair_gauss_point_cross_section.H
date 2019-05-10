/*!
\file beam_to_solid_volume_meshtying_pair_gauss_point_cross_section.H

\brief Meshtying element for meshtying between a beam and a 3D solid element using Gauss points
on the surface of the (circular) beam cross section.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"


// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying using Gauss point projection on the (circular) beam
   * cross section.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPairGaussPointCrossSection
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   private:
    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + solid::n_dof_> TYPE_BTS_VMT_AD;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairGaussPointCrossSection();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairGaussPointCrossSection(){};

    /**
     * \brief Initialize this pair. Check that the correct geometry pair is given in the input file.
     */
    void Init(const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
        std::vector<DRT::Element const*> elements) override;

    /**
     * \brief Evaluate this contact element pair.
     * @param forcevec1 (out) Force vector on element 1.
     * @param forcevec2 (out) Force vector on element 2.
     * @param stiffmat11 (out) Stiffness contributions on element 1 - element 1.
     * @param stiffmat12 (out) Stiffness contributions on element 1 - element 2.
     * @param stiffmat21 (out) Stiffness contributions on element 2 - element 1.
     * @param stiffmat22 (out) Stiffness contributions on element 2 - element 2.
     * @return True if pair is in contact.
     */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override;

    /**
     * \brief Create the geometry pair for this contact pair. We overload that function because this
     * pair requires explicitly that a cross section projection pair is created.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
        override;

   protected:
    /**
     * \brief Calculate the position on the beam, also taking into accout parameter coordinates on
     * the cross section.
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    void EvaluateBeamPosition(
        const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& integration_point,
        LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_beam, bool reference) const override;
  };
}  // namespace BEAMINTERACTION

#endif
