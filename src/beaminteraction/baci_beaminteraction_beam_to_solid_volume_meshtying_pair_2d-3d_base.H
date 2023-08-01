/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for 2D-3D beam-to-solid volume mesh tying.

\level 3
*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_BASE_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_BASE_H


#include "baci_beaminteraction_beam_to_solid_volume_meshtying_pair_base.H"

#include "baci_beam3_triad_interpolation_local_rotation_vectors.H"

// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Base class for 2D-3D beam to solid volume mesh tying
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPair2D3DBase
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPair2D3DBase() = default;

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPair2D3DBase(){};

    /**
     * \brief Create the geometry pair for this contact pair. We overload that function because this
     * pair requires explicitly that a cross section projection pair is created.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
        override;

   protected:
    /**
     * \brief Calculate the position on the beam, also taking into account parameter coordinates on
     * the cross section.
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    void EvaluateBeamPositionDouble(
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
        CORE::LINALG::Matrix<3, 1, double>& r_beam, bool reference) const override;

    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * @return RPC with the type of geometry pair for this beam contact pair.
     */
    inline Teuchos::RCP<
        GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, beam, solid>>
    CastGeometryPair() const
    {
      return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::
              GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, beam, solid>>(
          this->geometry_pair_, true);
    };

    /**
     * \brief Get the triad of the beam at the parameter coordinate xi
     * @param xi (in) Parameter coordinate on the beam
     * @param triad (out) Beam cross section triad
     * @param reference (in) If the triad in the reference or current configuration should be
     * returned
     */
    virtual void GetTriadAtXiDouble(
        const double xi, CORE::LINALG::Matrix<3, 3, double>& triad, const bool reference) const = 0;
  };
}  // namespace BEAMINTERACTION

#endif
