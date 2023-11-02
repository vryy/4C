/*----------------------------------------------------------------------*/
/*! \file

\brief Class for 2D-3D beam-to-solid volume mesh tying based on a plane beam element. This
simplifies the triad construction and torsion free beam elements can be used.

\level 3
*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_PLANE_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_PLANE_H


#include "baci_beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_base.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for 2D-3D beam-to-solid volume mesh tying based on a plane beam element. This
   * simplifies the triad construction and torsion free beam elements can be used.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPair2D3DPlane
      : public BeamToSolidVolumeMeshtyingPair2D3DBase<beam, solid>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPair2D3DBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPair2D3DPlane() = default;


    /*!
     *\brief things that need to be done in a separate loop before the actual evaluation loop
     *      over all contact pairs
     */
    void PreEvaluate() override;

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
    bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22) override;

   protected:
    /**
     * \brief Get the triad of the beam at the parameter coordinate xi (derived)
     */
    void GetTriadAtXiDouble(const double xi, CORE::LINALG::Matrix<3, 3, double>& triad,
        const bool reference) const override;
  };
}  // namespace BEAMINTERACTION

#endif
