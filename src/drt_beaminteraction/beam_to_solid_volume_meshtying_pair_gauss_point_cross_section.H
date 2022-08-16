/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a beam and a 3D solid element using Gauss points
on the surface of the (circular) beam cross section.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_H_


#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section_base.H"


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
      : public BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionBase<beam, solid>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairGaussPointCrossSection();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairGaussPointCrossSection(){};

    /*!
     *\brief things that need to be done in a separate loop before the actual evaluation loop
     *      over all contact pairs
     */
    virtual void PreEvaluate() override;

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

   protected:
    /**
     * \brief Get the triad of the beam at the parameter coordinate xi (derived)
     */
    void GetTriadAtXiDouble(
        const double xi, LINALG::Matrix<3, 3, double>& triad, const bool reference) const override;
  };
}  // namespace BEAMINTERACTION

#endif
