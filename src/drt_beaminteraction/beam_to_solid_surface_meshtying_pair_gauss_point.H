/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_H_


#include "beam_to_solid_surface_meshtying_pair_base.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for Gauss-point-to-segment beam to surface surface mesh tying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename beam, typename surface>
  class BeamToSolidSurfaceMeshtyingPairGaussPoint
      : public BeamToSolidSurfaceMeshtyingPairBase<
            Sacado::ELRFad::SLFad<double, beam::n_dof_ + surface::n_dof_>, beam, surface>
  {
   private:
    //! Type to be used for scalar AD variables.
    using scalar_type = Sacado::ELRFad::SLFad<double, beam::n_dof_ + surface::n_dof_>;

    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairGaussPoint();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairGaussPoint(){};

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
  };
}  // namespace BEAMINTERACTION

#endif
