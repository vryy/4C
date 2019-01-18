/*!
\file beam_to_solid_volume_meshtying_pair_gauss_point.H

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using Gauss point
projection.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"


// Forward declarations.
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolume;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying using Gauss point projection.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPairGaussPoint
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairGaussPoint();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairGaussPoint(){};

    /**
     * \brief Evaluate this contact element pair.
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override;

   private:
    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + solid::n_dof_> TYPE_BTS_VMT_AD;
  };
}  // namespace BEAMINTERACTION

#endif
