/*!
\file beam_to_solid_volume_meshtying_pair_mortar.H

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying using mortar shape functions for the contact
   * tractions.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   * @param mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar shape
   * functions.
   */
  template <typename beam, typename solid, typename mortar>
  class BeamToSolidVolumeMeshtyingPairMortar
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairMortar();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairMortar(){};

    /**
     * \brief Evaluate the mortar matrices $D$ and $M$ for this contact element pair.
     */
    void EvaluateDM(LINALG::SerialDenseMatrix* mortar_D, LINALG::SerialDenseMatrix* mortar_M);

   private:
    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + solid::n_dof_> TYPE_BTS_VMT_AD;
  };
}  // namespace BEAMINTERACTION

#endif
