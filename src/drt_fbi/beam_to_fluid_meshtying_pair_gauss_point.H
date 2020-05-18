/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a 1D beam and a 3D fluid element using Gauss point
projection.

\level 1
\maintainer Nora Hagmeyer
*/

#ifndef BEAM_TO_FLUID_MESHTYING_PAIR_GAUSS_POINT_H_
#define BEAM_TO_FLUID_MESHTYING_PAIR_GAUSS_POINT_H_


#include "beam_to_fluid_meshtying_pair_base.H"

namespace FBI
{
  class PairFactory;
}
namespace BEAMINTERACTION
{
  class BeamToFluidMeshtyingPairGPTS_TestSuite;

  /**
   * \brief Class for beam to fluid meshtying using Gauss point projection.
   *
   * \param[in] beam Type from GEOMETRYPAIR::ElementDiscretization representing the beam.
   * \param[in] volume Type from GEOMETRYPAIR::ElementDiscretization... representing the fluid.
   */
  template <typename beam, typename fluid>
  class BeamToFluidMeshtyingPairGaussPoint : public BeamToFluidMeshtyingPairBase<beam, fluid>
  {
    friend class FBI::PairFactory;
    friend class BeamToFluidMeshtyingPairGPTS_TestSuite;

   public:
    /**
     * \brief Destructor.
     */
    virtual ~BeamToFluidMeshtyingPairGaussPoint(){};

    /**
     * \brief Evaluate this contact element pair.
     *
     * \params[inout] forcevec1 (out) Force vector on element 1.
     * \params[inout] forcevec2 (out) Force vector on element 2.
     * \params[inout] stiffmat11 (out) Stiffness contributions on element 1 - element 1.
     * \params[inout] stiffmat12 (out) Stiffness contributions on element 1 - element 2.
     * \params[inout] stiffmat21 (out) Stiffness contributions on element 2 - element 1.
     * \params[inout] stiffmat22 (out) Stiffness contributions on element 2 - element 2.
     *
     * \returns True if pair is in contact.
     */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override;

   protected:
    //! Shortcut to base class.
    using base_class = BeamToFluidMeshtyingPairBase<beam, fluid>;

    //! Scalar type for FAD variables.
    using scalar_type = typename base_class::scalar_type;

   protected:
    /** \brief You will have to use the FBI::PairFactory
     *
     */

    BeamToFluidMeshtyingPairGaussPoint();
  };
}  // namespace BEAMINTERACTION

#endif
