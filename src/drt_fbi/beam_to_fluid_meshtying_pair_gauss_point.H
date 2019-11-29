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
   * \params[in] beam Type from GEOMETRYPAIR::ElementDiscretization representing the beam.
   * \params[in] volume Type from GEOMETRYPAIR::ElementDiscretization... representing the fluid.
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

   private:
    /// Type to be used for scalar AD variables. This can not be inherited from the base class.
    typedef Sacado::ELRFad::SLFad<double, beam::n_dof_ + fluid::n_dof_> TYPE_BTS_VMT_AD;

   protected:
    /** \brief You will have to use the FBI::PairFactory
     *
     */
    BeamToFluidMeshtyingPairGaussPoint();

    /**
     * \brief Evaluate the penalty force for the visualization
     */
    virtual void EvaluatePenaltyForce(LINALG::Matrix<3, 1, TYPE_BTS_VMT_AD>& force,
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point,
        LINALG::Matrix<3, 1, TYPE_BTS_VMT_AD> v_beam) const override;
  };
}  // namespace BEAMINTERACTION

#endif
