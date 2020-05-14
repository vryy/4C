/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingPairMortar_TestSuite;
}


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
    //! Define the unit test class as friend so it can set up a valid pair state for the test cases.
    friend BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar_TestSuite;

   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<beam, solid>;

    //! Type to be used for scalar AD variables.
    using scalar_type = typename base_class::scalar_type;

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
     * @param local_D (out) Local mortar matrix $D$.
     * @param local_M (out) Local mortar matrix $M$.
     * @param local_kappa (out) Local scaling vector.
     * @return True if pair is in contact.
     */
    bool EvaluateDM(LINALG::SerialDenseMatrix& local_D, LINALG::SerialDenseMatrix& local_M,
        LINALG::SerialDenseVector& local_kappa) override;

    /**
     * \brief This pair has no terms that need to be directly assembled into the global matrices.
     */
    inline bool IsAssemblyDirect() const override { return false; };

    /**
     * \brief This pair has terms that need to be indirectly (mortar like) assembled into the global
     * matrices.
     */
    inline bool IsAssemblyIndirect() const override { return true; };

    /**
     * \brief Add the visualization of this pair to the beam to solid vtu output writer. This will
     * add mortar specific data to the output.
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list (not used in this class).
     */
    void GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
        const Teuchos::ParameterList& visualization_params) const override;

   protected:
    /**
     * \brief For the mortar pairs it does not make sense to calculate forces at the integration
     * points.
     * @param r_beam (in) Position on the beam.
     * @param r_solid (in) Position on the solid.
     * @param force (out) Return 0 by default.
     */
    void EvaluatePenaltyForce(const LINALG::Matrix<3, 1, scalar_type>& r_beam,
        const LINALG::Matrix<3, 1, scalar_type>& r_solid,
        LINALG::Matrix<3, 1, scalar_type>& force) const override;
  };
}  // namespace BEAMINTERACTION

#endif
