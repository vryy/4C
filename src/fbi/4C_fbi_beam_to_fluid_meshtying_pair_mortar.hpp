/*! \file

\brief Meshtying element for meshtying between a 3D beam and a 3D fluid element using mortar shape
functions.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_PAIR_MORTAR_HPP
#define FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_PAIR_MORTAR_HPP


#include "4C_config.hpp"

#include "4C_fbi_beam_to_fluid_meshtying_pair_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to fluid meshtying using mortar shape functions for the contact
   * tractions.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param fluid Type from GEOMETRYPAIR::ElementDiscretization... representing the fluid.
   * @param mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar shape
   * functions.
   */
  template <typename beam, typename fluid, typename mortar>
  class BeamToFluidMeshtyingPairMortar : public BeamToFluidMeshtyingPairBase<beam, fluid>
  {
   private:
    //! Shortcut to base class.
    using base_class = BeamToFluidMeshtyingPairBase<beam, fluid>;

    //! Scalar type for FAD variables.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToFluidMeshtyingPairMortar();



    /**
     * \brief Evaluate the mortar matrices $D$ and $M$ for this meshtying element pair.
     * @param local_D (out) Local mortar matrix $D$.
     * @param local_M (out) Local mortar matrix $M$.
     * @param local_kappa (out) Local scaling vector.
     * @param local_constraint_offset (outl) Local constraint offset vector.
     * @return True if pair is in contact.
     */
    bool EvaluateDM(Core::LinAlg::SerialDenseMatrix& local_D,
        Core::LinAlg::SerialDenseMatrix& local_M, Core::LinAlg::SerialDenseVector& local_kappa,
        Core::LinAlg::SerialDenseVector& local_constraint_offset) override;

    /**
     * \brief This pair enforces constraints via a mortar-type method, which requires an own
     * assembly method (provided by the mortar manager).
     */
    inline bool IsAssemblyDirect() const override { return false; };

    /**
     * \brief Add the visualization of this pair to the vtu output writer. This will
     * add mortar specific data to the output.
     * @param visualization_writer Object that manages all visualization related data for beam
     * to fluid pairs
     * @param visualization_params Parameter list
     */
    void get_pair_visualization(
        Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
        Teuchos::ParameterList& visualization_params) const override;

   protected:
    virtual void evaluate_penalty_force(Core::LinAlg::Matrix<3, 1, scalar_type>& force,
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point,
        Core::LinAlg::Matrix<3, 1, scalar_type> v_beam) const;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
