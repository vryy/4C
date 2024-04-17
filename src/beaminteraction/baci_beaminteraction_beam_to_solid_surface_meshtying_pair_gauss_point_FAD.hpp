/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element. The
coupling terms are evaluated using FAD.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_FAD_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_FAD_HPP


#include "baci_config.hpp"

#include "baci_beaminteraction_beam_to_solid_surface_meshtying_pair_gauss_point_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BEAMINTERACTION
{
  /**
   * \brief Class for Gauss-point-to-segment beam to surface surface mesh tying.
   * @tparam scalar_type Scalar type to be used for FAD evaluation.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceMeshtyingPairGaussPointFAD
      : public BeamToSolidSurfaceMeshtyingPairGaussPointBase<scalar_type, beam, surface>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairGaussPointBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairGaussPointFAD();


    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<CORE::LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
