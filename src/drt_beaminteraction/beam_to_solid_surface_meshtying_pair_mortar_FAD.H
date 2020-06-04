/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element, coupling terms are
evaluated with FAD.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_FAD_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_FAD_H_


#include "beam_to_solid_surface_meshtying_pair_mortar_base.H"

#include "../drt_geometry_pair/geometry_pair_scalar_types.H"


// Forward declaration.
namespace INPAR
{
  namespace BEAMTOSOLID
  {
    enum class BeamToSolidMortarShapefunctions;
  }
}  // namespace INPAR

namespace BEAMINTERACTION
{
  /**
   * \brief Class for Mortar beam to surface surface mesh tying.
   * @tparam scalar_type Type for scalar variables.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   * @tparam mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar
   * shape functions.
   */
  template <typename scalar_type, typename beam, typename surface, typename mortar>
  class BeamToSolidSurfaceMeshtyingPairMortarFAD
      : public BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface, mortar>
  {
   private:
    //! Shortcut to the base class.
    using base_class =
        BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface, mortar>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairMortarFAD();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairMortarFAD(){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix. (derived)
     */
    void EvaluateAndAssemble(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Epetra_Vector& global_lambda) override;

    /**
     * \brief Evaluate the mortar matrices $D$ and $M$ for this contact element pair and assemble
     * them into global matrices. (derived)
     */
    void EvaluateAndAssembleDM(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_D,
        LINALG::SparseMatrix& global_M, Epetra_FEVector& global_constraint,
        Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active) override;
  };


  /**
   * \brief Factory function for beam-to-solid mortar FAD pairs.
   * @param surface_shape (in) Type of surface element.
   * @param mortar_shapefunction (in) Type of mortar shape function.
   * @return Pointer to the created pair.
   */
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
      const DRT::Element::DiscretizationType surface_shape,
      const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction);
}  // namespace BEAMINTERACTION

#endif
