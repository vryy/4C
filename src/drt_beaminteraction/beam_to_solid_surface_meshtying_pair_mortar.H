/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_surface_meshtying_pair_mortar_base.H"

#include "geometry_pair_scalar_types.H"


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
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   * @tparam mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar
   * shape functions.
   */
  template <typename beam, typename surface, typename mortar>
  class BeamToSolidSurfaceMeshtyingPairMortar
      : public BeamToSolidSurfaceMeshtyingPairMortarBase<
            GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>, beam, surface, mortar>
  {
   private:
    //! Type to be used for scalar AD variables.
    using scalar_type = GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>;

    //! Shortcut to the base class.
    using base_class =
        BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface, mortar>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairMortar();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairMortar(){};

    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling. (derived)
     */
    void EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_G_B,
        LINALG::SparseMatrix& global_G_S, LINALG::SparseMatrix& global_FB_L,
        LINALG::SparseMatrix& global_FS_L, Epetra_FEVector& global_constraint,
        Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

   protected:
    /**
     * \brief Evaluate the local mortar matrices for this contact element pair.
     */
    void EvaluateDM(LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
        LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double>& local_M,
        LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
        LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint) const;
  };

  /**
   * \brief Factory function for beam-to-solid mortar pairs.
   * @param surface_shape (in) Type of surface element.
   * @param mortar_shapefunction (in) Type of mortar shape function.
   * @return Pointer to the created pair.
   */
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BeamToSolidSurfaceMeshtyingPairMortarFactory(
      const DRT::Element::DiscretizationType surface_shape,
      const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction);
}  // namespace BEAMINTERACTION

#endif
