/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element, coupling terms are
evaluated with FAD.

\level 3
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
namespace LARGEROTATIONS
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}  // namespace LARGEROTATIONS


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
        const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector) override;

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
    //! Number of rotational Lagrange multipliers.
    unsigned int n_mortar_rot_;
  };


  /**
   * \brief Class for beam to solid surface rotational mesh tying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param surface Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   * @param mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar shape
   * functions for displacement coupling.
   *
   * @tparam scalar_type Scalar type for position DOFs.
   * @tparam beam Type of used line element.
   * @tparam scalar_type Type of used surface element.
   * @tparam scalar_type Type of used mortar interpolation.
   */
  template <typename scalar_type, typename beam, typename surface, typename mortar>
  class BeamToSolidSurfaceMeshtyingPairMortarRotationFAD
      : public BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface, mortar>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface, mortar>;

    //! FAD type to evaluate the rotational coupling terms. The first 3 entries are the values of
    //! psi_beam, the following entries are the discrete solid DOFs.
    using scalar_type_rot_1st = typename Sacado::Fad::SLFad<double, 3 + surface::n_dof_>;
    using scalar_type_rot_2nd = typename FADUTILS::HigherOrderFadType<2, scalar_type_rot_1st>::type;

    //! Number of rotational DOF for the SR beams;
    static const unsigned int n_dof_rot_ = 9;
    static const unsigned int n_dof_pair_ = n_dof_rot_ + surface::n_dof_;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairMortarRotationFAD() : base_class()
    {
      this->n_mortar_rot_ = mortar::n_dof_;
    };

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector) override;

    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling. (derived)
     */
    void EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_GB,
        LINALG::SparseMatrix& global_GS, LINALG::SparseMatrix& global_FB,
        LINALG::SparseMatrix& global_FS, Epetra_FEVector& global_constraint,
        Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

   private:
    /**
     * \brief Calculate the solid rotation vector on the surface.
     * @param xi (in) Parameter coordinate in the surface
     * @param q_solid_ref (in) Solid position vector in the reference configuration
     * @param q_solid (in) Solid position vector in the current configuration
     * @param quaternion_beam_ref (in) Reference rotation of the beam.
     * @param psi_solid (out) Rotation vector on solid surface.
     * @param element (in) Pointer to the face element.
     */
    template <typename scalar_type_rot>
    void GetSurfaceRotationVector(const LINALG::Matrix<3, 1, double>& xi,
        const LINALG::Matrix<surface::n_dof_, 1, double>& q_solid_ref,
        const LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot>& q_solid,
        const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
        LINALG::Matrix<3, 1, scalar_type_rot>& psi_solid,
        const DRT::Element* element = nullptr) const;
  };

  /**
   * \brief Factory function for beam-to-solid mortar FAD pairs.
   * @param surface_shape (in) Type of surface element.
   * @param mortar_shapefunction (in) Type of mortar shape function.
   * @return Pointer to the created pair.
   */
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
      const DRT::Element::DiscretizationType surface_shape,
      const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction,
      const bool rotational_coupling);
}  // namespace BEAMINTERACTION

#endif
