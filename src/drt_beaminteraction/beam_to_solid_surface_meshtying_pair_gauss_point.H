/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_H_


#include "beam_to_solid_surface_meshtying_pair_base.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for Gauss-point-to-segment beam to surface surface mesh tying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename beam, typename surface>
  class BeamToSolidSurfaceMeshtyingPairGaussPoint
      : public BeamToSolidSurfaceMeshtyingPairBase<
            Sacado::ELRFad::SLFad<double, beam::n_dof_ + surface::n_dof_>, beam, surface>
  {
   private:
    //! Type to be used for scalar AD variables.
    using scalar_type = Sacado::ELRFad::SLFad<double, beam::n_dof_ + surface::n_dof_>;

    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairGaussPoint();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairGaussPoint(){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix) override;
  };
}  // namespace BEAMINTERACTION

#endif
