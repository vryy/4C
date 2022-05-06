/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for contact between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_H_
#define BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_H_


#include "beam_to_solid_surface_contact_pair_base.H"

#include "beam3contact_defines.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to surface surface contact based on manual variation of the gap function.
   * @tparam scalar_type Type for scalar DOF values.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceContactPairGapVariation
      : public BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceContactPairGapVariation();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceContactPairGapVariation(){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;
  };

  /**
   * \brief Class for beam to surface surface contact based on variation of the penalty potential.
   * @tparam scalar_type Type for scalar DOF values.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceContactPairPotential
      : public BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceContactPairPotential();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceContactPairPotential(){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;
  };
}  // namespace BEAMINTERACTION

#endif
