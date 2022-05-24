/*----------------------------------------------------------------------*/
/*! \file

\brief Mesh tying element to couple points of two 3D beam elements together.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_BEAM_POINT_COUPLING_PAIR_H_
#define BEAM_TO_BEAM_POINT_COUPLING_PAIR_H_


#include "beam_contact_pair.H"

#include "Sacado.hpp"


// Forward declarations.
namespace LARGEROTATIONS
{
  template <unsigned int numnodes, typename T>
  class TriadInterpolationLocalRotationVectors;
}  // namespace LARGEROTATIONS


namespace BEAMINTERACTION
{
  /**
   * \brief Class for point-wise beam to beam mesh tying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   */
  template <typename beam>
  class BeamToBeamPointCouplingPair : public BeamContactPair
  {
   protected:
    //! FAD type for rotational coupling. The 6 dependent DOFs are the 3 rotational DOFs of each
    //! beam element.
    using scalar_type_rot = typename Sacado::Fad::SLFad<double, 6>;

    //! FAD type for positional coupling.
    using scalar_type_pos = typename Sacado::Fad::SLFad<double, 2 * beam::n_dof_>;

   public:
    /**
     * \brief Standard Constructor.
     *
     * @param penalty_parameter_rot (in) Penalty parameter for rotational coupling.
     * @param penalty_parameter_pos (in) Penalty parameter for positional coupling.
     * @param pos_in_parameterspace (in) Coupling positions in the beam parameter spaces.
     */
    BeamToBeamPointCouplingPair(double penalty_parameter_rot, double penalty_parameter_pos,
        std::array<double, 2> pos_in_parameterspace);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToBeamPointCouplingPair(){};

    /**
     * \brief Setup the beam coupling pair.
     */
    void Setup() override;

    /**
     * \brief Things that need to be done in a separate loop before the actual evaluation loop over
     * all contact pairs. (derived)
     */
    void PreEvaluate() override{};

    /**
     * \brief Evaluate this contact element pair.
     */
    bool Evaluate(LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22) override
    {
      return false;
    }

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     *
     * @param discret (in) Pointer to the discretization.
     * @param force_vector (in / out) Global force vector.
     * @param stiffness_matrix (in / out) Global stiffness matrix.
     * @param displacement_vector (in) Global displacement vector.
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

    /**
     * \brief No need to update pair state vectors, as everything is done in the EvaluateAndAssemble
     * call.
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override{};

    /**
     * \brief This pair is always active.
     */
    inline bool GetContactFlag() const override { return true; }

    /**
     * \brief Get number of active contact point pairs on this element pair. Not yet implemented.
     */
    unsigned int GetNumAllActiveContactPointPairs() const override
    {
      dserror("GetNumAllActiveContactPointPairs not yet implemented!");
      return 0;
    };

    /**
     * \brief Get coordinates of all active contact points on element1. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement1(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("GetAllActiveContactPointCoordsElement1 not yet implemented!");
    }

    /**
     * \brief Get coordinates of all active contact points on element2. Not yet implemented.
     */
    void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("GetAllActiveContactPointCoordsElement2 not yet implemented!");
    }

    /**
     * \brief Get all (scalar) contact forces of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactForces(std::vector<double>& forces) const override
    {
      dserror("GetAllActiveContactForces not yet implemented!");
    }

    /**
     * \brief Get all (scalar) gap values of this contact pair. Not yet implemented.
     */
    void GetAllActiveContactGaps(std::vector<double>& gaps) const override
    {
      dserror("GetAllActiveContactGaps not yet implemented!");
    }

    /**
     * \brief Get energy of penalty contact. Not yet implemented.
     */
    double GetEnergy() const override
    {
      dserror("GetEnergy not implemented yet!");
      return 0.0;
    }

    /**
     * \brief Print information about this beam contact element pair to screen.
     */
    void Print(std::ostream& out) const override;

    /**
     * \brief Print this beam contact element pair to screen.
     */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;

   private:
    /**
     * \brief Evaluate the positional coupling terms and directly assemble them into the global
     * force vector and stiffness matrix.
     *
     * @param discret (in) Pointer to the discretization.
     * @param force_vector (in / out) Global force vector.
     * @param stiffness_matrix (in / out) Global stiffness matrix.
     * @param displacement_vector (in) Global displacement vector.
     */
    void EvaluateAndAssemblePositionalCoupling(
        const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) const;

    /**
     * \brief Evaluate the rotational coupling terms and directly assemble them into the global
     * force vector and stiffness matrix.
     *
     * @param discret (in) Pointer to the discretization.
     * @param force_vector (in / out) Global force vector.
     * @param stiffness_matrix (in / out) Global stiffness matrix.
     * @param displacement_vector (in) Global displacement vector.
     */
    void EvaluateAndAssembleRotationalCoupling(
        const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) const;

   private:
    //! Number of rotational DOF for the SR beams;
    static const unsigned int n_dof_rot_ = 9;

    //! Number of dimensions for each rotation.
    const unsigned int rot_dim_ = 3;

    //! Penalty parameter for positional coupling.
    double penalty_parameter_pos_;

    //! Penalty parameter for rotational coupling.
    double penalty_parameter_rot_;

    //! Coupling point positions in the element parameter spaces.
    std::array<double, 2> position_in_parameterspace_;
  };  // namespace BEAMINTERACTION
}  // namespace BEAMINTERACTION

#endif
