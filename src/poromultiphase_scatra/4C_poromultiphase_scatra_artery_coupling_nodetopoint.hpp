/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for coupling 1D node to point in 3D (non-conforming) between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

\level 3

    *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODETOPOINT_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODETOPOINT_HPP


#include "4C_config.hpp"

#include "4C_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroMultiPhaseScaTra
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeToPoint : public PoroMultiPhaseScaTraArtCouplNonConforming
  {
   public:
    //! constructor
    PoroMultiPhaseScaTraArtCouplNodeToPoint(Teuchos::RCP<Core::FE::Discretization> arterydis,
        Teuchos::RCP<Core::FE::Discretization> contdis,
        const Teuchos::ParameterList& couplingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname);

    //! set-up of global system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_cont,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmap_art) override;

    //! setup the strategy
    void Setup() override;

    //! apply mesh movement (on artery elements)
    void ApplyMeshMovement() override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> blood_vessel_volume_fraction() override;

    //! Evaluate the 1D-3D coupling
    void evaluate(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override;

    /*!
     * @brief set the artery diameter in material to be able to use it on 1D discretization
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void set_artery_diam_in_material() override
    {
      FOUR_C_THROW(
          "Function 'set_artery_diam_in_material()' not possible for node-to-point coupling");
    };

    /*!
     * @brief reset the integrated diameter vector to zero
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void reset_integrated_diam_to_zero() override
    {
      FOUR_C_THROW(
          "Function 'reset_integrated_diam_to_zero()' not possible for node-to-point coupling");
    };

    /*!
     * @brief evaluate additional linearization of (integrated) element diameter dependent terms
     * (Hagen-Poiseuille)
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void evaluate_additional_linearizationof_integrated_diam() override
    {
      FOUR_C_THROW(
          "Function 'evaluate_additional_linearizationof_integrated_diam()' not possible for "
          "node-to-point coupling");
    };

    /*!
     * @brief get the segment lengths of element 'artelegid'
     * \note segment length is set to zero since we have no segments in node-to-point coupling
     */
    std::vector<double> get_ele_segment_lengths(const int artelegid) override { return {0.0}; };

   private:
    //! print out the coupling method
    void print_out_coupling_method() const override;

    //! preevaluate the coupling pairs
    void pre_evaluate_coupling_pairs();

    //! Output Coupling pairs
    void output_coupling_pairs() const;
  };
}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif