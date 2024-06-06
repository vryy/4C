/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for surface-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_SURFBASED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_SURFBASED_HPP


#include "4C_config.hpp"

#include "4C_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroMultiPhaseScaTra
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplSurfBased : public PoroMultiPhaseScaTraArtCouplNonConforming
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplSurfBased(Teuchos::RCP<Discret::Discretization> arterydis,
        Teuchos::RCP<Discret::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

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
    void Evaluate(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override;

   private:
    //! pre-evaluate the pairs
    void pre_evaluate_coupling_pairs();

    //! calculate the volume fraction occupied by blood vessels
    void calculate_blood_vessel_volume_fraction();

    /*!
     * @brief et the artery diameter in material to be able to use it on 1D discretization
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void set_artery_diam_in_material() override{};

    /*!
     * @brief reset the integrated diameter vector to zero
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void reset_integrated_diam_to_zero() override{};

    /*!
     * @brief evaluate additional linearization of (integrated) element diameter dependent terms
     * (Hagen-Poiseuille)
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void evaluate_additional_linearizationof_integrated_diam() override{};

    //! get the segment lengths of element 'artelegid'
    std::vector<double> get_ele_segment_lengths(const int artelegid) override { return {2.0}; };

    //! print out the coupling method
    void print_out_coupling_method() const override;

    /*!
     * @brief check if pair is not active
     * @param[in] coupling_pair: the coupling pair which is checked
     * @return true if pair is not active, false otherwise
     */
    static bool is_not_active(
        const Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScatraArteryCouplingPairBase>
            coupling_pair);
  };
}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
