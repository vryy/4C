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

namespace POROMULTIPHASESCATRA
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplSurfBased : public PoroMultiPhaseScaTraArtCouplNonConforming
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplSurfBased(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    //! set-up of global system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art) override;

    //! setup the strategy
    void Setup() override;

    //! apply mesh movement (on artery elements)
    void ApplyMeshMovement() override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction() override;

    //! Evaluate the 1D-3D coupling
    void Evaluate(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override;

   private:
    //! pre-evaluate the pairs
    void PreEvaluateCouplingPairs();

    //! calculate the volume fraction occupied by blood vessels
    void CalculateBloodVesselVolumeFraction();

    /*!
     * @brief et the artery diameter in material to be able to use it on 1D discretization
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void SetArteryDiamInMaterial() override{};

    /*!
     * @brief reset the integrated diameter vector to zero
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void ResetIntegratedDiamToZero() override{};

    /*!
     * @brief evaluate additional linearization of (integrated) element diameter dependent terms
     * (Hagen-Poiseuille)
     * \note nothing to do for surface-based formulation since varying diameter not yet possible
     */
    void EvaluateAdditionalLinearizationofIntegratedDiam() override{};

    //! get the segment lengths of element 'artelegid'
    std::vector<double> GetEleSegmentLengths(const int artelegid) override { return {2.0}; };

    //! print out the coupling method
    void PrintOutCouplingMethod() const override;

    /*!
     * @brief check if pair is not active
     * @param[in] coupling_pair: the coupling pair which is checked
     * @return true if pair is not active, false otherwise
     */
    static bool IsNotActive(
        const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
            coupling_pair);
  };
}  // namespace POROMULTIPHASESCATRA


FOUR_C_NAMESPACE_CLOSE

#endif
