/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for coupling 1D node to point in 3D (non-conforming) between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

\level 3

    *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODETOPOINT_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODETOPOINT_HPP


#include "baci_config.hpp"

#include "baci_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROMULTIPHASESCATRA
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeToPoint : public PoroMultiPhaseScaTraArtCouplNonConforming
  {
   public:
    //! constructor
    PoroMultiPhaseScaTraArtCouplNodeToPoint(Teuchos::RCP<DRT::Discretization> arterydis,
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

    /*!
     * @brief set the artery diameter in material to be able to use it on 1D discretization
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void SetArteryDiamInMaterial() override
    {
      dserror("Function 'SetArteryDiamInMaterial()' not possible for node-to-point coupling");
    };

    /*!
     * @brief reset the integrated diameter vector to zero
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void ResetIntegratedDiamToZero() override
    {
      dserror("Function 'ResetIntegratedDiamToZero()' not possible for node-to-point coupling");
    };

    /*!
     * @brief evaluate additional linearization of (integrated) element diameter dependent terms
     * (Hagen-Poiseuille)
     * \note not possible for node-to-point formulation since varying diameter not yet possible
     */
    void EvaluateAdditionalLinearizationofIntegratedDiam() override
    {
      dserror(
          "Function 'EvaluateAdditionalLinearizationofIntegratedDiam()' not possible for "
          "node-to-point coupling");
    };

    /*!
     * @brief get the segment lengths of element 'artelegid'
     * \note segment length is set to zero since we have no segments in node-to-point coupling
     */
    std::vector<double> GetEleSegmentLengths(const int artelegid) override { return {0.0}; };

   private:
    //! print out the coupling method
    void PrintOutCouplingMethod() const override;

    //! preevaluate the coupling pairs
    void PreEvaluateCouplingPairs();

    //! Output Coupling pairs
    void OutputCouplingPairs() const;
  };
}  // namespace POROMULTIPHASESCATRA


FOUR_C_NAMESPACE_CLOSE

#endif