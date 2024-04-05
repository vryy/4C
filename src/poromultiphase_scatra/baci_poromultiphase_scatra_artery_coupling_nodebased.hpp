/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for node-based coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP

#include "baci_config.hpp"

#include "baci_coupling_adapter.hpp"
#include "baci_poromultiphase_scatra_artery_coupling_base.hpp"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  class Monolithic;
}  // namespace FSI
namespace CORE::LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace CORE::LINALG


namespace POROMULTIPHASESCATRA
{
  //! Node based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeBased : public PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplNodeBased(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& meshtyingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    //! Evaluate the 1D-3D coupling
    void Evaluate(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override
    {
      // nothing to do here, is done in SetupSystem for this type of coupling
    }

    //! set-up linear system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art) override;

    /*!
     * @brief setup global vector
     *
     * @param[out]  vec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void SetupVector(Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override;

    /*!
     * @brief extract single field vectors
     *
     * @param[out]  globalvec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void ExtractSingleFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont,
        Teuchos::RCP<const Epetra_Vector>& vec_art) override;

    //! check if initial fields on coupled DOFs are equal
    void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override;

    //! access artery (1D) dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    //! access full dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const override;

    //! init the strategy
    void Init() override;

    //! setup the strategy
    void Setup() override;

    //! apply mesh movement (on artery elements)
    void ApplyMeshMovement() override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction() override;

    //! print out the coupling method
    void PrintOutCouplingMethod() const override;

   private:
    //! set-up of global rhs vector of coupled problem
    void SetupRHS(Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_cont,
        Teuchos::RCP<const Epetra_Vector> rhs_art);

    //! set-up of global matrix of coupled problem
    void SetupMatrix(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art);

    /*!
     * @brief setup map extractor for artery mesh tying
     *
     * full map -> all DOFs
     * maps(0)  -> coupled DOFs
     * maps(1)  -> uncoupled DOFs
     *
     * @param[in]   mapextractor the map extractor to setup
     * @param[in]   dis discretization
     * @param[in]   coupleddofs vector with DOFs to couple
     */
    void SetupMapExtractor(Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mapextractor,
        Teuchos::RCP<DRT::Discretization> dis, const std::vector<int>& coupleddofs);

    /*!
     * @brief check if dirichlet BC is defined on coupled dofs, which is not possible
     *
     * @param[in]   dis discretizatiom
     * @param[in]   coupleddofmap map with coupled DOFs
     */
    void CheckDbcOnCoupledDofs(
        Teuchos::RCP<DRT::Discretization> dis, const Teuchos::RCP<const Epetra_Map>& coupleddofmap);

    //! name of the condition
    const std::string condname_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_;

    //! extractors for continous field and artery field, maps(0) -> Coupled Dofs, maps(1) uncoupled
    //! Dofs
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> contfieldex_;
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> artex_;

    //! coupling adapter
    Teuchos::RCP<CORE::ADAPTER::Coupling> artcontfieldcoup_;

    //! needed for matrix transforms
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> sbbtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> sbitransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> sibtransform_;
  };

}  // namespace POROMULTIPHASESCATRA


BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_H
