/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for node-based coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
}

namespace FSI
{
  class Monolithic;
}  // namespace FSI
namespace Core::LinAlg
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace Core::LinAlg


namespace PoroMultiPhaseScaTra
{
  //! Node based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeBased : public PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplNodeBased(Teuchos::RCP<Discret::Discretization> arterydis,
        Teuchos::RCP<Discret::Discretization> contdis,
        const Teuchos::ParameterList& meshtyingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname);

    //! Evaluate the 1D-3D coupling
    void Evaluate(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override
    {
      // nothing to do here, is done in SetupSystem for this type of coupling
    }

    //! set-up linear system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_cont,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmap_art) override;

    /*!
     * @brief setup global vector
     *
     * @param[out]  vec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void setup_vector(Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override;

    /*!
     * @brief extract single field vectors
     *
     * @param[out]  globalvec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void extract_single_field_vectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont,
        Teuchos::RCP<const Epetra_Vector>& vec_art) override;

    //! check if initial fields on coupled DOFs are equal
    void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override;

    //! access artery (1D) dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    //! access full dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() const override;

    //! init the strategy
    void Init() override;

    //! setup the strategy
    void Setup() override;

    //! apply mesh movement (on artery elements)
    void ApplyMeshMovement() override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> blood_vessel_volume_fraction() override;

    //! print out the coupling method
    void print_out_coupling_method() const override;

   private:
    //! set-up of global rhs vector of coupled problem
    void setup_rhs(Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_cont,
        Teuchos::RCP<const Epetra_Vector> rhs_art);

    //! set-up of global matrix of coupled problem
    void setup_matrix(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_cont,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat_art);

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
    void setup_map_extractor(Teuchos::RCP<Core::LinAlg::MultiMapExtractor> mapextractor,
        Teuchos::RCP<Discret::Discretization> dis, const std::vector<int>& coupleddofs);

    /*!
     * @brief check if dirichlet BC is defined on coupled dofs, which is not possible
     *
     * @param[in]   dis discretizatiom
     * @param[in]   coupleddofmap map with coupled DOFs
     */
    void check_dbc_on_coupled_dofs(Teuchos::RCP<Discret::Discretization> dis,
        const Teuchos::RCP<const Epetra_Map>& coupleddofmap);

    //! name of the condition
    const std::string condname_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! extractors for continous field and artery field, maps(0) -> Coupled Dofs, maps(1) uncoupled
    //! Dofs
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> contfieldex_;
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> artex_;

    //! coupling adapter
    Teuchos::RCP<Core::Adapter::Coupling> artcontfieldcoup_;

    //! needed for matrix transforms
    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform> sbbtransform_;
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform> sbitransform_;
    Teuchos::RCP<Core::LinAlg::MatrixColTransform> sibtransform_;
  };

}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
