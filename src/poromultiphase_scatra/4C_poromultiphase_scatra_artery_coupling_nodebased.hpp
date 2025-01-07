// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FSI
{
  class Monolithic;
}  // namespace FSI


namespace PoroMultiPhaseScaTra
{
  //! Node based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeBased : public PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    PoroMultiPhaseScaTraArtCouplNodeBased(std::shared_ptr<Core::FE::Discretization> arterydis,
        std::shared_ptr<Core::FE::Discretization> contdis,
        const Teuchos::ParameterList& meshtyingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname);

    //! Evaluate the 1D-3D coupling
    void evaluate(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        std::shared_ptr<Core::LinAlg::Vector<double>> rhs) override
    {
      // nothing to do here, is done in SetupSystem for this type of coupling
    }

    //! set-up linear system of equations of coupled problem
    void setup_system(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_cont,
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_art,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_art,
        std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_cont,
        std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_art) override;

    /*!
     * @brief setup global vector
     *
     * @param[out]  vec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void setup_vector(std::shared_ptr<Core::LinAlg::Vector<double>> vec,
        std::shared_ptr<const Core::LinAlg::Vector<double>> vec_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>> vec_art) override;

    /*!
     * @brief extract single field vectors
     *
     * @param[out]  globalvec combined vector containing both artery and continuous field quantities
     * @param[in]   vec_cont vector containing quantities from continuous field
     * @param[in]   vec_art vector containing quantities from artery field
     */
    void extract_single_field_vectors(std::shared_ptr<const Core::LinAlg::Vector<double>> globalvec,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& vec_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& vec_art) override;

    //! check if initial fields on coupled DOFs are equal
    void check_initial_fields(std::shared_ptr<const Core::LinAlg::Vector<double>> vec_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>> vec_art) override;

    //! access artery (1D) dof row map
    std::shared_ptr<const Epetra_Map> artery_dof_row_map() const override;

    //! access full dof row map
    std::shared_ptr<const Epetra_Map> dof_row_map() const override;

    //! init the strategy
    void init() override;

    //! setup the strategy
    void setup() override;

    //! apply mesh movement (on artery elements)
    void apply_mesh_movement() override;

    //! access to blood vessel volume fraction
    std::shared_ptr<const Core::LinAlg::Vector<double>> blood_vessel_volume_fraction() override;

    //! print out the coupling method
    void print_out_coupling_method() const override;

   private:
    //! set-up of global rhs vector of coupled problem
    void setup_rhs(std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_art);

    //! set-up of global matrix of coupled problem
    void setup_matrix(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_cont,
        Core::LinAlg::SparseMatrix& sysmat_art);

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
    void setup_map_extractor(Core::LinAlg::MultiMapExtractor& mapextractor,
        Core::FE::Discretization& dis, const std::vector<int>& coupleddofs);

    /*!
     * @brief check if dirichlet BC is defined on coupled dofs, which is not possible
     *
     * @param[in]   dis discretizatiom
     * @param[in]   coupleddofmap map with coupled DOFs
     */
    void check_dbc_on_coupled_dofs(
        Core::FE::Discretization& dis, const std::shared_ptr<const Epetra_Map>& coupleddofmap);

    //! name of the condition
    const std::string condname_;

    //! dof row map split in (field) blocks
    std::shared_ptr<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! extractors for continuous field and artery field, maps(0) -> Coupled Dofs, maps(1) uncoupled
    //! Dofs
    std::shared_ptr<Core::LinAlg::MultiMapExtractor> contfieldex_;
    std::shared_ptr<Core::LinAlg::MultiMapExtractor> artex_;

    //! coupling adapter
    std::shared_ptr<Coupling::Adapter::Coupling> artcontfieldcoup_;

    //! needed for matrix transforms
    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> sbbtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> sbitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> sibtransform_;
  };

}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
