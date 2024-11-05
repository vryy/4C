// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    PoroMultiPhaseScaTraArtCouplSurfBased(std::shared_ptr<Core::FE::Discretization> arterydis,
        std::shared_ptr<Core::FE::Discretization> contdis,
        const Teuchos::ParameterList& couplingparams, const std::string& condname,
        const std::string& artcoupleddofname, const std::string& contcoupleddofname);

    //! set-up of global system of equations of coupled problem
    void setup_system(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_cont,
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat_art,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_cont,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_art,
        std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_cont,
        std::shared_ptr<const Core::LinAlg::MapExtractor> dbcmap_art) override;

    //! setup the strategy
    void setup() override;

    //! apply mesh movement (on artery elements)
    void apply_mesh_movement() override;

    //! access to blood vessel volume fraction
    std::shared_ptr<const Core::LinAlg::Vector<double>> blood_vessel_volume_fraction() override;

    //! Evaluate the 1D-3D coupling
    void evaluate(std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        std::shared_ptr<Core::LinAlg::Vector<double>> rhs) override;

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
  };
}  // namespace PoroMultiPhaseScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
