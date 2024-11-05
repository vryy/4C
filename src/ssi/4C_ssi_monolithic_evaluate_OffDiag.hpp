// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSI_MONOLITHIC_EVALUATE_OFFDIAG_HPP
#define FOUR_C_SSI_MONOLITHIC_EVALUATE_OFFDIAG_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class SSIStructureWrapper;
}  // namespace Adapter

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace ScaTra
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace ScaTra

namespace SSI
{
  namespace Utils
  {
    class SSIMeshTying;
  }

  class ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraStructureOffDiagCoupling(
        std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_structure,
        std::shared_ptr<const Epetra_Map> full_map_structure,
        std::shared_ptr<const SSI::Utils::SSIMeshTying> ssi_structure_meshtying,
        std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra,
        std::shared_ptr<Adapter::SSIStructureWrapper> structure);

    virtual ~ScatraStructureOffDiagCoupling() = default;

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void evaluate_off_diag_block_scatra_structure_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> scatrastructureblock);

    //! evaluation contributions to off-diagonal manifold scatra-structure block of global system
    //! matrix
    virtual void evaluate_off_diag_block_scatra_manifold_structure_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> scatramanifoldstructureblock);

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void evaluate_off_diag_block_scatra_structure_interface(
        Core::LinAlg::SparseOperator& scatrastructureinterface);

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    virtual void evaluate_off_diag_block_structure_scatra_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> structurescatradomain) const;

   protected:
    //! copy slave side symmetric contributions to the scatra structure interface linearization
    //! entries to master side scaled by -1.0
    void copy_slave_to_master_scatra_structure_symmetric_interface_contributions(
        std::shared_ptr<const Core::LinAlg::SparseOperator> slavematrix,
        std::shared_ptr<Core::LinAlg::SparseOperator>& mastermatrix);

    //! evaluate symmetric contributions to the scatra structure interface linearization on the
    //! slave side
    void evaluate_scatra_structure_symmetric_interface_contributions_slave_side(
        std::shared_ptr<Core::LinAlg::SparseOperator> slavematrix);

    //! evaluate non-symmetric contributions to the scatra structure interface linearization
    void evaluate_scatra_structure_non_symmetric_interface_contributions_slave_side(
        std::shared_ptr<Core::LinAlg::SparseOperator> slavematrix,
        std::shared_ptr<Core::LinAlg::SparseOperator> mastermatrix);

    //! map extractor associated with all degrees of freedom inside structural field
    std::shared_ptr<const Epetra_Map> full_map_structure() const { return full_map_structure_; }

    //! scatra discretization
    std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra_field() const { return scatra_; }

   private:
    //! map extractor associated with all degrees of freedom inside structure field
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside structural field
    std::shared_ptr<const Epetra_Map> full_map_structure_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! scatra discretization
    std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra_;

    //! structure problem
    std::shared_ptr<Adapter::SSIStructureWrapper> structure_;

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    std::shared_ptr<const SSI::Utils::SSIMeshTying> ssi_structure_meshtying_;
  };

  class ScatraManifoldStructureOffDiagCoupling : public ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraManifoldStructureOffDiagCoupling(
        std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_structure,
        std::shared_ptr<const Epetra_Map> full_map_structure,
        std::shared_ptr<const SSI::Utils::SSIMeshTying> ssi_structure_meshtying,
        std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra,
        std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra_manifold,
        std::shared_ptr<Adapter::SSIStructureWrapper> structure);

    void evaluate_off_diag_block_scatra_manifold_structure_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> scatramanifoldstructureblock) override;

   private:
    //! scatra manifold discretization
    std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra_manifold_;
  };

  /*!
   * fixme: This class is only introduced since the ssti framework is not yet restructured as the
   * ssi framework in the sense that e.g. the mesh tying contributions are still added within the
   * assembly. Once this is changed the structure can be adapted according to the ssi framework and
   * this class is redundant.
   */
  class ScatraStructureOffDiagCouplingSSTI : public ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraStructureOffDiagCouplingSSTI(
        std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_structure,
        std::shared_ptr<const Epetra_Map> full_map_scatra,
        std::shared_ptr<const Epetra_Map> full_map_structure,
        std::shared_ptr<const SSI::Utils::SSIMeshTying> ssi_structure_meshtying,
        std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra,
        std::shared_ptr<Adapter::SSIStructureWrapper> structure);

    void evaluate_off_diag_block_structure_scatra_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> structurescatradomain) const override;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    std::shared_ptr<const Epetra_Map> full_map_scatra_;
  };
}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
