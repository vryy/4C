// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP
#define FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
}  // namespace Adapter

namespace Core::LinAlg
{
  class SparseOperator;
  class MultiMapExtractor;
}  // namespace Core::LinAlg

namespace ScaTra
{
  class MeshtyingStrategyS2I;
}

namespace SSI
{
  namespace Utils
  {
    class SSIMeshTying;
  }
}  // namespace SSI

namespace SSTI
{
  class ThermoStructureOffDiagCoupling
  {
   public:
    //! constructor
    explicit ThermoStructureOffDiagCoupling(
        std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmapstructure,
        std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmapthermo,
        std::shared_ptr<const Epetra_Map> full_map_structure,
        std::shared_ptr<const Epetra_Map> full_map_thermo,
        std::shared_ptr<const SSI::Utils::SSIMeshTying> ssti_structure_meshtying,
        std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        std::shared_ptr<Adapter::SSIStructureWrapper> structure,
        std::shared_ptr<Adapter::ScaTraBaseAlgorithm> thermo);

    //! derivative of structure residuals w.r.t. thermo dofs in domain
    void evaluate_off_diag_block_structure_thermo_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> structurethermodomain);

    //! derivative of thermo residuals w.r.t. structure dofs in domain
    void evaluate_off_diag_block_thermo_structure_domain(
        std::shared_ptr<Core::LinAlg::SparseOperator> thermostructuredomain);

    //! derivative of thermo residuals w.r.t. structure dofs on interface
    void evaluate_off_diag_block_thermo_structure_interface(
        Core::LinAlg::SparseOperator& thermostructureinterface);

   private:
    void copy_slave_to_master_thermo_structure_interface(
        std::shared_ptr<const Core::LinAlg::SparseOperator> slavematrix,
        std::shared_ptr<Core::LinAlg::SparseOperator>& mastermatrix);

    void evaluate_thermo_structure_interface_slave_side(
        std::shared_ptr<Core::LinAlg::SparseOperator> slavematrix);

    //! map extractor associated with all degrees of freedom inside structure field
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmapstructure_;

    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with all degrees of freedom inside structural field
    std::shared_ptr<const Epetra_Map> full_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    std::shared_ptr<const Epetra_Map> full_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    std::shared_ptr<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! SSTI structure meshtying object containing coupling adapters, converters and maps
    std::shared_ptr<const SSI::Utils::SSIMeshTying> ssti_structure_meshtying_;

    //! structure problem
    std::shared_ptr<Adapter::SSIStructureWrapper> structure_;

    //! thermo problem
    std::shared_ptr<Adapter::ScaTraBaseAlgorithm> thermo_;
  };
}  // namespace SSTI

FOUR_C_NAMESPACE_CLOSE

#endif
