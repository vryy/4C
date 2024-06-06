/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSTI

\level 2

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP
#define FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
  class CouplingSlaveConverter;
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
  namespace UTILS
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
        Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapstructure,
        Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying,
        Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<Adapter::SSIStructureWrapper> structure,
        Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo);

    //! derivative of structure residuals w.r.t. thermo dofs in domain
    void evaluate_off_diag_block_structure_thermo_domain(
        Teuchos::RCP<Core::LinAlg::SparseOperator> structurethermodomain);

    //! derivative of thermo residuals w.r.t. structure dofs in domain
    void evaluate_off_diag_block_thermo_structure_domain(
        Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructuredomain);

    //! derivative of thermo residuals w.r.t. structure dofs on interface
    void evaluate_off_diag_block_thermo_structure_interface(
        Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructureinterface);

   private:
    void copy_slave_to_master_thermo_structure_interface(
        Teuchos::RCP<const Core::LinAlg::SparseOperator> slavematrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& mastermatrix);

    void evaluate_thermo_structure_interface_slave_side(
        Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix);

    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapstructure_;

    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    Teuchos::RCP<const Epetra_Map> full_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! SSTI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying_;

    //! structure problem
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure_;

    //! thermo problem
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo_;
  };
}  // namespace SSTI

FOUR_C_NAMESPACE_CLOSE

#endif
