/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSTI

\level 2

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP
#define FOUR_C_SSTI_MONOLITHIC_EVALUATE_OFFDIAG_HPP

#include "baci_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class Coupling;
  class CouplingSlaveConverter;
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
}  // namespace ADAPTER

namespace CORE::LINALG
{
  class SparseOperator;
  class MultiMapExtractor;
}  // namespace CORE::LINALG

namespace SCATRA
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
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmapstructure,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! derivative of structure residuals w.r.t. thermo dofs in domain
    void EvaluateOffDiagBlockStructureThermoDomain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structurethermodomain);

    //! derivative of thermo residuals w.r.t. structure dofs in domain
    void EvaluateOffDiagBlockThermoStructureDomain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> thermostructuredomain);

    //! derivative of thermo residuals w.r.t. structure dofs on interface
    void EvaluateOffDiagBlockThermoStructureInterface(
        Teuchos::RCP<CORE::LINALG::SparseOperator> thermostructureinterface);

   private:
    void CopySlaveToMasterThermoStructureInterface(
        Teuchos::RCP<const CORE::LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& mastermatrix);

    void EvaluateThermoStructureInterfaceSlaveSide(
        Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix);

    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmapstructure_;

    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    Teuchos::RCP<const Epetra_Map> full_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! SSTI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying_;

    //! structure problem
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure_;

    //! thermo problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
  };
}  // namespace SSTI

FOUR_C_NAMESPACE_CLOSE

#endif
