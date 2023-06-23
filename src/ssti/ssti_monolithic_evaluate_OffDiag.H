/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSTI

\level 2

 *----------------------------------------------------------------------*/

#ifndef SSTI_MONOLITHIC_EVALUATE_OFFDIAG_H
#define SSTI_MONOLITHIC_EVALUATE_OFFDIAG_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class Coupling;
  class CouplingSlaveConverter;
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
}  // namespace ADAPTER

namespace LINALG
{
  class SparseOperator;
  class MultiMapExtractor;
}  // namespace LINALG

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
        Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapstructure,
        Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! derivative of structure residuals w.r.t. thermo dofs in domain
    void EvaluateOffDiagBlockStructureThermoDomain(
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain);

    //! derivative of thermo residuals w.r.t. structure dofs in domain
    void EvaluateOffDiagBlockThermoStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain);

    //! derivative of thermo residuals w.r.t. structure dofs on interface
    void EvaluateOffDiagBlockThermoStructureInterface(
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface);

   private:
    void CopySlaveToMasterThermoStructureInterface(
        Teuchos::RCP<const LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator>& mastermatrix);

    void EvaluateThermoStructureInterfaceSlaveSide(
        Teuchos::RCP<LINALG::SparseOperator> slavematrix);

    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapstructure_;

    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    Teuchos::RCP<const Epetra_Map> full_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! SSTI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying_;

    //! structure problem
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;

    //! thermo problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
  };
}  // namespace SSTI
#endif  // SSTI_MONOLITHIC_EVALUATE_OFFDIAG_H
