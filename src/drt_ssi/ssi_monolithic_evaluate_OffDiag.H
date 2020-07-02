/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSI

\level 2


 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_EVALUATE_OffDiag_H
#define SSI_MONOLITHIC_EVALUATE_OffDiag_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
  class Coupling;
}  // namespace ADAPTER

namespace LINALG
{
  class SparseOperator;
  class MultiMapExtractor;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
}  // namespace SCATRA

namespace SSI
{
  class ScatraStructureOffDiagCoupling
  {
   public:
    //! constructor
    explicit ScatraStructureOffDiagCoupling(
        Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_scatra,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const ADAPTER::Coupling> interface_coupling_structure,
        Teuchos::RCP<const Epetra_Map> interface_map_scatra,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure);

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void EvaluateOffDiagBlockScatraStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator>& scatrastructureblock);

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void EvaluateOffDiagBlockScatraStructureInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrastructureinterface);

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    void EvaluateOffDiagBlockStructureScatraDomain(
        Teuchos::RCP<LINALG::SparseOperator>& structurescatradomain) const;

   private:
    //! copy slave entries to master side scaled by -1.0
    void CopySlaveToMasterScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator>& mastermatrix);

    //! copy slave entries to master side scaled by -1.0
    void EvaluateScatraStructureInterfaceSlaveSide(
        Teuchos::RCP<LINALG::SparseOperator> slavematrix);

    //! map extractor associated with all degrees of freedom inside structure field
    const Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside scatra field
    const Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    const Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! scatra-scatra interface coupling adapter for structure field
    const Teuchos::RCP<const ADAPTER::Coupling> interface_coupling_structure_;

    //! map associated with all degrees of freedom on scatra interface
    const Teuchos::RCP<const Epetra_Map> interface_map_scatra_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    const Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! scatra problem
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! structure problem
    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;
  };
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_EVALUATE_OffDiag_H
