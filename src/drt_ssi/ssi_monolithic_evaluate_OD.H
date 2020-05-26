/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of OD blocks for monolithic SSI

\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_EVALUATE_OD_H
#define SSI_MONOLITHIC_EVALUATE_OD_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
}  // namespace ADAPTER

namespace LINALG
{
  class SparseOperator;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
}  // namespace SCATRA

namespace SSI
{
  class ScatraStructureODCoupling
  {
   public:
    //! constructor
    explicit ScatraStructureODCoupling(const Teuchos::RCP<const Epetra_Map>& full_map_scatra,
        const Teuchos::RCP<const Epetra_Map>& full_map_structure,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I>& meshtying_strategy_s2i,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>& scatra,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper>& structure);

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void EvaluateODBlockScatraStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator>& scatrastructureblock);

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void EvaluateODBlockScatraStructureInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrastructureinterface);

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    void EvaluateODBlockStructureScatraDomain(
        Teuchos::RCP<LINALG::SparseOperator>& structurescatradomain) const;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    const Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    const Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! structure problem
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;
  };
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_EVALUATE_OD_H
