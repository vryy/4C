/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of OD blocks for monolithic STI

\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/
#ifndef STI_MONOLITHIC_EVALUATE_OD_H
#define STI_MONOLITHIC_EVALUATE_OD_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
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

namespace STI
{
  //! base class for evaluation of scatra-thermo OD blocks
  class ScatraThermoODCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoODCoupling(Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    virtual ~ScatraThermoODCoupling() = default;

    //! evaluation of domain contributions to scatra-thermo OD block
    void EvaluateODBlockScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblock);

    //! evaluation of interface contributions to scatra-thermo OD block
    virtual void EvaluateODBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) = 0;

    //! evaluation of domain contributions to thermo-scatra OD block
    void EvaluateODBlockThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablock);

    //! evaluation of interface contributions to thermo-scatra OD block
    virtual void EvaluateODBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) = 0;

   protected:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_;

    //! map extractor associated with all degrees of freedom inside temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with all degrees of freedom inside scatra field
    const Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside thermo field
    const Teuchos::RCP<const Epetra_Map> full_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! ScaTra subproblem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! Thermo subproblem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
  };

  //! evaluation of scatra-thermo OD blocks for matching nodes
  class ScatraThermoODCouplingMatchingNodes : public ScatraThermoODCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoODCouplingMatchingNodes(Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    virtual ~ScatraThermoODCouplingMatchingNodes() = default;

    //! evaluation of interface contributions to scatra-thermo OD block
    void EvaluateODBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) override;

    //! evaluation of interface contributions to thermo-scatra OD block
    void EvaluateODBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) override;

   private:
    //! copy slave entries to master side scaled by -1.0
    void CopySlaveToMasterScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator>& mastermatrix);

    //! evaluate condition on slave side
    void EvaluateScatraThermoInterfaceSlaveSide(Teuchos::RCP<LINALG::SparseOperator> slavematrix);
  };

  //! evaluation of scatra-thermo OD blocks for standard Mortar on scatra discretization and
  //! condensed Bubnov Mortar on thermo discretization
  class ScatraThermoODCouplingMortarStandard : public ScatraThermoODCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoODCouplingMortarStandard(Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    virtual ~ScatraThermoODCouplingMortarStandard() = default;

    //! evaluation of interface contributions to scatra-thermo OD block
    void EvaluateODBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) override;

    //! evaluation of interface contributions to thermo-scatra OD block
    void EvaluateODBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) override;
  };

}  // namespace STI
#endif  // #ifndef STI_MONOLITHIC_EVALUATE_OD_H
