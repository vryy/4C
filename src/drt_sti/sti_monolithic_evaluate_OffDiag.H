/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic STI

\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/
#ifndef STI_MONOLITHIC_EVALUATE_OffDiag_H
#define STI_MONOLITHIC_EVALUATE_OffDiag_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

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
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace STI
{
  //! base class for evaluation of scatra-thermo off-diagonal blocks
  class ScatraThermoOffDiagCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoOffDiagCoupling(Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        const Teuchos::RCP<const Epetra_Map> interface_map_scatra,
        const Teuchos::RCP<const Epetra_Map> interface_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    virtual ~ScatraThermoOffDiagCoupling() = default;

    //! evaluation of domain contributions to scatra-thermo OD block
    void EvaluateOffDiagBlockScatraThermoDomain(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblock);

    //! evaluation of interface contributions to scatra-thermo OD block
    virtual void EvaluateOffDiagBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) = 0;

    //! evaluation of domain contributions to thermo-scatra OD block
    void EvaluateOffDiagBlockThermoScatraDomain(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablock);

    //! evaluation of interface contributions to thermo-scatra OD block
    virtual void EvaluateOffDiagBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) = 0;

   protected:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> BlockMaps() { return blockmaps_; };

    //! map extractor associated with all degrees of freedom inside temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> BlockMapThermo() { return blockmapthermo_; };

    //! map extractor associated with degrees of freedom on interface of temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> BlockMapThermoInterface()
    {
      return blockmapthermointerface_;
    };

    //! map extractor associated with all degrees of freedom inside scatra field
    const Teuchos::RCP<const Epetra_Map> FullMapScatra() { return full_map_scatra_; };

    //! map extractor associated with all degrees of freedom inside thermo field
    const Teuchos::RCP<const Epetra_Map> FullMapThermo() { return full_map_thermo_; };

    //! map associated with all degrees of freedom on thermo interface
    const Teuchos::RCP<const Epetra_Map> InterfaceMapScaTra() { return interface_map_scatra_; };

    //! map associated with all degrees of freedom on scatra interface
    const Teuchos::RCP<const Epetra_Map> InterfaceMapThermo() { return interface_map_thermo_; };

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> MeshtyingStrategyScaTra()
    {
      return meshtying_strategy_scatra_;
    };

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> MeshtyingStrategyThermo()
    {
      return meshtying_strategy_thermo_;
    };

    //! ScaTra subproblem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_->ScaTraField(); };

    //! Thermo subproblem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ThermoField() { return thermo_->ScaTraField(); };

   private:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_;

    //! map extractor associated with all degrees of freedom inside temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo_;

    //! map extractor associated with degrees of freedom on interface of temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface_;

    //! map extractor associated with all degrees of freedom inside scatra field
    const Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside thermo field
    const Teuchos::RCP<const Epetra_Map> full_map_thermo_;

    //! map associated with all degrees of freedom on thermo interface
    const Teuchos::RCP<const Epetra_Map> interface_map_scatra_;

    //! map associated with all degrees of freedom on scatra interface
    const Teuchos::RCP<const Epetra_Map> interface_map_thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;

    //! ScaTra subproblem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! Thermo subproblem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
  };

  //! evaluation of scatra-thermo off-diagonal blocks for matching nodes
  class ScatraThermoOffDiagCouplingMatchingNodes : public ScatraThermoOffDiagCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoOffDiagCouplingMatchingNodes(
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        const Teuchos::RCP<const Epetra_Map> interface_map_scatra,
        const Teuchos::RCP<const Epetra_Map> interface_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    ~ScatraThermoOffDiagCouplingMatchingNodes() override = default;

    //! evaluation of interface contributions to scatra-thermo off-diagonal block
    void EvaluateOffDiagBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) override;

    //! evaluation of interface contributions to thermo-scatra off-diagonal block
    void EvaluateOffDiagBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) override;

   private:
    //! copy slave entries to master side scaled by -1.0
    void CopySlaveToMasterScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator>& mastermatrix);

    //! evaluate condition on slave side
    void EvaluateScatraThermoInterfaceSlaveSide(Teuchos::RCP<LINALG::SparseOperator> slavematrix);
  };

  //! evaluation of scatra-thermo off-diagonal blocks for standard Mortar on scatra discretization
  //! and condensed Bubnov Mortar on thermo discretization
  class ScatraThermoOffDiagCouplingMortarStandard : public ScatraThermoOffDiagCoupling
  {
   public:
    //! constructor
    explicit ScatraThermoOffDiagCouplingMortarStandard(
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface,
        const Teuchos::RCP<const Epetra_Map> full_map_scatra,
        const Teuchos::RCP<const Epetra_Map> full_map_thermo,
        const Teuchos::RCP<const Epetra_Map> interface_map_scatra,
        const Teuchos::RCP<const Epetra_Map> interface_map_thermo,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
        Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo);

    //! destructor
    ~ScatraThermoOffDiagCouplingMortarStandard() override = default;

    //! evaluation of interface contributions to scatra-thermo off-diagonal block
    void EvaluateOffDiagBlockScatraThermoInterface(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface) override;

    //! evaluation of interface contributions to thermo-scatra off-diagonal block
    void EvaluateOffDiagBlockThermoScatraInterface(
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface) override;
  };

}  // namespace STI
#endif  // #ifndef STI_MONOLITHIC_EVALUATE_OffDiag_H
