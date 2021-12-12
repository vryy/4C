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
  class SSIStructureWrapper;
  class Coupling;
}  // namespace ADAPTER

namespace LINALG
{
  class MultiMapExtractor;
  class SparseOperator;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace SSI
{
  namespace UTILS
  {
    class SSIStructureMeshTying;
  }

  class ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraStructureOffDiagCoupling(
        Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure);

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void EvaluateOffDiagBlockScatraStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureblock);

    //! evaluation contributions to off-diagonal manifold scatra-structure block of global system
    //! matrix
    virtual void EvaluateOffDiagBlockScatraManifoldStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructureblock);

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void EvaluateOffDiagBlockScatraStructureInterface(
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface);

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    virtual void EvaluateOffDiagBlockStructureScatraDomain(
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) const;

   protected:
    //! copy slave side symmetric contributions to the scatra structure interface linearization
    //! entries to master side scaled by -1.0
    void CopySlaveToMasterScatraStructureSymmetricInterfaceContributions(
        Teuchos::RCP<const LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator>& mastermatrix);

    //! evaluate symmetric contributions to the scatra structure interface linearization on the
    //! slave side
    void EvaluateScatraStructureSymmetricInterfaceContributionsSlaveSide(
        Teuchos::RCP<LINALG::SparseOperator> slavematrix);

    //! evaluate non-symmetric contributions to the scatra structure interface linearization
    void EvaluateScatraStructureNonSymmetricInterfaceContributionsSlaveSide(
        Teuchos::RCP<LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<LINALG::SparseOperator> mastermatrix);

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> FullMapStructure() const { return full_map_structure_; }

    //! scatra discretization
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const { return scatra_; }

   private:
    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! scatra discretization
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    //! structure problem
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying_;
  };

  class ScatraManifoldStructureOffDiagCoupling : public ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraManifoldStructureOffDiagCoupling(
        Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure);

    void EvaluateOffDiagBlockScatraManifoldStructureDomain(
        Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructureblock) override;

   private:
    //! scatra manifold discretization
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold_;
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
        Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_scatra,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure);

    void EvaluateOffDiagBlockStructureScatraDomain(
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) const override;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Epetra_Map> full_map_scatra_;
  };
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_EVALUATE_OffDiag_H
