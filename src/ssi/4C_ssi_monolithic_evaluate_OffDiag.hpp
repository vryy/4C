/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSI

\level 2


 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SSI_MONOLITHIC_EVALUATE_OFFDIAG_HPP
#define FOUR_C_SSI_MONOLITHIC_EVALUATE_OFFDIAG_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class SSIStructureWrapper;
  class Coupling;
}  // namespace ADAPTER

namespace CORE::LINALG
{
  class MultiMapExtractor;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace SSI
{
  namespace UTILS
  {
    class SSIMeshTying;
  }

  class ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraStructureOffDiagCoupling(
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

    virtual ~ScatraStructureOffDiagCoupling() = default;

    //! evaluate domain contributions to off-diagonal scatra-structure block of global system matrix
    void evaluate_off_diag_block_scatra_structure_domain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructureblock);

    //! evaluation contributions to off-diagonal manifold scatra-structure block of global system
    //! matrix
    virtual void evaluate_off_diag_block_scatra_manifold_structure_domain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatramanifoldstructureblock);

    //! evaluate interface contributions to  off-diagonal scatra-structure block of global system
    //! matrix
    void evaluate_off_diag_block_scatra_structure_interface(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructureinterface);

    //! evaluate domain contributions to off-diagonal structure-scatra block of global system matrix
    virtual void evaluate_off_diag_block_structure_scatra_domain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structurescatradomain) const;

   protected:
    //! copy slave side symmetric contributions to the scatra structure interface linearization
    //! entries to master side scaled by -1.0
    void copy_slave_to_master_scatra_structure_symmetric_interface_contributions(
        Teuchos::RCP<const CORE::LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& mastermatrix);

    //! evaluate symmetric contributions to the scatra structure interface linearization on the
    //! slave side
    void evaluate_scatra_structure_symmetric_interface_contributions_slave_side(
        Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix);

    //! evaluate non-symmetric contributions to the scatra structure interface linearization
    void evaluate_scatra_structure_non_symmetric_interface_contributions_slave_side(
        Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<CORE::LINALG::SparseOperator> mastermatrix);

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure() const { return full_map_structure_; }

    //! scatra discretization
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> sca_tra_field() const { return scatra_; }

   private:
    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! scatra discretization
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    //! structure problem
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure_;

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying_;
  };

  class ScatraManifoldStructureOffDiagCoupling : public ScatraStructureOffDiagCoupling
  {
   public:
    explicit ScatraManifoldStructureOffDiagCoupling(
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

    void evaluate_off_diag_block_scatra_manifold_structure_domain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatramanifoldstructureblock) override;

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
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
        Teuchos::RCP<const Epetra_Map> full_map_scatra,
        Teuchos::RCP<const Epetra_Map> full_map_structure,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

    void evaluate_off_diag_block_structure_scatra_domain(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structurescatradomain) const override;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Epetra_Map> full_map_scatra_;
  };
}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
