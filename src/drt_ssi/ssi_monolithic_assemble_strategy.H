/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H
#define SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "../linalg/linalg_sparseoperator.H"
#include "../drt_adapter/adapter_coupling.H"
#include "ssi_utils.H"

namespace LINALG
{
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
}  // namespace LINALG

namespace SSI
{
  namespace UTILS
  {
    class SSISlaveSideConverter;
  }
  class SSIMono;

  /*!
  We have three options how the global system matrix and the sub matrices are arranged:
  1) System matrix: sparse
    ->Scatra matrix sparse
    ->Structure matrix sparse
  2) System matrix: block
    2a) Scatra matrix block
    ->Structure matrix sparse
    2b) Scatra matrix sparse
    ->Structure matrix sparse

  The inheritance hierarchy is appropriate*/
  class AssembleStrategyBase
  {
   public:
    //! constructor
    explicit AssembleStrategyBase(
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_scatra,
        Teuchos::RCP<const Epetra_Vector> rhs_structure,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_manifold_side,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_scatra_side);

    //! assemble ScaTra-ScaTra-Block into system matrix
    virtual void AssembleScatraScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix) = 0;

    //! assemble ScaTra-Structure-Block into system matrix
    virtual void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) = 0;

    //! assemble Structure-Structure-Block into system matrix
    virtual void AssembleStructureStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix) = 0;

    //! assemble Structure-ScaTra-Block into system matrix
    virtual void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix) = 0;

    //! assemble ScaTra Manifold-ScaTra Manifold-Block into system matrix
    virtual void AssembleScatramanifoldScatramanifold(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix) = 0;

    //! assemble ScaTra Manifold-Structure-Block into system matrix
    virtual void AssembleScatramanifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix) = 0;

    //! assemble ScaTra Manifold-ScaTra-Block into system matrix
    virtual void AssembleScatramanifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix) = 0;

    //! assemble ScaTra-ScaTra Manifold-Block into system matrix
    virtual void AssembleScatraScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix) = 0;

   protected:
    //! solve additional scatra field on manifolds
    bool IsScaTraManifold() const { return is_scatra_manifold_; }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

   private:
    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps_;
  };

  //======================================================================================================
  // SSI problem is organized in sub matrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategyBlock(
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void AssembleScatraScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix) override = 0;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override = 0;

    void AssembleStructureStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix) override = 0;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix) override = 0;

    void AssembleScatramanifoldScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix) override =
        0;

    void AssembleScatramanifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix) override = 0;

    void AssembleScatramanifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix) override = 0;

    void AssembleScatraScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix) override = 0;

   protected:
    //! position of scatra blocks in system matrix
    Teuchos::RCP<std::vector<int>> BlockPositionScaTra() const { return block_position_scatra_; }

    //! position of scatra manifold blocks in system matrix
    Teuchos::RCP<std::vector<int>> BlockPositionScaTraManifold() const
    {
      return block_position_scatra_manifold_;
    }

    //! position of structure block in system matrix
    int PositionStructure() const { return position_structure_; };

   private:
    //! position of scatra blocks in system matrix
    Teuchos::RCP<std::vector<int>> block_position_scatra_;

    //! position of scatra manifold blocks in system matrix
    Teuchos::RCP<std::vector<int>> block_position_scatra_manifold_;

    //! position of structure block in system matrix
    int position_structure_;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse structure sub matrix and block scatra sub matrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockBlock(
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void AssembleScatraScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructureStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix) override;

    void AssembleScatramanifoldScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix) override;

    void AssembleScatramanifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix) override;

    void AssembleScatramanifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix) override;

    void AssembleScatraScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse sub matrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockSparse(
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void AssembleScatraScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructureStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix) override;

    void AssembleScatramanifoldScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix) override;

    void AssembleScatramanifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix) override;

    void AssembleScatramanifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix) override;

    void AssembleScatraScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  //======================================================================================================
  // SSI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategySparse(
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void AssembleScatraScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructureStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix) override;

    void AssembleScatramanifoldScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix) override;

    void AssembleScatramanifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix) override;

    void AssembleScatramanifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix) override;

    void AssembleScatraScatramanifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold,
      LINALG::MatrixType matrixtype_ssi, LINALG::MatrixType matrixtype_scatra);

}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H
