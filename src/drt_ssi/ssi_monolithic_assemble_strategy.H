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
    explicit AssembleStrategyBase(const SSI::SSIMono& ssi_mono);

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_scatra,
        Teuchos::RCP<const Epetra_Vector> rhs_structure,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_manifold_side,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_scatra_side);

    //! assemble ScaTra-Block into system matrix
    virtual void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block into system matrix
    virtual void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) = 0;

    //! assemble Structure-Block into system matrix
    virtual void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block into system matrix
    virtual void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) = 0;

    //! assemble ScaTra Manifold-Block into system matrix
    virtual void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifolddomain) = 0;

    //! assemble ScaTra Manifold-Structure-Block into system matrix
    virtual void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldstructuredomain) = 0;

    //! assemble ScaTra Manifold-Scatra-Block into system matrix
    virtual void AssembleScaTraManifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldscatradomain) = 0;

    //! assemble ScaTra-Manifold Scatra-Block into system matrix
    virtual void AssembleScatraScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifolddomain) = 0;

   protected:
    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& SSIMono() const { return ssi_mono_; }

   private:
    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& ssi_mono_;
  };

  //======================================================================================================
  // SSI problem is organized in submatrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategyBlock(const SSI::SSIMono& ssi_mono);

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override = 0;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override = 0;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override = 0;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override = 0;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifolddomain) override = 0;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldstructuredomain) override = 0;

    void AssembleScaTraManifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldscatradomain) override = 0;

    void AssembleScatraScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifolddomain) override = 0;

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
  // SSI problem is organized in sparse structure submatrix and block scatra submatrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockBlock(const SSI::SSIMono& ssi_mono);

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldstructuredomain) override;

    void AssembleScaTraManifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldscatradomain) override;

    void AssembleScatraScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifolddomain) override;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse submatrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockSparse(const SSI::SSIMono& ssi_mono);

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldstructuredomain) override;

    void AssembleScaTraManifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldscatradomain) override;

    void AssembleScatraScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifolddomain) override;
  };

  //======================================================================================================
  // SSI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategySparse(const SSI::SSIMono& ssi_mono);

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldstructuredomain) override;

    void AssembleScaTraManifoldScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> manifoldscatradomain) override;

    void AssembleScatraScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatramanifolddomain) override;
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSI::AssembleStrategyBase> BuildAssembleStrategy(const SSI::SSIMono& ssi_mono,
      LINALG::MatrixType matrixtype_ssi, LINALG::MatrixType matrixtype_scatra);

}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H
