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

namespace LINALG
{
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
}  // namespace LINALG

namespace SSI
{
  class SSI_Mono;
}

namespace SSI
{
  /*!
  We have three options how the global system matrix and the submatrices are arranged:
  1) System matrix: sparse
    ->Scatra matrix sparse
    ->Structure matrix sparse
  2) System matrix: block
    2a) Scatra matrix block
    ->Structure matrix sparse
    2b) Scatra matrix sparse
    ->Strucutre matrix sparse

  The inheritance hierarchy is appropriate*/
  class AssembleStrategyBase
  {
   public:
    //! constructor
    explicit AssembleStrategyBase(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono);

    //! destructor
    virtual ~AssembleStrategyBase() = default;

    //! write 1.0 on main diagonal of slave side dofs
    virtual void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    virtual void ApplyStructuralDBCSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator>& systemmatrix) = 0;

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector>& RHS, Teuchos::RCP<Epetra_Vector> RHSscatra,
        Teuchos::RCP<const Epetra_Vector> RHSstructure);

    //! assemble ScaTra-Block into system matrix
    virtual void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    virtual void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) = 0;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    virtual void AssembleScatraStructureInterface(
        Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside) = 0;

    //! assemble Structure-Block into system matrix
    virtual void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    virtual void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) = 0;

   protected:
    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure);

    //! assemble scatra-structure block (domain contributions) into system matrix for meshtying
    void AssembleScatraStructureDomainMeshtying(LINALG::SparseMatrix& systemmatrix_scatra_structure,
        const LINALG::SparseMatrix& scatrastructuredomain, bool zero);

    //! assemble structure block (domain contributions) into system matrix for meshtying
    void AssembleStructureDomainMeshtying(LINALG::SparseMatrix& systemmatrix_structure,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain, bool zero);

    //! assemble structure-scatra block (domain contributions) into system matrix for meshtying
    void AssembleStructureScatraDomainMeshtying(LINALG::SparseMatrix& systemmatrix_structure_scatra,
        const LINALG::SparseMatrix& structurescatradomain, bool zero);

    //! monolithic algorithm for scalar-structure interaction
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono_;
  };

  //======================================================================================================
  // SSI problem is organized in submatrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    //! constructor
    explicit AssembleStrategyBlock(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono);

    ~AssembleStrategyBlock() override = default;

    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    void ApplyStructuralDBCSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator>& systemmatrix) override;

    //! assemble ScaTraBlock into system matrix
    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override = 0;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside) override = 0;

    //! assemble StructureBlock into system matrix
    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override = 0;

   protected:
    //! cast system matrix to block matrix
    void CastSystemMatrixBlock(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>& systemmatrix_block);
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse structure submatrix and block scatra submatrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    //! constructor
    explicit AssembleStrategyBlockBlock(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono);

    ~AssembleStrategyBlockBlock() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse submatrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    //! constructor
    explicit AssembleStrategyBlockSparse(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono);

    ~AssembleStrategyBlockSparse() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;
  };

  //======================================================================================================
  // SSI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    //! constructor
    explicit AssembleStrategySparse(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono);

    ~AssembleStrategySparse() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void ApplyStructuralDBCSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator>& systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

   protected:
    //! cast system matrix to sparse matrix
    void CastSystemMatrixSparse(Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix>& systemmatrix_sparse);
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, LINALG::MatrixType matrixtype_ssi,
      LINALG::MatrixType matrixtype_scatra);

}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H
