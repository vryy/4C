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
  We have three options how the global system matrix and the submatrices are arranged:
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
    explicit AssembleStrategyBase(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

    //! write 1.0 on main diagonal of slave side dofs
    virtual void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<const Epetra_Vector> rhs_scatra,
        Teuchos::RCP<const Epetra_Vector> rhs_structure,
        Teuchos::RCP<const Epetra_Vector> rhs_manifold);

    //! assemble ScaTra-Block into system matrix
    virtual void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block  into system matrix
    virtual void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) = 0;

    //! assemble Structure-Block into system matrix
    virtual void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block into system matrix
    virtual void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) = 0;

    //! assemble ScaTra Manifold-Block into system matrix
    virtual void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifolddomain) = 0;

    //! assemble ScaTra Manifold-Structure-Block into system matrix
    virtual void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain) = 0;

   protected:
    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure);

    //! assemble xxx-structure block into system matrix for meshtying
    void AssembleXXXStructureMeshtying(LINALG::SparseMatrix& systemmatrix_xxx_structure,
        const LINALG::SparseMatrix& xxx_structurematrix);

    //! assemble structure block  into system matrix for meshtying
    void AssembleStructureMeshtying(LINALG::SparseMatrix& systemmatrix_structure,
        Teuchos::RCP<LINALG::SparseMatrix> structurematrix);

    //! assemble structure-xxx block into system matrix for meshtying
    void AssembleStructureXXXMeshtying(LINALG::SparseMatrix& systemmatrix_structure_xxx,
        const LINALG::SparseMatrix& structures_xxx_matrix);

    //! monolithic algorithm for scalar-structure interaction
    const Teuchos::RCP<const SSI::SSIMono> ssi_mono_;

    //! interface maps (condensed: inner and master side, slave side, 3 domain intersection slave
    //! side)
    //@{
    Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const { return mapstructurecondensed_; }
    Teuchos::RCP<const Epetra_Map> MapStructureSlave() const { return mapstructureslave_; }
    Teuchos::RCP<const Epetra_Map> MapStructureSlave3DomainIntersection() const
    {
      return mapstructureslave3domainintersection_;
    }
    //@}

    // Slave structure converter for structure coupling
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter() const
    {
      return slave_side_converter_->InterfaceCouplingAdapterStructureSlaveConverter();
    }

    //! Slave structure converter for structure coupling with 3 domain intersection
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter3DomainIntersection() const
    {
      return slave_side_converter_
          ->InterfaceCouplingAdapterStructureSlaveConverter3DomainIntersection();
    }

    //! indicates, if meshtying includes three entities
    bool Meshtying3DomainIntersection() const { return meshtying_3_domain_intersection_; }

   private:
    //! interface maps (condensed: inner and master side, slave side, 3 domain intersection slave
    //! side)
    //@{
    Teuchos::RCP<const Epetra_Map> mapstructurecondensed_;
    Teuchos::RCP<const Epetra_Map> mapstructureslave_;
    Teuchos::RCP<const Epetra_Map> mapstructureslave3domainintersection_;
    //@}

    //! slave side converter
    Teuchos::RCP<SSI::UTILS::SSISlaveSideConverter> slave_side_converter_;

    //! indicates, if meshtying includes three entities
    bool meshtying_3_domain_intersection_;
  };

  //======================================================================================================
  // SSI problem is organized in submatrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategyBlock(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override = 0;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override = 0;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override = 0;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override = 0;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override = 0;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifolddomain) override = 0;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain) override = 0;

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
    explicit AssembleStrategyBlockBlock(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain) override;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse submatrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockSparse(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain) override;
  };

  //======================================================================================================
  // SSI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategySparse(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleScaTraManifold(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifolddomain) override;

    void AssembleScaTraManifoldStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain) override;
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSI::SSIMono> ssi_mono, LINALG::MatrixType matrixtype_ssi,
      LINALG::MatrixType matrixtype_scatra);

}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_ASSEMBLE_STRATEGY_H
