/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#ifndef SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H
#define SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "ssti_utils.H"

#include "../drt_adapter/adapter_coupling.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
  enum class MatrixType;
}  // namespace LINALG


namespace SSTI
{
  class SSTIMono;

  /*!
  We have three options how the global system matrix and the sub matrices are arranged:
  1) System matrix: sparse
    ->Scatra + Thermo matrix sparse
    ->Structure matrix sparse
  2) System matrix: block
    2a) Scatra + Thermo matrix block
    ->Structure matrix sparse
    2b) Scatra + Thermo matrix sparse
    ->Structure matrix sparse

  The inheritance hierarchy is appropriate*/
  class AssembleStrategyBase
  {
   public:
    AssembleStrategyBase(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    //! write 1.0 on main diagonal of slave side dofs
    virtual void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    virtual void ApplyStructuralDBCSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector> RHS, Teuchos::RCP<const Epetra_Vector> RHSscatra,
        Teuchos::RCP<const Epetra_Vector> RHSstructure,
        Teuchos::RCP<const Epetra_Vector> RHSthermo);

    //! assemble ScaTra-Block into system matrix
    virtual void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    virtual void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface) = 0;

    //! assemble ScaTra-Thermo-Block (domain contributions) into system matrix
    virtual void AssembleScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain) = 0;

    virtual void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface) = 0;

    //! assemble Structure-Block into system matrix
    virtual void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    virtual void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain) = 0;

    //! assemble Thermo-ScaTra-Block into system matrix
    virtual void AssembleThermoScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface) = 0;

    //! assemble Thermo-Structure-Block into system matrix
    virtual void AssembleThermoStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleStructureThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain) = 0;

   protected:
    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure);

    //! assemble x-structure block into system matrix for meshtying
    void AssembleXXXStructureMeshtying(LINALG::SparseMatrix& systemmatrix_x_structure,
        const LINALG::SparseMatrix& x_structurematrix);

    //! assemble structure block  into system matrix for meshtying
    void AssembleStructureMeshtying(LINALG::SparseMatrix& systemmatrix_structure,
        Teuchos::RCP<const LINALG::SparseMatrix> structurematrix);

    //! assemble structure-x block into system matrix for meshtying
    void AssembleStructureXXXMeshtying(LINALG::SparseMatrix& systemmatrix_structure_x,
        const LINALG::SparseMatrix& structures_x_matrix);

    //! Meshtying adapters
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingThermo() const
    {
      return ssti_mono_->MeshtyingThermo();
    }
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingScatra() const
    {
      return ssti_mono_->MeshtyingScatra();
    }
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> SSTIStructureMeshtying() const
    {
      return ssti_mono_->SSTIStructureMeshTying();
    }
    //@}

    //! SSTI mono maps
    const Teuchos::RCP<SSTI::SSTIMapsMono> AllMaps() const { return ssti_mono_->AllMaps(); }

    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper> StructureField() const
    {
      return ssti_mono_->StructureField();
    }

    //! interior and master side dofs of structure meshtying
    Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const;

    //! Slave side dofs of structure meshtying
    Teuchos::RCP<const Epetra_Map> MapStructureSlave() const;
    Teuchos::RCP<const Epetra_Map> MapStructureSlave3DomainIntersection() const;
    //@}

    //! Slave structure converter for structure coupling
    //@{
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter() const;
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter3DomainIntersection() const;
    //@}

    //! flags indicating meshtying
    //@{
    bool InterfaceMeshtying() const { return ssti_mono_->InterfaceMeshtying(); }
    bool Meshtying3DomainIntersection() const;
    //@}

   private:
    //! monolithic algorithm for scalar-structure-thermo interaction
    const Teuchos::RCP<const SSTI::SSTIMono> ssti_mono_;
  };

  //======================================================================================================
  //! SSTI problem is organized in sub matrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    AssembleStrategyBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override = 0;

    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override = 0;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface) override = 0;

    void AssembleScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain) override = 0;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface) override = 0;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override = 0;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override = 0;

    void AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain) override = 0;

    void AssembleThermoScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface) override = 0;

    void AssembleThermoStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface) override = 0;

    void AssembleStructureThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain) override = 0;

   protected:
    //! position of scatra blocks in system matrix
    Teuchos::RCP<std::vector<int>> BlockPositionScaTra() const { return block_position_scatra_; };

    //! position of thermo blocks in system matrix
    Teuchos::RCP<std::vector<int>> BlockPositionThermo() const { return block_position_thermo_; };

    //! position of structure block in system matrix
    int PositionStructure() const { return position_structure_; };

   private:
    //! position of scatra blocks in system matrix
    Teuchos::RCP<std::vector<int>> block_position_scatra_;

    //! position of thermo blocks in system matrix
    Teuchos::RCP<std::vector<int>> block_position_thermo_;

    //! position of structure block in system matrix
    int position_structure_;
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse structure sub matrix and block scatra sub matrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface);
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse sub matrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockSparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface);
  };

  //======================================================================================================
  //! SSTI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    AssembleStrategySparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructure(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface);
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSTI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, LINALG::MatrixType matrixtype_ssti,
      LINALG::MatrixType matrixtype_scatra);


}  // namespace SSTI
#endif  // #ifndef SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H
