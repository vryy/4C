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
  We have three options how the global system matrix and the submatrices are arranged:
  1) System matrix: sparse
    ->Scatra + Thermo matrix sparse
    ->Structure matrix sparse
  2) System matrix: block
    2a) Scatra + Thermo matrix block
    ->Structure matrix sparse
    2b) Scatra + Thermo matrix sparse
    ->Strucutre matrix sparse

  The inheritance hierarchy is appropriate*/
  class AssembleStrategyBase
  {
   public:
    AssembleStrategyBase(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    virtual ~AssembleStrategyBase() = default;

    //! write 1.0 on main diagonal of slave side dofs
    virtual void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    virtual void ApplyStructuralDBCSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix) = 0;

    //! assemble RHS
    void AssembleRHS(Teuchos::RCP<Epetra_Vector> RHS, Teuchos::RCP<Epetra_Vector> RHSscatra,
        Teuchos::RCP<const Epetra_Vector> RHSstructure,
        Teuchos::RCP<const Epetra_Vector> RHSthermo);

    //! assemble ScaTra-Block into system matrix
    virtual void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    virtual void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) = 0;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    virtual void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside) = 0;

    virtual void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) = 0;

    //! assemble Structure-Block into system matrix
    virtual void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    virtual void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) = 0;

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

    void AssembleStructureThermoDomainMeshtying(LINALG::SparseMatrix& systemmatrix_structure_thermo,
        const LINALG::SparseMatrix& structurethermodomain, bool zero);

    void AssembleThermoStructureDomainMeshtying(LINALG::SparseMatrix& systemmatrix_thermo_structure,
        const LINALG::SparseMatrix& thermostructuredomain, bool zero);

    const Teuchos::RCP<SSTI::SSTIStructuralMeshtying> StructuralMeshtying() const
    {
      return ssti_mono_->StructuralMeshtying();
    }

    bool InterfaceMeshtying() const { return ssti_mono_->InterfaceMeshtying(); }

    bool Meshtying3DomainIntersection() const
    {
      return StructuralMeshtying()->Meshtying3DomainIntersection();
    }

    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingThermo() const
    {
      return ssti_mono_->MeshtyingThermo();
    }

    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingScatra() const
    {
      return ssti_mono_->MeshtyingScatra();
    }

    const Teuchos::RCP<SSTI::SSTIMapsMono> AllMaps() const { return ssti_mono_->AllMaps(); }

    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper> StructureField() const
    {
      return ssti_mono_->StructureField();
    }

    Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const
    {
      return StructuralMeshtying()->MapStructureCondensed();
    }
    Teuchos::RCP<const Epetra_Map> MapStructureSlave() const
    {
      return StructuralMeshtying()->MapsInterfaceStructure()->Map(0);
    }

   private:
    //! monolithic algorithm for scalar-structure-thermo interaction
    const Teuchos::RCP<const SSTI::SSTIMono> ssti_mono_;
  };

  //======================================================================================================
  //! SSTI problem is organized in submatrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    AssembleStrategyBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    ~AssembleStrategyBlock() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override = 0;

    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override = 0;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override = 0;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside) override = 0;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override = 0;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override = 0;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override = 0;

    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override = 0;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override = 0;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override = 0;

    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override = 0;

    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override = 0;

    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override = 0;

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
  //! SSTI problem is organized in sparse structure submatrix and block scatra submatrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    ~AssembleStrategyBlockBlock() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse submatrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockSparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    ~AssembleStrategyBlockSparse() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;
  };

  //======================================================================================================
  //! SSTI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    AssembleStrategySparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    ~AssembleStrategySparse() override = default;

    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;
  };

  //! build specifc assemble strategy
  Teuchos::RCP<SSTI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, LINALG::MatrixType matrixtype_ssti,
      LINALG::MatrixType matrixtype_scatra);


}  // namespace SSTI
#endif  // #ifndef SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H
