/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#ifndef SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H
#define SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "../drt_adapter/adapter_coupling.H"

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
    //! constructor
    AssembleStrategyBase(
        Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter);

    //! destructor
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

    //! structural coupling converter
    ADAPTER::CouplingSlaveConverter converter_;

    //! monolithic algorithm for scalar-structure interaction
    const Teuchos::RCP<const SSTI::SSTIMono> ssti_mono_;
  };

  //======================================================================================================
  //! SSTI problem is organized in submatrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    //! constructor
    AssembleStrategyBlock(
        Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter);

    //! destructor
    ~AssembleStrategyBlock() override = default;

    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    //! assemble ScaTraBlock into system matrix
    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override = 0;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside) override = 0;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override = 0;

    //! assemble StructureBlock into system matrix
    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override = 0;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override = 0;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override = 0;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override = 0;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override = 0;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override = 0;

    //! assemble Thermo-Block into system matrix
    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override = 0;

   protected:
    //! cast system matrix to block matrix
    void CastSystemMatrixBlock(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>& systemmatrix_block);
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse structure submatrix and block scatra submatrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    //! constructor
    AssembleStrategyBlockBlock(
        Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter);

    //! destructor
    ~AssembleStrategyBlockBlock() override = default;

    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    //! assemble ScaTraBlock into system matrix
    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    //! assemble StructureBlock into system matrix
    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse submatrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    //! constructor
    AssembleStrategyBlockSparse(
        Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter);

    //! destructor
    ~AssembleStrategyBlockSparse() override = default;

    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    //! assemble ScaTraBlock into system matrix
    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    //! assemble StructureBlock into system matrix
    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;
  };

  //======================================================================================================
  //! SSTI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    //! constructor
    AssembleStrategySparse(
        Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter);

    //! destructor
    ~AssembleStrategySparse() override = default;

    //! write 1.0 on main diagonal of slave side dofs
    void ApplyMeshtyingSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    //! apply structural Dirichlet boundary conditions on system matrix
    void ApplyStructuralDBCSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> systemmatrix) override;

    //! assemble ScaTraBlock into system matrix
    void AssembleScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatradomain) override;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    void AssembleScatraStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain) override;

    //! assemble ScaTra-Structure-Block (interface contributions) into system matrix
    void AssembleScatraStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface) override;

    void AssembleScatraThermoInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface) override;

    //! assemble StructureBlock into system matrix
    void AssembleStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix> structuredomain) override;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    void AssembleStructureScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermodomain) override;

    void AssembleThermoScatraDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain) override;

    void AssembleThermoScatraInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain) override;

    //! assemble Thermo-Block into system matrix
    void AssembleThermoStructureInterface(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface) override;

    //! assemble Thermo-Block into system matrix
    void AssembleStructureThermoDomain(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseOperator> structurethermodomain) override;

   protected:
    //! cast system matrix to sparse matrix
    void CastSystemMatrixSparse(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<LINALG::SparseMatrix>& systemmatrix_sparse);
  };

  //! build specifc assemble strategy
  Teuchos::RCP<SSTI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSTI::SSTIMono> ssti_mono,
      const ADAPTER::CouplingSlaveConverter& converter, LINALG::MatrixType matrixtype_ssti,
      LINALG::MatrixType matrixtype_scatra);


}  // namespace SSTI
#endif  // #ifndef SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_H
