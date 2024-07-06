/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_HPP
#define FOUR_C_SSTI_MONOLITHIC_ASSEMBLE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_ssti_monolithic.hpp"
#include "4C_ssti_utils.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
  enum class MatrixType;
}  // namespace Core::LinAlg


namespace SSTI
{

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
    /**
     * Virtual destructor.
     */
    virtual ~AssembleStrategyBase() = default;

    AssembleStrategyBase(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    //! write 1.0 on main diagonal of slave side dofs
    virtual void apply_meshtying_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) = 0;

    //! apply structural Dirichlet boundary conditions on system matrix
    virtual void apply_structural_dbc_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) = 0;

    //! assemble RHS
    void assemble_rhs(Teuchos::RCP<Epetra_Vector> RHS, Teuchos::RCP<const Epetra_Vector> RHSscatra,
        Teuchos::RCP<const Epetra_Vector> RHSstructure,
        Teuchos::RCP<const Epetra_Vector> RHSthermo);

    //! assemble ScaTra-Block into system matrix
    virtual void assemble_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain) = 0;

    //! assemble ScaTra-Structure-Block (domain contributions) into system matrix
    virtual void assemble_scatra_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface) = 0;

    //! assemble ScaTra-Thermo-Block (domain contributions) into system matrix
    virtual void assemble_scatra_thermo_domain(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain) = 0;

    virtual void assemble_scatra_thermo_interface(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface) = 0;

    //! assemble Structure-Block into system matrix
    virtual void assemble_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain) = 0;

    //! assemble Structure-ScaTra-Block (domain contributions) into system matrix
    virtual void assemble_structure_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void assemble_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain) = 0;

    //! assemble Thermo-ScaTra-Block into system matrix
    virtual void assemble_thermo_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface) = 0;

    //! assemble Thermo-Structure-Block into system matrix
    virtual void assemble_thermo_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface) = 0;

    //! assemble Thermo-Block into system matrix
    virtual void assemble_structure_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain) = 0;

   protected:
    //! write 1.0 on main diagonal of slave side dofs
    void apply_meshtying_sys_mat(Core::LinAlg::SparseMatrix& systemmatrix_structure);

    //! assemble x-structure block into system matrix for meshtying
    void assemble_xxx_structure_meshtying(Core::LinAlg::SparseMatrix& systemmatrix_x_structure,
        const Core::LinAlg::SparseMatrix& x_structurematrix);

    //! assemble structure block  into system matrix for meshtying
    void assemble_structure_meshtying(Core::LinAlg::SparseMatrix& systemmatrix_structure,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain);

    //! assemble structure-x block into system matrix for meshtying
    void assemble_structure_xxx_meshtying(Core::LinAlg::SparseMatrix& systemmatrix_structure_x,
        const Core::LinAlg::SparseMatrix& structures_x_matrix);

    //! Meshtying adapters
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_thermo() const
    {
      return ssti_mono_->meshtying_thermo();
    }
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_scatra() const
    {
      return ssti_mono_->meshtying_scatra();
    }
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying() const
    {
      return ssti_mono_->ssti_structure_mesh_tying();
    }
    //@}

    //! SSTI mono maps
    Teuchos::RCP<SSTI::SSTIMapsMono> all_maps() const { return ssti_mono_->all_maps(); }

    Teuchos::RCP<Adapter::SSIStructureWrapper> structure_field() const
    {
      return ssti_mono_->structure_field();
    }

    //! flag indicating meshtying
    bool interface_meshtying() const { return ssti_mono_->interface_meshtying(); }

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

    void apply_meshtying_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override = 0;

    void apply_structural_dbc_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override;

    void assemble_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain) override = 0;

    void assemble_scatra_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface) override = 0;

    void assemble_scatra_thermo_domain(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain) override = 0;

    void assemble_scatra_thermo_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface) override = 0;

    void assemble_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain) override = 0;

    void assemble_structure_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain) override = 0;

    void assemble_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain) override = 0;

    void assemble_thermo_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface) override = 0;

    void assemble_thermo_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface) override = 0;

    void assemble_structure_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain) override = 0;

   protected:
    //! position of scatra blocks in system matrix
    std::vector<int> block_position_sca_tra() const { return block_position_scatra_; };

    //! position of thermo blocks in system matrix
    std::vector<int> block_position_thermo() const { return block_position_thermo_; };

    //! position of structure block in system matrix
    int position_structure() const { return position_structure_; };

   private:
    //! position of scatra blocks in system matrix
    std::vector<int> block_position_scatra_;

    //! position of thermo blocks in system matrix
    std::vector<int> block_position_thermo_;

    //! position of structure block in system matrix
    int position_structure_;
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse structure sub matrix and block scatra sub matrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void apply_meshtying_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override;

    void assemble_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain) override;

    void assemble_scatra_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface) override;

    void assemble_scatra_thermo_domain(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain) override;

    void assemble_scatra_thermo_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface) override;

    void assemble_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain) override;

    void assemble_structure_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain) override;

    void assemble_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain) override;

    void assemble_thermo_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface) override;

    void assemble_thermo_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface) override;

    void assemble_structure_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void assemble_thermo_scatra_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface);
  };

  // *********************************************************************************************
  //! SSTI problem is organized in sparse sub matrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    AssembleStrategyBlockSparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void apply_meshtying_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override;

    void assemble_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain) override;

    void assemble_scatra_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface) override;

    void assemble_scatra_thermo_domain(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain) override;

    void assemble_scatra_thermo_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface) override;

    void assemble_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain) override;

    void assemble_structure_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain) override;

    void assemble_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain) override;

    void assemble_thermo_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface) override;

    void assemble_thermo_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface) override;

    void assemble_structure_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void assemble_thermo_scatra_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface);
  };

  //======================================================================================================
  //! SSTI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    AssembleStrategySparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono);

    void apply_meshtying_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override;

    void apply_structural_dbc_system_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix) override;

    void assemble_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain) override;

    void assemble_scatra_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface) override;

    void assemble_scatra_thermo_domain(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain) override;

    void assemble_scatra_thermo_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface) override;

    void assemble_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain) override;

    void assemble_structure_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain) override;

    void assemble_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain) override;

    void assemble_thermo_scatra(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface) override;

    void assemble_thermo_structure(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface) override;

    void assemble_structure_thermo(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain) override;

   private:
    //! assemble interface contribution from thermo-scatra block
    void assemble_thermo_scatra_interface(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
        Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface);
  };

  //! build specific assemble strategy
  Teuchos::RCP<SSTI::AssembleStrategyBase> BuildAssembleStrategy(
      Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, Core::LinAlg::MatrixType matrixtype_ssti,
      Core::LinAlg::MatrixType matrixtype_scatra);

}  // namespace SSTI
FOUR_C_NAMESPACE_CLOSE

#endif
