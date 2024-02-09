/*----------------------------------------------------------------------*/
/*! \file
\brief Mesh tying strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef BACI_SSI_MONOLITHIC_MESHTYING_STRATEGY_HPP
#define BACI_SSI_MONOLITHIC_MESHTYING_STRATEGY_HPP

#include "baci_config.hpp"

#include "baci_linalg_sparseoperator.hpp"

#include <Teuchos_RCPDecl.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace ADAPTER
{
  class CouplingSlaveConverter;
}

namespace SSI
{
  namespace UTILS
  {
    class SSIMaps;
    class SSIMeshTying;
  }  // namespace UTILS

  class SSIMono;

  //! base functionality for scatra structure interaction mesh tying
  class MeshtyingStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~MeshtyingStrategyBase() = default;

    //! constructor
    explicit MeshtyingStrategyBase(bool is_scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying);

    /*!
     * @brief apply mesh tying to structure matrix
     *
     * @param[out] ssi_structure_matrix  structure matrix including mesh tying constraints
     * @param[in] structure_matrix       structure matrix from structure problem
     * @param[in] do_uncomplete          flag indicating if we need to uncomplete the matrix before
     *                                   adding something
     */
    void ApplyMeshtyingToStructureMatrix(CORE::LINALG::SparseMatrix& ssi_structure_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseMatrix> structure_matrix, bool do_uncomplete);

    /*!
     * @brief apply mesh tying to scatra manifold structure matrix
     *
     * @param[in,out] manifold_structure_matrix  scalar transport on manifold structure matrix
     * @param[in]     do_uncomplete              flag indicating if we need to uncomplete the matrix
     *                                           before adding something
     */
    virtual void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) = 0;

    /*!
     * @brief apply mesh tying to scatra structure matrix
     *
     * @param[in,out] scatra_structure_matrix  scatra structure matrix
     * @param[in]     do_uncomplete            flag indicating if we need to uncomplete the matrix
     *                                         before adding something
     */
    virtual void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix, bool do_uncomplete) = 0;

    /*!
     * @brief apply mesh tying to structure right hand side vector
     *
     * @param[in] structure_rhs  structure right hand side vector without mesh tying contributions
     * @return structure right hand side vector including mesh tying contributions
     */
    Epetra_Vector ApplyMeshtyingToStructureRHS(Teuchos::RCP<const Epetra_Vector> structure_rhs);

    /*!
     * @brief apply mesh tying to the structure scatra matrix
     *
     * @param[in,out] structure_scatra_matrix  structure scatra matrix
     * @param[in]     do_uncomplete            flag indicating if we need to uncomplete the matrix
     *                                         before adding something
     */
    virtual void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix, bool do_uncomplete) = 0;

   protected:
    /*!
     * @brief apply mesh tying to structure-xxx block
     *
     * @param[out] ssi_structure_xxx_matrix  structure xxx matrix block including mesh tying
     *                                       constraints
     * @param[in] structure_xxx_matrix       structure xxx matrix block
     */
    void ApplyMeshtyingToStructureXXX(CORE::LINALG::SparseMatrix& ssi_structure_xxx_matrix,
        const CORE::LINALG::SparseMatrix& structure_xxx_matrix);

    /*!
     * @brief apply mesh tying to xxx-structure block
     *
     * @param[out] ssi_xxx_structure_matrix  xxx structure matrix block including mesh tying
     *                                       constraints
     * @param[in] xxx_structure_matrix       xxx structure matrix block
     */
    void ApplyMeshtyingToXXXStructure(CORE::LINALG::SparseMatrix& ssi_xxx_structure_matrix,
        const CORE::LINALG::SparseMatrix& xxx_structure_matrix);

    //! solve additional scatra field on manifolds
    bool IsScaTraManifold() const { return is_scatra_manifold_; }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> SSIStructureMeshtying() const
    {
      return ssi_structure_meshtying_;
    }

    //! scatra structure contribution matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> temp_scatra_struct_mat_;

    //! scatra-manifold structure system matrix used to apply mesh tying to this matrix block
    Teuchos::RCP<CORE::LINALG::SparseOperator> temp_scatramanifold_struct_mat_;

    //! structure scatra system matrix used to apply mesh tying to this matrix block
    Teuchos::RCP<CORE::LINALG::SparseOperator> temp_struct_scatra_mat_;

   private:
    /*!
     * @brief finalize mesh tying by applying pseudo dirichlet conditions to the slave side, i.e.
     * write 1.0 on main diagonal of slave side dofs
     *
     * @param[in,out] ssi_structure_matrix  structure matrix with mesh tying constraints
     */
    void FinalizeMeshtyingStructureMatrix(CORE::LINALG::SparseMatrix& ssi_structure_matrix);

    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps_;

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying_;
  };

  //! SSI problem is represented by one sparse matrix
  class MeshtyingStrategySparse : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategySparse(bool is_scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix,
        bool do_uncomplete) override;
  };

  //! SSI problem is composed of sub matrices
  class MeshtyingStrategyBlock : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBlock(bool is_scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix,
        bool do_uncomplete) override;

   protected:
    //! position of scatra blocks in system matrix
    const std::vector<int>& BlockPositionScaTra() const { return block_position_scatra_; }

    //! position of scatra manifold blocks in system matrix
    const std::vector<int>& BlockPositionScaTraManifold() const
    {
      return block_position_scatra_manifold_;
    }

    //! position of structure block in system matrix
    int PositionStructure() const { return position_structure_; };

   private:
    //! position of scatra blocks in system matrix
    const std::vector<int> block_position_scatra_;

    //! position of scatra manifold blocks in system matrix
    std::vector<int> block_position_scatra_manifold_;

    //! position of structure block in system matrix
    const int position_structure_;
  };

  //! build specific mesh tying strategy
  Teuchos::RCP<SSI::MeshtyingStrategyBase> BuildMeshtyingStrategy(bool is_scatra_manifold,
      CORE::LINALG::MatrixType matrixtype_scatra, Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
      Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying);
}  // namespace SSI
BACI_NAMESPACE_CLOSE

#endif  // SSI_MONOLITHIC_MESHTYING_STRATEGY_H
