/*----------------------------------------------------------------------*/
/*! \file
\brief Mesh tying strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_MESHTYING_STRATEGY_H
#define SSI_MONOLITHIC_MESHTYING_STRATEGY_H

#include "Teuchos_RCPDecl.hpp"
#include "vector"

#include "../linalg/linalg_sparseoperator.H"

namespace LINALG
{
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace ADAPTER
{
  class CouplingSlaveConverter;
}

namespace SSI
{
  namespace UTILS
  {
    class SSIMaps;
    class SSISlaveSideConverter;
    class SSIStructureMeshTying;
  }  // namespace UTILS

  class SSIMono;

  //! base functionality for scatra structure interaction mesh tying
  class MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBase(bool is_scatra_manifold, bool meshtying_3_domain_intersection,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);

    /*!
     * @brief apply mesh tying to structure matrix
     *
     * @param[out] ssi_structure_matrix  structure matrix including mesh tying constraints
     * @param[in] structure_matrix       structure matrix from structure problem
     */
    void ApplyMeshtyingToStructureMatrix(LINALG::SparseMatrix& ssi_structure_matrix,
        Teuchos::RCP<const LINALG::SparseMatrix> structure_matrix);

    /*!
     * @brief apply mesh tying to scatra manifold structure matrix
     *
     * @param[in,out] manifold_structure_matrix  scalar transport on manifold structure matrix
     * @param[in]     do_uncomplete              flag indicating if we need to uncomplete the matrix
     *                                           before adding something
     */
    virtual void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix, bool do_uncomplete) = 0;

    /*!
     * @brief apply mesh tying to scatra structure matrix
     *
     * @param[in,out] scatra_structure_matrix  scatra structure matrix
     * @param[in]     do_uncomplete            flag indicating if we need to uncomplete the matrix
     *                                         before adding something
     */
    virtual void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix, bool do_uncomplete) = 0;

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
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix, bool do_uncomplete) = 0;

   protected:
    /*!
     * @brief apply mesh tying to structure-xxx block
     *
     * @param[out] ssi_structure_xxx_matrix  structure xxx matrix block including mesh tying
     *                                       constraints
     * @param[in] structure_xxx_matrix       structure xxx matrix block
     */
    void ApplyMeshtyingToStructureXXX(LINALG::SparseMatrix& ssi_structure_xxx_matrix,
        const LINALG::SparseMatrix& structure_xxx_matrix);

    /*!
     * @brief apply mesh tying to xxx-structure block
     *
     * @param[out] ssi_xxx_structure_matrix  xxx structure matrix block including mesh tying
     *                                       constraints
     * @param[in] xxx_structure_matrix       xxx structure matrix block
     */
    void ApplyMeshtyingToXXXStructure(LINALG::SparseMatrix& ssi_xxx_structure_matrix,
        const LINALG::SparseMatrix& xxx_structure_matrix);

    //! solve additional scatra field on manifolds
    bool IsScaTraManifold() const { return is_scatra_manifold_; }

    //! structural map extractor for coupling of two entities (0: other degrees of freedom, 1:
    //! slave-side degrees of freedom, 2: master-side degrees of freedom)
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsCoupStruct() const;

    //! Slave side dofs of structural meshtying
    Teuchos::RCP<const Epetra_Map> MapStructureSlave() const;
    Teuchos::RCP<const Epetra_Map> MapStructureSlave3DomainIntersection() const;

    //! structural map extractor for 3 domains intersecting (0: other degrees of freedom, 1:
    //! slave-side degrees  of freedom, 2: master-side degrees of freedom)
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsCoupStruct3DomainIntersection() const;

    //! return map for interior and master-side structural degrees of freedom
    Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const;

    //! return the slave structure converter for structure coupling
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter() const;

    //! return the slave structure converter for structure coupling with 3 domain intersection
    ADAPTER::CouplingSlaveConverter& StructureSlaveConverter3DomainIntersection() const;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> SSIStructureMeshtying() const
    {
      return ssi_structure_meshtying_;
    }

    //! scatra structure contribution matrix
    Teuchos::RCP<LINALG::SparseOperator> temp_scatra_struct_mat_;

    //! scatra-manifold structure system matrix used to apply mesh tying to this matrix block
    Teuchos::RCP<LINALG::SparseOperator> temp_scatramanifold_struct_mat_;

    //! structure scatra system matrix used to apply mesh tying to this matrix block
    Teuchos::RCP<LINALG::SparseOperator> temp_struct_scatra_mat_;

   private:
    //! indicates, if meshtying includes three entities
    bool Meshtying3DomainIntersection() const { return meshtying_3_domain_intersection_; }

   private:
    /*!
     * @brief finalize mesh tying by applying pseudo dirichlet conditions to the slave side, i.e.
     * write 1.0 on main diagonal of slave side dofs
     *
     * @param[in,out] ssi_structure_matrix  structure matrix with mesh tying constraints
     */
    void FinalizeMeshtyingStructureMatrix(LINALG::SparseMatrix& ssi_structure_matrix);

    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    //! indicates, if meshtying includes three entities
    const bool meshtying_3_domain_intersection_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps_;

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying_;
  };

  //! SSI problem is represented by one sparse matrix
  class MeshtyingStrategySparse : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategySparse(bool is_scatra_manifold, bool meshtying_3_domain_intersection,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix, bool do_uncomplete) override;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix, bool do_uncomplete) override;
  };

  //! SSI problem is composed of sub matrices
  class MeshtyingStrategyBlock : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBlock(bool is_scatra_manifold, bool meshtying_3_domain_intersection,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override = 0;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix,
        bool do_uncomplete) override = 0;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix,
        bool do_uncomplete) override = 0;

   protected:
    //! position of scatra blocks in system matrix
    Teuchos::RCP<std::vector<int>> BlockPositionScaTra() { return block_position_scatra_; }

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

  //! SSI problem is composed of sparse structure sub matrix and sparse scatra sub matrix
  class MeshtyingStrategyBlockSparse : public MeshtyingStrategyBlock
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBlockSparse(bool is_scatra_manifold,
        bool meshtying_3_domain_intersection, Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix, bool do_uncomplete) override;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix, bool do_uncomplete) override;
  };

  //! SSI problem is composed of sparse structure sub matrix and block scatra sub matrix
  class MeshtyingStrategyBlockBlock : public MeshtyingStrategyBlock
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBlockBlock(bool is_scatra_manifold,
        bool meshtying_3_domain_intersection, Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);

    void ApplyMeshtyingToScatraManifoldStructure(
        Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix, bool do_uncomplete) override;

    void ApplyMeshtyingToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix, bool do_uncomplete) override;
  };

  //! build specific mesh tying strategy
  Teuchos::RCP<SSI::MeshtyingStrategyBase> BuildMeshtyingStrategy(bool is_scatra_manifold,
      LINALG::MatrixType matrixtype_ssi, LINALG::MatrixType matrixtype_scatra,
      bool meshtying_3_domain_intersection, Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
      Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying);
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_MESHTYING_STRATEGY_H