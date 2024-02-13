/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef BACI_SSI_MONOLITHIC_DBC_HANDLER_HPP
#define BACI_SSI_MONOLITHIC_DBC_HANDLER_HPP

#include "baci_config.hpp"

#include "baci_linalg_sparseoperator.hpp"

#include <Teuchos_RCPDecl.hpp>

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class SSIStructureWrapper;
}

namespace DRT
{
  namespace UTILS
  {
    class LocsysManager;
  }
}  // namespace DRT

namespace CORE::LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace SSI
{
  class SSIMono;

  namespace UTILS
  {
    class SSIMaps;
  }

  //! base class of Dirichlet boundary condition handler for monolithic scatra-structure interaction
  class DBCHandlerBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~DBCHandlerBase() = default;

    //! constructor
    explicit DBCHandlerBase(bool is_scatra_manifold, Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

    /*!
     * @brief apply Dirichlet boundary conditions to global right hand side vector
     *
     * @param[in,out] rhs  global right hand side vector
     */
    void ApplyDBCToRHS(Teuchos::RCP<Epetra_Vector> rhs);

    /*!
     * @brief apply Dirichlet boundary conditions to global system matrix
     *
     * @param[in,out] system_matrix  system matrix
     */
    void ApplyDBCToSystemMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix);

    /*!
     * @brief apply structure Dirichlet boundary conditions to system matrix
     *
     * @param[in,out] system_matrix  system matrix
     */
    void ApplyStructureDBCToSystemMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix);

   protected:
    //! solve additional scatra field on manifolds
    bool IsScaTraManifold() const { return is_scatra_manifold_; }

    //! access to scalar transport field
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const { return scatra_; }

    //! access to scalar transport on manifold field
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraManifoldField() const { return scatra_manifold_; }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

    //! access to structural field
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> StructureField() const { return structure_; }

   private:
    /*!
     * @brief apply structure Dirichlet boundary conditions to system matrix in case a rotation of
     * the coordinate system has been defined
     *
     * @param[in,out] system_matrix        system matrix
     * @param[in] dbcmap_structure         Dirichlet boundary condition map of structure problem
     * @param[in] locsysmanager_structure  manager of the local coordinate system for the structure
     *                                     problem
     */
    virtual void ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) = 0;

    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    /// scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    /// scalar transport on manifold field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps_;

    /// structure field solver
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure_;
  };

  //! Dirichlet boundary condition handler class for monolithic scalar-structure interaction if the
  //! global system matrix is a sparse matrix
  class DBCHandlerSparse : public DBCHandlerBase
  {
   public:
    //! constructor
    explicit DBCHandlerSparse(bool is_scatra_manifold,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

   private:
    void ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) override;
  };

  //! Dirichlet boundary condition handler class for monolithic scalar-structure interaction if the
  //! global system matrix is a block matrix
  class DBCHandlerBlock : public DBCHandlerBase
  {
   public:
    //! constructor
    explicit DBCHandlerBlock(bool is_scatra_manifold, Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
        Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);

   private:
    void ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) override;

    //! position of structure block in system matrix
    int PositionStructure() const { return position_structure_; };

    //! position of structure block in system matrix
    const int position_structure_;
  };

  /*!
   * @brief build specific Dirichlet boundary condition handler
   *
   * @param[in] ssi_mono        monolithic algorithm for scalar-structure interaction
   * @param[in] matrixtype_ssi  matrix type of scalar-structure interaction system matrix
   * @return Dirichlet boundary condition handler
   */
  Teuchos::RCP<SSI::DBCHandlerBase> BuildDBCHandler(bool is_scatra_manifold,
      CORE::LINALG::MatrixType matrixtype_ssi, Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
      Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
      Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure);
}  // namespace SSI
BACI_NAMESPACE_CLOSE

#endif  // SSI_MONOLITHIC_DBC_HANDLER_H
