/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_DBC_HANDLER_H
#define SSI_MONOLITHIC_DBC_HANDLER_H

#include "Teuchos_RCPDecl.hpp"

#include "../linalg/linalg_sparseoperator.H"

namespace DRT
{
  namespace UTILS
  {
    class LocsysManager;
  }
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;

  //! base class of Dirichlet boundary condition handler for monolithic scatra-structure interaction
  class DBCHandlerBase
  {
   public:
    //! constructor
    explicit DBCHandlerBase(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

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
    void ApplyDBCToSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> system_matrix);

    /*!
     * @brief apply structure Dirichlet boundary conditions to system matrix
     *
     * @param[in,out] system_matrix  system matrix
     */
    void ApplyStructureDBCToSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> system_matrix);

   protected:
    //! return the monolithic algorithm for scalar-structure interaction
    Teuchos::RCP<const SSI::SSIMono> SSIMono() const { return ssi_mono_; }

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
        Teuchos::RCP<LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) = 0;

    //! monolithic algorithm for scalar-structure interaction
    const Teuchos::RCP<const SSI::SSIMono> ssi_mono_;
  };

  //! Dirichlet boundary condition handler class for monolithic scalar-structure interaction if the
  //! global system matrix is a sparse matrix
  class DBCHandlerSparse : public DBCHandlerBase
  {
   public:
    //! constructor
    explicit DBCHandlerSparse(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

   private:
    void ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) override;
  };

  //! Dirichlet boundary condition handler class for monolithic scalar-structure interaction if the
  //! global system matrix is a block matrix
  class DBCHandlerBlock : public DBCHandlerBase
  {
   public:
    //! constructor
    explicit DBCHandlerBlock(const Teuchos::RCP<const SSI::SSIMono> ssi_mono);

   private:
    void ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        Teuchos::RCP<LINALG::SparseOperator> system_matrix,
        const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
        Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure) override;

    //! position of structure block in system matrix
    int PositionStructure() const { return position_structure_; };

    //! position of structure block in system matrix
    int position_structure_;
  };

  /*!
   * @brief build specific Dirichlet boundary condition handler
   *
   * @param[in] ssi_mono        monolithic algorithm for scalar-structure interaction
   * @param[in] matrixtype_ssi  matrix type of scalar-structure interaction system matrix
   * @return Dirichlet boundary condition handler
   */
  Teuchos::RCP<SSI::DBCHandlerBase> BuildDBCHandler(
      Teuchos::RCP<const SSI::SSIMono> ssi_mono, LINALG::MatrixType matrixtype_ssi);
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_DBC_HANDLER_H
