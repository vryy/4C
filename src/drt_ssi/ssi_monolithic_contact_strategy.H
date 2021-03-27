/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_CONTACT_STRATEGY_H
#define SSI_MONOLITHIC_CONTACT_STRATEGY_H

#include "Epetra_Vector.h"
#include "Teuchos_RCPDecl.hpp"

namespace LINALG
{
  enum class MatrixType;
  class SparseOperator;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;

  //! base functionality for scatra structure contact interaction
  class ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategyBase(const SSI::SSIMono& ssi_mono);

    //! apply contact contributions to the scatra residual
    void ApplyContactToScatraResidual(Teuchos::RCP<Epetra_Vector> scatra_residual);

    //! apply contact contributions to scatra sub matrix
    virtual void ApplyContactToScatraScatra(
        Teuchos::RCP<LINALG::SparseOperator> scatra_scatra_matrix) = 0;

    //! apply contact contributions to scatra-structure sub matrix
    virtual void ApplyContactToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix) = 0;

    //! apply contact contributions to structure-scatra sub matrix
    virtual void ApplyContactToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix) = 0;

   protected:
    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& SSIMono() const { return ssi_mono_; }

   private:
    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& ssi_mono_;
  };

  //! SSI (sub) matrices are sparse matrices
  class ContactStrategySparse : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategySparse(const SSI::SSIMono& ssi_mono);

    void ApplyContactToScatraScatra(
        Teuchos::RCP<LINALG::SparseOperator> scatra_scatra_matrix) override;

    void ApplyContactToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix) override;

    void ApplyContactToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! SSI (sub) matrices are block matrices
  class ContactStrategyBlock : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategyBlock(const SSI::SSIMono& ssi_mono);

    void ApplyContactToScatraScatra(
        Teuchos::RCP<LINALG::SparseOperator> scatra_scatra_matrix) override;

    void ApplyContactToScatraStructure(
        Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix) override;

    void ApplyContactToStructureScatra(
        Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! build specific contact strategy
  Teuchos::RCP<SSI::ContactStrategyBase> BuildContactStrategy(
      const SSI::SSIMono& ssi_mono, LINALG::MatrixType matrixtype_scatra);
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_CONTACT_STRATEGY_H