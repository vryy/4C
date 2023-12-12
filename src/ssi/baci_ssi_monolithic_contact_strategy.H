/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef BACI_SSI_MONOLITHIC_CONTACT_STRATEGY_H
#define BACI_SSI_MONOLITHIC_CONTACT_STRATEGY_H

#include <Epetra_Vector.h>
#include <Teuchos_RCPDecl.hpp>

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  class CoNitscheStrategySsi;
}

namespace CORE::LINALG
{
  enum class MatrixType;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace SSI
{
  class SSIMono;

  namespace UTILS
  {
    class SSIMaps;
  }

  //! base functionality for scatra structure contact interaction
  class ContactStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~ContactStrategyBase() = default;

    //! constructor
    explicit ContactStrategyBase(
        Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    //! apply contact contributions to the scatra residual
    void ApplyContactToScatraResidual(Teuchos::RCP<Epetra_Vector> scatra_residual);

    //! apply contact contributions to scatra sub matrix
    virtual void ApplyContactToScatraScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) = 0;

    //! apply contact contributions to scatra-structure sub matrix
    virtual void ApplyContactToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) = 0;

    //! apply contact contributions to structure-scatra sub matrix
    virtual void ApplyContactToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) = 0;

   protected:
    //! return contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> CoNitscheStrategySsi() const
    {
      return contact_strategy_nitsche_;
    }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> SSIMaps() const { return ssi_maps_; }

   private:
    //! store contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_strategy_nitsche_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps_;
  };

  //! SSI (sub) matrices are sparse matrices
  class ContactStrategySparse : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategySparse(
        Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    void ApplyContactToScatraScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) override;

    void ApplyContactToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) override;

    void ApplyContactToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! SSI (sub) matrices are block matrices
  class ContactStrategyBlock : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategyBlock(
        Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    void ApplyContactToScatraScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) override;

    void ApplyContactToScatraStructure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) override;

    void ApplyContactToStructureScatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! build specific contact strategy
  Teuchos::RCP<SSI::ContactStrategyBase> BuildContactStrategy(
      Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
      Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, CORE::LINALG::MatrixType matrixtype_scatra);
}  // namespace SSI
BACI_NAMESPACE_CLOSE

#endif  // SSI_MONOLITHIC_CONTACT_STRATEGY_H
