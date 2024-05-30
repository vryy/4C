/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic/partitioning SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_CONTACT_STRATEGY_HPP
#define FOUR_C_SSI_CONTACT_STRATEGY_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class NitscheStrategySsi;
}

namespace CORE::LINALG
{
  enum class MatrixType;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace SSI
{
  class SsiMono;

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
    explicit ContactStrategyBase(Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    //! apply contact contributions to the scatra residual
    void apply_contact_to_scatra_residual(Teuchos::RCP<Epetra_Vector> scatra_residual);

    //! apply contact contributions to scatra sub matrix
    virtual void apply_contact_to_scatra_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) = 0;

    //! apply contact contributions to scatra-structure sub matrix
    virtual void apply_contact_to_scatra_structure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) = 0;

    //! apply contact contributions to structure-scatra sub matrix
    virtual void apply_contact_to_structure_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) = 0;

   protected:
    //! return contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::NitscheStrategySsi> nitsche_strategy_ssi() const
    {
      return contact_strategy_nitsche_;
    }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps() const { return ssi_maps_; }

   private:
    //! store contact nitsche strategy for ssi problems
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_strategy_nitsche_;

    //! this object holds all maps relevant to monolithic/partitioning scalar transport - structure
    //! interaction
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps_;
  };

  //! SSI (sub) matrices are sparse matrices
  class ContactStrategySparse : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategySparse(
        Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    void apply_contact_to_scatra_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) override;

    void apply_contact_to_scatra_structure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) override;

    void apply_contact_to_structure_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! SSI (sub) matrices are block matrices
  class ContactStrategyBlock : public ContactStrategyBase
  {
   public:
    //! constructor
    explicit ContactStrategyBlock(
        Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
        Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

    void apply_contact_to_scatra_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix) override;

    void apply_contact_to_scatra_structure(
        Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix) override;

    void apply_contact_to_structure_scatra(
        Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix) override;
  };

  //! build specific contact strategy
  Teuchos::RCP<SSI::ContactStrategyBase> BuildContactStrategy(
      Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
      Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, CORE::LINALG::MatrixType matrixtype_scatra);
}  // namespace SSI

FOUR_C_NAMESPACE_CLOSE

#endif
