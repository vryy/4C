/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_STRATEGY_TSI_H
#define BACI_CONTACT_NITSCHE_STRATEGY_TSI_H

#include "baci_config.H"

#include "baci_contact_nitsche_strategy.H"

#include <utility>

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class Coupling;
}

namespace CONTACT
{
  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class NitscheStrategyTsi : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategyTsi(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(
              DofRowMap, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    //! Shared data constructor
    NitscheStrategyTsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, DofRowMap, NodeRowMap, params, std::move(interface), dim, comm,
              alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override
    {
      dserror("not implemented");
    }

    void Integrate(const CONTACT::ParamsInterface& cparams) override;

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     derived from contact abstract strategy.
     The Nitsche strategy does not have
      */
    void Setup(bool redistributed, bool init) override;

    void UpdateTraceIneqEtimates() override;
    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    void SetParentState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const override;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    //! [derived]
    bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) override
    {
      if (fix_redistribution_) return false;
      return CONTACT::AbstractStrategy::RedistributeContact(dis, vel);
    }

    void EnableRedistribution() { fix_redistribution_ = false; }

    // don't want = operator and cctor
    NitscheStrategyTsi operator=(const NitscheStrategyTsi& old) = delete;
    NitscheStrategyTsi(const NitscheStrategyTsi& old) = delete;

   protected:
    // create an appropriate vector for the RHS
    Teuchos::RCP<Epetra_FEVector> SetupRhsBlockVec(
        const enum CONTACT::VecBlockType& bt) const override;

    // create an appropriate matrix block
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SetupMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bt) override;

    // complete matrix block with correct maps
    void CompleteMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc) override;

    // do not reditribute (during constructor phase)
    bool fix_redistribution_;

    Teuchos::RCP<Epetra_Vector> curr_state_temp_;

    Teuchos::RCP<Epetra_FEVector> ft_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ktt_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ktd_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kdt_;
  };
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_STRATEGY_TSI_H
