/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_TSI_H
#define CONTACT_NITSCHE_STRATEGY_TSI_H

#include <utility>

#include "contact_nitsche_strategy.H"

namespace ADAPTER
{
  class Coupling;
}

namespace CONTACT
{
  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   For a more general documentation of the involved functions refer to CoAbstract Strategy.

   */
  class CoNitscheStrategyTsi : public CoNitscheStrategy
  {
   public:
    //! Standard constructor
    CoNitscheStrategyTsi(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : CoNitscheStrategy(
              DofRowMap, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    //! Shared data constructor
    CoNitscheStrategyTsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : CoNitscheStrategy(data_ptr, DofRowMap, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor) override
    {
      dserror("not implemented");
    }

    void Integrate(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     derived from contact abstract strategy.
     The Nitsche strategy does not have
      */
    void Setup(bool redistributed, bool init) override;

    void UpdateTraceIneqEtimates() override;

    ~CoNitscheStrategyTsi() override = default;

    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    void SetParentState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(const enum DRT::UTILS::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    //! [derived]
    bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) override
    {
      if (fix_redistribution_) return false;
      return CONTACT::CoAbstractStrategy::RedistributeContact(dis, vel);
    }

    void EnableRedistribution() { fix_redistribution_ = false; }

    // don't want = operator and cctor
    CoNitscheStrategyTsi operator=(const CoNitscheStrategyTsi& old) = delete;
    CoNitscheStrategyTsi(const CoNitscheStrategyTsi& old) = delete;

   protected:
    // create an appropriate vector for the RHS
    Teuchos::RCP<Epetra_FEVector> SetupRhsBlockVec(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    // create an appropriate matrix block
    Teuchos::RCP<LINALG::SparseMatrix> SetupMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) override;

    // complete matrix block with correct maps
    void CompleteMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<LINALG::SparseMatrix> kc) override;

    // do not reditribute (during constructor phase)
    bool fix_redistribution_;

    Teuchos::RCP<Epetra_Vector> curr_state_temp_;

    Teuchos::RCP<Epetra_FEVector> ft_;
    Teuchos::RCP<LINALG::SparseMatrix> ktt_;
    Teuchos::RCP<LINALG::SparseMatrix> ktd_;
    Teuchos::RCP<LINALG::SparseMatrix> kdt_;
  };
}  // namespace CONTACT
#endif  //  CONTACT_NITSCHE_STRATEGY_TSI_H
