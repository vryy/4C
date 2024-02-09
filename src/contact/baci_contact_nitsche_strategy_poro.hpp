/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_STRATEGY_PORO_HPP
#define BACI_CONTACT_NITSCHE_STRATEGY_PORO_HPP

#include "baci_config.hpp"

#include "baci_contact_nitsche_strategy.hpp"

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
  class NitscheStrategyPoro : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategyPoro(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::Interface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(
              DofRowMap, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          no_penetration_(params.get<bool>("CONTACTNOPEN"))
    {
    }

    //! Shared data constructor
    NitscheStrategyPoro(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, DofRowMap, NodeRowMap, params, std::move(interface), dim, comm,
              alphaf, maxdof),
          no_penetration_(params.get<bool>("CONTACTNOPEN"))
    {
    }

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override;

    //  void Integrate(CONTACT::ParamsInterface& cparams);
    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    void SetParentState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bp) const override;

    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bp) const;


    // Flag for Poro No Penetration Condition
    bool HasPoroNoPenetration() const override { return no_penetration_; }

    // don't want = operator and cctor
    NitscheStrategyPoro operator=(const NitscheStrategyPoro& old) = delete;
    NitscheStrategyPoro(const NitscheStrategyPoro& old) = delete;

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

    bool no_penetration_;

    Teuchos::RCP<Epetra_FEVector> fp_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kpp_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kpd_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kdp_;
  };
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_STRATEGY_PORO_H
