/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_poro.hpp"

#include "baci_contact_interface.hpp"
#include "baci_contact_nitsche_utils.hpp"
#include "baci_contact_paramsinterface.hpp"
#include "baci_coupling_adapter.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_so3_plast_ssn.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Operator.h>

BACI_NAMESPACE_OPEN

void CONTACT::NitscheStrategyPoro::ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
    const int iter, bool predictor)
{
  if (predictor) return;

  CONTACT::NitscheStrategy::ApplyForceStiffCmt(dis, kt, f, step, iter, predictor);

  // Evaluation for all interfaces
  fp_ = CreateRhsBlockPtr(CONTACT::VecBlockType::porofluid);
  kpp_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::porofluid_porofluid);
  kpd_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::porofluid_displ);
  kdp_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::displ_porofluid);
  //    for (int i = 0; i < (int) interface_.size(); ++i)
  //    {
  //      for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
  //      {
  //        MORTAR::Element* mele
  //        =dynamic_cast<MORTAR::Element*>(interface_[i]->Discret().gElement(
  //            interface_[i]->Discret().ElementColMap()->GID(e)));
  //        mele->GetNitscheContainer().ClearAll();
  //      }
  //    }
}

void CONTACT::NitscheStrategyPoro::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_svelocity)
  {
    SetParentState(statename, vec);
  }
  else
    CONTACT::NitscheStrategy::SetState(statename, vec);
}

void CONTACT::NitscheStrategyPoro::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  //
  if (statename == MORTAR::state_fvelocity || statename == MORTAR::state_fpressure)
  {
    Teuchos::RCP<DRT::Discretization> dis = GLOBAL::Problem::Instance()->GetDis("porofluid");
    if (dis == Teuchos::null) dserror("didn't get my discretization");

    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(), true));
    CORE::LINALG::Export(vec, *global);


    // set state on interfaces
    for (const auto& interface : interface_)
    {
      DRT::Discretization& idiscret = interface->Discret();

      for (int j = 0; j < interface->Discret().ElementColMap()->NumMyElements(); ++j)
      {
        const int gid = interface->Discret().ElementColMap()->GID(j);

        auto* ele = dynamic_cast<MORTAR::Element*>(idiscret.gElement(gid));

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;

        if (ele->ParentSlaveElement())  // if this pointer is nullptr, this parent is impermeable
        {
          // this gets values in local order
          ele->ParentSlaveElement()->LocationVector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*global, myval, lm);

          std::vector<double> vel;
          std::vector<double> pres;

          for (int n = 0; n < ele->ParentSlaveElement()->NumNode(); ++n)
          {
            for (unsigned dim = 0; dim < 3; ++dim)
            {
              vel.push_back(myval[n * 4 + dim]);
            }
            pres.push_back(myval[n * 4 + 3]);
          }

          ele->MoData().ParentPFPres() = pres;
          ele->MoData().ParentPFVel() = vel;
          ele->MoData().ParentPFDof() = lm;
        }
      }
    }
  }
  else
    CONTACT::NitscheStrategy::SetParentState(statename, vec);
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategyPoro::SetupRhsBlockVec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::porofluid:
      return Teuchos::rcp(
          new Epetra_FEVector(*GLOBAL::Problem::Instance()->GetDis("porofluid")->DofRowMap()));
    default:
      return CONTACT::NitscheStrategy::SetupRhsBlockVec(bt);
  }
}

Teuchos::RCP<const Epetra_Vector> CONTACT::NitscheStrategyPoro::GetRhsBlockPtr(
    const enum CONTACT::VecBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case CONTACT::VecBlockType::porofluid:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fp_), 0));
    default:
      return CONTACT::NitscheStrategy::GetRhsBlockPtr(bp);
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::NitscheStrategyPoro::SetupMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_porofluid:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    case CONTACT::MatBlockType::porofluid_displ:
    case CONTACT::MatBlockType::porofluid_porofluid:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *GLOBAL::Problem::Instance()->GetDis("porofluid")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::NitscheStrategy::SetupMatrixBlockPtr(bt);
  }
}

void CONTACT::NitscheStrategyPoro::CompleteMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *GLOBAL::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // col map
                  *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::porofluid_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *GLOBAL::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::porofluid_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::NitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::NitscheStrategyPoro::GetMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case CONTACT::MatBlockType::porofluid_porofluid:
      return kpp_;
    case CONTACT::MatBlockType::porofluid_displ:
      return kpd_;
    case CONTACT::MatBlockType::displ_porofluid:
      return kdp_;
    default:
      return CONTACT::NitscheStrategy::GetMatrixBlockPtr(bp, nullptr);
  }
}

BACI_NAMESPACE_CLOSE
