/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_tsi.hpp"

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

void CONTACT::NitscheStrategyTsi::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_temperature)
  {
    double inf_delta = 0.;
    if (curr_state_temp_ == Teuchos::null)
    {
      curr_state_temp_ = Teuchos::rcp(new Epetra_Vector(vec));
      inf_delta = 1.e12;
    }
    else
    {
      Epetra_Vector delta(vec);
      delta.Update(-1., *curr_state_temp_, 1.);
      delta.NormInf(&inf_delta);
    }
    if (inf_delta < 1.e-16)
      return;
    else
    {
      curr_state_eval_ = false;
      (*curr_state_temp_) = vec;
      SetParentState(statename, vec);
    }
  }
  else
    CONTACT::NitscheStrategy::SetState(statename, vec);
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::NitscheStrategyTsi::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_temperature)
  {
    Teuchos::RCP<DRT::Discretization> disT = GLOBAL::Problem::Instance()->GetDis("thermo");
    if (disT.is_null()) dserror("set state temperature, but no thermo-discretization???");

    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*disT->DofColMap(), true));
    CORE::LINALG::Export(vec, *global);

    // set state on interfaces
    for (auto& interface : interface_)
    {
      DRT::Discretization& idiscret = interface->Discret();

      for (int j = 0; j < idiscret.ElementColMap()->NumMyElements(); ++j)
      {
        int gid = idiscret.ElementColMap()->GID(j);

        DRT::Element* e = idiscret.gElement(gid);
        if (e == nullptr) dserror("basic element not found");

        auto* ele = dynamic_cast<MORTAR::Element*>(idiscret.gElement(gid));
        DRT::Element* ele_parentT = disT->gElement(ele->ParentElementId());

        std::vector<int> lm, lmowner, lmstride;
        ele_parentT->LocationVector(*disT, lm, lmowner, lmstride);

        std::vector<double> myval;
        DRT::UTILS::ExtractMyValues(*global, myval, lm);

        ele->MoData().ParentTemp() = myval;
        ele->MoData().ParentTempDof() = lm;
      }
    }
  }
  else
    CONTACT::NitscheStrategy::SetParentState(statename, vec);
}

void CONTACT::NitscheStrategyTsi::Setup(bool redistributed, bool init)
{
  CONTACT::NitscheStrategy::Setup(redistributed, init);

  curr_state_temp_ = Teuchos::null;
}

void CONTACT::NitscheStrategyTsi::UpdateTraceIneqEtimates()
{
  auto NitWgt =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(Params(), "NITSCHE_WEIGHTING");
  for (auto& interface : interface_)
  {
    for (int e = 0; e < interface->Discret().ElementColMap()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<MORTAR::Element*>(
          interface->Discret().gElement(interface->Discret().ElementColMap()->GID(e)));
      if (NitWgt == INPAR::CONTACT::NitWgt_slave && !mele->IsSlave()) continue;
      if (NitWgt == INPAR::CONTACT::NitWgt_master && mele->IsSlave()) continue;
      mele->EstimateNitscheTraceMaxEigenvalueCombined();
    }
  }
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategyTsi::SetupRhsBlockVec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::temp:
      return Teuchos::rcp(
          new Epetra_FEVector(*GLOBAL::Problem::Instance()->GetDis("thermo")->DofRowMap()));
    default:
      return CONTACT::NitscheStrategy::SetupRhsBlockVec(bt);
  }
}

Teuchos::RCP<const Epetra_Vector> CONTACT::NitscheStrategyTsi::GetRhsBlockPtr(
    const enum CONTACT::VecBlockType& bt) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case CONTACT::VecBlockType::temp:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(ft_), 0));
    default:
      return CONTACT::NitscheStrategy::GetRhsBlockPtr(bt);
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::NitscheStrategyTsi::SetupMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_temp:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    case CONTACT::MatBlockType::temp_displ:
    case CONTACT::MatBlockType::temp_temp:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *GLOBAL::Problem::Instance()->GetDis("thermo")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::NitscheStrategy::SetupMatrixBlockPtr(bt);
  }
}

void CONTACT::NitscheStrategyTsi::CompleteMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *GLOBAL::Problem::Instance()->GetDis("thermo")->DofRowMap(),     // col map
                  *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::temp_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *GLOBAL::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *GLOBAL::Problem::Instance()->GetDis("thermo")->DofRowMap(),     // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::temp_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::NitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::NitscheStrategyTsi::GetMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case CONTACT::MatBlockType::temp_temp:
      return ktt_;
    case CONTACT::MatBlockType::temp_displ:
      return ktd_;
    case CONTACT::MatBlockType::displ_temp:
      return kdt_;
    default:
      return CONTACT::NitscheStrategy::GetMatrixBlockPtr(bt, cparams);
  }
}


void CONTACT::NitscheStrategyTsi::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::NitscheStrategy::Integrate(cparams);

  ft_ = CreateRhsBlockPtr(CONTACT::VecBlockType::temp);
  ktt_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::temp_temp);
  ktd_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::temp_displ);
  kdt_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::displ_temp);
}

BACI_NAMESPACE_CLOSE
