/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy_tsi.H"

#include "contact_interface.H"
#include "contact_nitsche_utils.H"
#include "contact_paramsinterface.H"

#include <Epetra_FEVector.h>
#include <Epetra_Operator.h>

#include "adapter_coupling.H"

#include "lib_discret.H"
#include "lib_globalproblem.H"

#include "so3_plast_ssn.H"

#include "linalg_utils_sparse_algebra_manipulation.H"

void CONTACT::CoNitscheStrategyTsi::SetState(
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
    CONTACT::CoNitscheStrategy::SetState(statename, vec);
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategyTsi::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_temperature)
  {
    Teuchos::RCP<DRT::Discretization> disT = DRT::Problem::Instance()->GetDis("thermo");
    if (disT.is_null()) dserror("set state temperature, but no thermo-discretization???");

    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*disT->DofColMap(), true));
    LINALG::Export(vec, *global);

    // set state on interfaces
    for (auto& interface : interface_)
    {
      DRT::Discretization& idiscret = interface->Discret();

      for (int j = 0; j < idiscret.ElementColMap()->NumMyElements(); ++j)
      {
        int gid = idiscret.ElementColMap()->GID(j);

        DRT::Element* e = idiscret.gElement(gid);
        if (e == nullptr) dserror("basic element not found");

        auto* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(gid));
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
    CONTACT::CoNitscheStrategy::SetParentState(statename, vec);
}

void CONTACT::CoNitscheStrategyTsi::Setup(bool redistributed, bool init)
{
  CONTACT::CoNitscheStrategy::Setup(redistributed, init);

  curr_state_temp_ = Teuchos::null;
}

void CONTACT::CoNitscheStrategyTsi::UpdateTraceIneqEtimates()
{
  auto NitWgt =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(Params(), "NITSCHE_WEIGHTING");
  for (auto& interface : interface_)
  {
    for (int e = 0; e < interface->Discret().ElementColMap()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<MORTAR::MortarElement*>(
          interface->Discret().gElement(interface->Discret().ElementColMap()->GID(e)));
      if (NitWgt == INPAR::CONTACT::NitWgt_slave && !mele->IsSlave()) continue;
      if (NitWgt == INPAR::CONTACT::NitWgt_master && mele->IsSlave()) continue;
      mele->EstimateNitscheTraceMaxEigenvalueCombined();
    }
  }
}

Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategyTsi::SetupRhsBlockVec(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  switch (bt)
  {
    case DRT::UTILS::VecBlockType::temp:
      return Teuchos::rcp(
          new Epetra_FEVector(*DRT::Problem::Instance()->GetDis("thermo")->DofRowMap()));
    default:
      return CONTACT::CoNitscheStrategy::SetupRhsBlockVec(bt);
  }
}

Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategyTsi::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case DRT::UTILS::VecBlockType::temp:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(ft_), 0));
    default:
      return CONTACT::CoNitscheStrategy::GetRhsBlockPtr(bt);
  }
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategyTsi::SetupMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_temp:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("structure")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
    case DRT::UTILS::MatBlockType::temp_displ:
    case DRT::UTILS::MatBlockType::temp_temp:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("thermo")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::CoNitscheStrategy::SetupMatrixBlockPtr(bt);
  }
}

void CONTACT::CoNitscheStrategyTsi::CompleteMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(*DRT::Problem::Instance()->GetDis("thermo")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),           // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::temp_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("thermo")->DofRowMap(),     // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::temp_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::CoNitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategyTsi::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case DRT::UTILS::MatBlockType::temp_temp:
      return ktt_;
    case DRT::UTILS::MatBlockType::temp_displ:
      return ktd_;
    case DRT::UTILS::MatBlockType::displ_temp:
      return kdt_;
    default:
      return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bt, cparams);
  }
}


void CONTACT::CoNitscheStrategyTsi::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  ft_ = CreateRhsBlockPtr(DRT::UTILS::VecBlockType::temp);
  ktt_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::temp_temp);
  ktd_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::temp_displ);
  kdt_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_temp);
}
