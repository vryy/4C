/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy_tsi.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include "../linalg/linalg_sparsematrix.H"
#include "contact_interface.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mortar/mortar_element.H"
#include "contact_nitsche_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "contact_paramsinterface.H"
#include "../drt_so3/so3_plast/so3_ssn_plast.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_fsi/fsi_matrixtransform.H"



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

  return;
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
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      DRT::Discretization& idiscret_ = interface_[i]->Discret();

      for (int j = 0; j < interface_[i]->Discret().ElementColMap()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->Discret().ElementColMap()->GID(j);

        DRT::Element* e = idiscret_.gElement(gid);
        if (e == NULL) dserror("basic element not found");

        MORTAR::MortarElement* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));
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

  return;
}

void CONTACT::CoNitscheStrategyTsi::Setup(bool redistributed, bool init)
{
  CONTACT::CoNitscheStrategy::Setup(redistributed, init);

  curr_state_temp_ = Teuchos::null;
  return;
}

void CONTACT::CoNitscheStrategyTsi::UpdateTraceIneqEtimates()
{
  INPAR::CONTACT::NitscheWeighting NitWgt =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(Params(), "NITSCHE_WEIGHTING");
  for (int i = 0; i < (int)interface_.size(); ++i)
    for (int e = 0; e < interface_[i]->Discret().ElementColMap()->NumMyElements(); ++e)
    {
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(
          interface_[i]->Discret().gElement(interface_[i]->Discret().ElementColMap()->GID(e)));
      if (NitWgt == INPAR::CONTACT::NitWgt_slave && mele->IsSlave() == false) continue;
      if (NitWgt == INPAR::CONTACT::NitWgt_master && mele->IsSlave() == true) continue;
      mele->EstimateNitscheTraceMaxEigenvalueCombined();
    }
  return;
}

Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategyTsi::SetupRhsBlockVec(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  switch (bt)
  {
    case DRT::UTILS::block_temp:
      return Teuchos::rcp(
          new Epetra_FEVector(*DRT::Problem::Instance()->GetDis("thermo")->DofRowMap()));
      break;
    default:
      return CONTACT::CoNitscheStrategy::SetupRhsBlockVec(bt);
      break;
  }
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategyTsi::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  if (curr_state_eval_ == false) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case DRT::UTILS::block_temp:
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
    case DRT::UTILS::block_displ_temp:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("structure")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      break;
    case DRT::UTILS::block_temp_displ:
    case DRT::UTILS::block_temp_temp:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("thermo")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      break;
    default:
      return CONTACT::CoNitscheStrategy::SetupMatrixBlockPtr(bt);
      break;
  }
  return Teuchos::null;
}

void CONTACT::CoNitscheStrategyTsi::CompleteMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case DRT::UTILS::block_displ_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(*DRT::Problem::Instance()->GetDis("thermo")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),           // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::block_temp_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("thermo")->DofRowMap(),     // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::block_temp_temp:
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
  if (curr_state_eval_ == false) dserror("you didn't evaluate this contact state first");

  switch (bt)
  {
    case DRT::UTILS::block_temp_temp:
      return ktt_;
      break;
    case DRT::UTILS::block_temp_displ:
      return ktd_;
      break;
    case DRT::UTILS::block_displ_temp:
      return kdt_;
      break;
    default:
      return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bt, cparams);
      break;
  }
}


void CONTACT::CoNitscheStrategyTsi::Integrate(CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  ft_ = CreateRhsBlockPtr(DRT::UTILS::block_temp);
  ktt_ = CreateMatrixBlockPtr(DRT::UTILS::block_temp_temp);
  ktd_ = CreateMatrixBlockPtr(DRT::UTILS::block_temp_displ);
  kdt_ = CreateMatrixBlockPtr(DRT::UTILS::block_displ_temp);

  return;
}
