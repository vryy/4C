/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy_poro.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mortar/mortar_element.H"
#include "contact_interface.H"
#include "contact_nitsche_utils.H"
#include "contact_paramsinterface.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_so3/so3_plast/so3_ssn_plast.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_cut/cut_output.H"
void CONTACT::CoNitscheStrategyPoro::ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
    const int iter, bool predictor)
{
  if (predictor) return;

  CONTACT::CoNitscheStrategy::ApplyForceStiffCmt(dis, kt, f, step, iter, predictor);

  // Evaluation for all interfaces
  fp_ = CreateRhsBlockPtr(DRT::UTILS::block_porofluid);
  kpp_ = CreateMatrixBlockPtr(DRT::UTILS::block_porofluid_porofluid);
  kpd_ = CreateMatrixBlockPtr(DRT::UTILS::block_porofluid_displ);
  kdp_ = CreateMatrixBlockPtr(DRT::UTILS::block_displ_porofluid);
  //    for (int i = 0; i < (int) interface_.size(); ++i)
  //    {
  //      for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
  //      {
  //        MORTAR::MortarElement* mele
  //        =dynamic_cast<MORTAR::MortarElement*>(interface_[i]->Discret().gElement(
  //            interface_[i]->Discret().ElementColMap()->GID(e)));
  //        mele->GetNitscheContainer().ClearAll();
  //      }
  //    }
}

void CONTACT::CoNitscheStrategyPoro::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_svelocity)
  {
    SetParentState(statename, vec);
  }
  else
    CONTACT::CoNitscheStrategy::SetState(statename, vec);
  return;
}

void CONTACT::CoNitscheStrategyPoro::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  //
  if (statename == MORTAR::state_fvelocity || statename == MORTAR::state_fpressure)
  {
    Teuchos::RCP<DRT::Discretization> dis = DRT::Problem::Instance()->GetDis("porofluid");
    if (dis == Teuchos::null) dserror("didn't get my discretization");

    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(), true));
    LINALG::Export(vec, *global);


    // set state on interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      DRT::Discretization& idiscret_ = interface_[i]->Discret();

      for (int j = 0; j < interface_[i]->Discret().ElementColMap()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->Discret().ElementColMap()->GID(j);

        MORTAR::MortarElement* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;

        if (ele->ParentSlaveElement())  // if this pointer is NULL, this parent is impermeable
        {
          // this gets values in local order
          ele->ParentSlaveElement()->LocationVector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*global, myval, lm);

          std::vector<double> vel;
          std::vector<double> pres;

          for (int n = 0; n < ele->ParentSlaveElement()->NumNode(); ++n)
          {
            for (uint dim = 0; dim < 3; ++dim)
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
    CONTACT::CoNitscheStrategy::SetParentState(statename, vec);

  return;
}

Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategyPoro::SetupRhsBlockVec(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  switch (bt)
  {
    case DRT::UTILS::block_porofluid:
      return Teuchos::rcp(
          new Epetra_FEVector(*DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap()));
      break;
    default:
      return CONTACT::CoNitscheStrategy::SetupRhsBlockVec(bt);
      break;
  }
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategyPoro::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bp) const
{
  if (curr_state_eval_ == false) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::block_porofluid:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fp_), 0));
    default:
      return CONTACT::CoNitscheStrategy::GetRhsBlockPtr(bp);
  }
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategyPoro::SetupMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt)
{
  switch (bt)
  {
    case DRT::UTILS::block_displ_porofluid:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("structure")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      break;
    case DRT::UTILS::block_porofluid_displ:
    case DRT::UTILS::block_porofluid_porofluid:
      return Teuchos::rcp(
          new LINALG::SparseMatrix(*Teuchos::rcpFromRef<const Epetra_Map>(
                                       *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap()),
              100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      break;
    default:
      return CONTACT::CoNitscheStrategy::SetupMatrixBlockPtr(bt);
      break;
  }
  return Teuchos::null;
}

void CONTACT::CoNitscheStrategyPoro::CompleteMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case DRT::UTILS::block_displ_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::block_porofluid_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::block_porofluid_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::CoNitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategyPoro::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bp) const
{
  if (curr_state_eval_ == false) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::block_porofluid_porofluid:
      return kpp_;
      break;
    case DRT::UTILS::block_porofluid_displ:
      return kpd_;
      break;
    case DRT::UTILS::block_displ_porofluid:
      return kdp_;
      break;
    default:
      return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bp);
      break;
  }
}
