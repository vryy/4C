/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "contact_nitsche_strategy_ssi.H"
#include "contact_interface.H"

#include "lib_globalproblem.H"

#include "mortar_element.H"

#include "linalg_utils_sparse_algebra_manipulation.H"

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  fs_ = CreateRhsBlockPtr(DRT::UTILS::VecBlockType::scatra);
  kss_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra);
  ksd_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ);
  kds_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::EvaluateReferenceState()
{
  // initialize an estimate of TraceHE
  InitTraceHE();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::InitTraceHE()
{
  for (const auto& interface : interface_)
  {
    for (int e = 0; e < interface->Discret().ElementColMap()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<MORTAR::MortarElement*>(
          interface->Discret().gElement(interface->Discret().ElementColMap()->GID(e)));
      // set the initial TraceHE to the edge length size of the element
      mele->TraceHE() = mele->MaxEdgeSize();
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  switch (statename)
  {
    case MORTAR::state_elch:
    case MORTAR::state_scalar:
    {
      double inf_delta = 0.0;
      if (curr_state_scalar_ == Teuchos::null)
      {
        curr_state_scalar_ = Teuchos::rcp(new Epetra_Vector(vec));
        inf_delta = 1.0e12;
      }
      else
      {
        auto delta(vec);
        delta.Update(-1.0, *curr_state_scalar_, 1.0);
        delta.NormInf(&inf_delta);
      }
      if (inf_delta < 1.0e-16)
        return;
      else
      {
        curr_state_eval_ = false;
        *curr_state_scalar_ = vec;
        SetParentState(statename, vec);
      }
      break;
    }
    default:
    {
      CONTACT::CoNitscheStrategy::SetState(statename, vec);
      break;
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  switch (statename)
  {
    case MORTAR::state_elch:
    case MORTAR::state_scalar:
    {
      auto scatra_dis = DRT::Problem::Instance()->GetDis("scatra");
      if (scatra_dis == Teuchos::null) dserror("didn't get scatra discretization");

      auto scatra_dofcolmap = Teuchos::rcp(new Epetra_Vector(*scatra_dis->DofColMap(), true));
      CORE::LINALG::Export(vec, *scatra_dofcolmap);

      // set state on interfaces
      for (const auto& interface : interface_)
      {
        // get the interface discretization
        const auto& interface_dis = interface->Discret();

        // loop over all interface column elements owned by current processor
        for (int j = 0; j < interface_dis.ElementColMap()->NumMyElements(); ++j)
        {
          const int gid = interface_dis.ElementColMap()->GID(j);

          auto* interface_ele = interface_dis.gElement(gid);
          if (interface_ele == nullptr) dserror("Did not find element.");

          auto* mortar_ele = dynamic_cast<MORTAR::MortarElement*>(interface_ele);
          auto* mortar_parent_ele = scatra_dis->gElement(mortar_ele->ParentElementId());

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          if (mortar_parent_ele == nullptr)
            dserror("Did not get parent element to extract scalar values");
          else
            mortar_parent_ele->LocationVector(*scatra_dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*scatra_dofcolmap, myval, lm);

          mortar_ele->MoData().ParentScalar() = myval;
          mortar_ele->MoData().ParentScalarDof() = lm;
        }
      }
      break;
    }
    default:
    {
      CONTACT::CoNitscheStrategy::SetParentState(statename, vec);
      break;
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategySsi::SetupRhsBlockVec(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  switch (bt)
  {
    case DRT::UTILS::VecBlockType::elch:
    case DRT::UTILS::VecBlockType::scatra:
      return Teuchos::rcp(
          new Epetra_FEVector(*DRT::Problem::Instance()->GetDis("scatra")->DofRowMap()));
    default:
      return CONTACT::CoNitscheStrategy::SetupRhsBlockVec(bt);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategySsi::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::VecBlockType::elch:
    case DRT::UTILS::VecBlockType::scatra:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fs_), 0));
    default:
      return CONTACT::CoNitscheStrategy::GetRhsBlockPtr(bp);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::CoNitscheStrategySsi::SetupMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_elch:
    case DRT::UTILS::MatBlockType::displ_scatra:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *DRT::Problem::Instance()->GetDis("structure")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    case DRT::UTILS::MatBlockType::elch_displ:
    case DRT::UTILS::MatBlockType::elch_elch:
    case DRT::UTILS::MatBlockType::scatra_displ:
    case DRT::UTILS::MatBlockType::scatra_scatra:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *DRT::Problem::Instance()->GetDis("scatra")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::CoNitscheStrategy::SetupMatrixBlockPtr(bt);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsi::CompleteMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_elch:
    case DRT::UTILS::MatBlockType::displ_scatra:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(*DRT::Problem::Instance()->GetDis("scatra")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),           // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::elch_displ:
    case DRT::UTILS::MatBlockType::scatra_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("scatra")->DofRowMap(),     // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::elch_elch:
    case DRT::UTILS::MatBlockType::scatra_scatra:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::CoNitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::CoNitscheStrategySsi::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::MatBlockType::elch_elch:
    case DRT::UTILS::MatBlockType::scatra_scatra:
      return kss_;
    case DRT::UTILS::MatBlockType::elch_displ:
    case DRT::UTILS::MatBlockType::scatra_displ:
      return ksd_;
    case DRT::UTILS::MatBlockType::displ_elch:
    case DRT::UTILS::MatBlockType::displ_scatra:
      return kds_;
    default:
      return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bp, nullptr);
  }
}