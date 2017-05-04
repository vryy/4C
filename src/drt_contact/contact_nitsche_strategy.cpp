/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_strategy.cpp

\brief Nitsche contact solving strategy

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/

#include "contact_nitsche_strategy.H"
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


/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator     seitz 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategy::ApplyForceStiffCmt(
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
    Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter,
    bool predictor)
{
  // mortar initialization and evaluation
  SetState(MORTAR::state_new_displacement, *dis);

  // just a Nitsche-version
  Teuchos::RCP<Epetra_FEVector> fc=Teuchos::rcp(new Epetra_FEVector(f->Map()));
  Teuchos::RCP<LINALG::SparseMatrix> kc =
      Teuchos::rcp(new LINALG::SparseMatrix(
          (dynamic_cast<Epetra_CrsMatrix*>(&(*kt->EpetraOperator())))->RowMap(),
          100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // Evaluation for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    interface_[i]->Initialize();
    interface_[i]->Evaluate(0,step_,iter_);
    for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
    {
      MORTAR::MortarElement* mele =dynamic_cast<MORTAR::MortarElement*>(interface_[i]->Discret().gElement(
          interface_[i]->Discret().ElementColMap()->GID(e)));
      mele->GetNitscheContainer().AssembleRHS(mele,DRT::UTILS::block_displ,fc);
      mele->GetNitscheContainer().AssembleMatrix(mele,DRT::UTILS::block_displ_displ,kc);
      mele->GetNitscheContainer().Clear();
    }
  }
  if(fc->GlobalAssemble(Add,false)!=0) dserror("GlobalAssemble failed");
  // add negative contact force here since the time integrator handed me a rhs!
  if (f->Update(-1.,*fc,1.)) dserror("update went wrong");
  dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true,Add);
  kt->UnComplete();
  kt->Add(*kc,false,1.,1.);
  kt->Complete();

  return;
}


/*----------------------------------------------------------------------*
 |  read restart information for contact                     seitz 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategy::DoReadRestart(
    IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis,
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = DRT::INPUT::IntegralValue<int>(Params(),
      "RESTART_WITH_CONTACT");
  if (restartwithcontact) dserror("not supported for nitsche contact");

  // set restart displacement state
  SetState(MORTAR::state_new_displacement, *dis);
  SetState(MORTAR::state_old_displacement, *dis);

  // Evaluation for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
    interface_[i]->Initialize();

  if (friction_)
  {
    for (int i=0;i<(int)interface_.size();++i)
    {
      interface_[i]->EvaluateNodalNormals();
      interface_[i]->ExportNodalNormals();
    }
    StoreToOld(MORTAR::StrategyBase::n_old);
  }

  if (DRT::INPUT::IntegralValue<int>(Params(),"NITSCHE_PENALTY_ADAPTIVE"))
    UpdateTraceIneqEtimates();

  return;
}

void CONTACT::CoNitscheStrategy::SetState(const enum MORTAR::StateType& statename,
    const Epetra_Vector& vec)
{
  if (statename==MORTAR::state_new_displacement)
  {
    double inf_delta=0.;
    if (curr_state_==Teuchos::null)
    {
      curr_state_=Teuchos::rcp(new Epetra_Vector(vec));
      inf_delta=1.e12;
    }
    else
    {
      Epetra_Vector delta(vec);
      delta.Update(-1.,*curr_state_,1.);
      delta.NormInf(&inf_delta);
    }
    if (inf_delta<1.e-12)
      return;
    else
    {
      curr_state_eval_=false;
      (*curr_state_)=vec;
      CoAbstractStrategy::SetState(statename,vec);
      SetParentState(statename,vec);
    }
  }
  else
  {
    curr_state_eval_=false;
    CoAbstractStrategy::SetState(statename,vec);
  }
  return;
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategy::SetParentState(const enum MORTAR::StateType& statename,
    const Epetra_Vector& vec)
{
  Teuchos::RCP<DRT::Discretization> dis = DRT::Problem::Instance()->GetDis("structure");
  if (dis==Teuchos::null)
    dserror("didn't get my discretization");
  if (statename==MORTAR::state_new_displacement)
  {
    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(),true));
    LINALG::Export(vec,*global);

    //set state on interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      DRT::Discretization& idiscret_ = interface_[i]->Discret();

      for (int j=0; j < interface_[i]->Discret().ElementColMap()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->Discret().ElementColMap()->GID(j);

        MORTAR::MortarElement* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;

        //this gets values in local order
        ele->ParentElement()->LocationVector(*dis,lm,lmowner,lmstride);

        std::vector<double> myval;
        DRT::UTILS::ExtractMyValues(*global,myval,lm);

        ele->MoData().ParentDisp() = myval;
        ele->MoData().ParentDof() = lm;
      }
    }
  }
  return;
}

void CONTACT::CoNitscheStrategy::Evaluate(
    CONTACT::ParamsInterface& cparams,
    const std::vector<Teuchos::RCP<const Epetra_Vector> >* eval_vec)
{
  PreEvaluate(cparams);
  const enum MORTAR::ActionType& act = cparams.GetActionType();
  switch(act)
  {
  case MORTAR::eval_none: break;
  case MORTAR::eval_force:
  case MORTAR::eval_force_stiff:
    Integrate(cparams);
    break;
  case MORTAR::eval_run_post_evaluate: break;
  case MORTAR::eval_reset:
    SetState(MORTAR::state_new_displacement, *((*eval_vec)[0])); break;
  case MORTAR::eval_recover: break;
  case MORTAR::eval_run_pre_compute_x: break;
  case MORTAR::eval_run_post_iterate: break;
  default: dserror("unknown action"); break;
  }
  return;
}

void CONTACT::CoNitscheStrategy::Integrate(CONTACT::ParamsInterface& cparams)
{
  // we already did this displacement state
  if (curr_state_eval_==true)
    return;

  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // Evaluation for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    interface_[i]->IParams().set<double>("TIMESTEP",cparams.GetDeltaTime());
    interface_[i]->Initialize();
    interface_[i]->Evaluate(0,step_,iter_);

    //store required integration time
    inttime_ += Interfaces()[i]->Inttime();
  }

  // check the parallel distribution
  CheckParallelDistribution(t_start);

  // now we also did this state
  curr_state_eval_=true;

  // ... and we can assemble the matric and rhs
  fc_=CreateRhsBlockPtr(DRT::UTILS::block_displ);
  kc_=CreateMatrixBlockPtr(DRT::UTILS::block_displ_displ);

  return;
}
Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategy::CreateRhsBlockPtr(
     const enum DRT::UTILS::VecBlockType& bt) const
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");

  Teuchos::RCP<Epetra_FEVector> fc=Teuchos::rcp(new Epetra_FEVector(curr_state_->Map()));
  for (int i = 0; i < (int) interface_.size(); ++i)
    for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
    {
      MORTAR::MortarElement* mele =dynamic_cast<MORTAR::MortarElement*>(interface_[i]->Discret().gElement(
          interface_[i]->Discret().ElementColMap()->GID(e)));
      mele->GetNitscheContainer().AssembleRHS(mele,bt,fc);
    }
  if(fc->GlobalAssemble(Add,false)!=0) dserror("GlobalAssemble failed");

  return fc;
}

Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategy::GetRhsBlockPtr(
     const enum DRT::UTILS::VecBlockType& bt) const
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");
  if (bt==DRT::UTILS::block_displ )
    return Teuchos::rcp(new Epetra_Vector(Copy,*(fc_),0));

  return Teuchos::null;
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategy::CreateMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt)
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");

  Teuchos::RCP<LINALG::SparseMatrix> kc = Teuchos::rcp(new LINALG::SparseMatrix(
      *Teuchos::rcpFromRef<const Epetra_Map>(*DRT::Problem::Instance()->
          GetDis("structure")->DofRowMap()),100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  for (int i = 0; i < (int) interface_.size(); ++i)
    for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
    {
      MORTAR::MortarElement* mele =dynamic_cast<MORTAR::MortarElement*>(interface_[i]->Discret().gElement(
          interface_[i]->Discret().ElementColMap()->GID(e)));
      mele->GetNitscheContainer().AssembleMatrix(mele,bt,kc);
    }
  if(dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true,Add))
    dserror("GlobalAssemble(...) failed");

  return kc;
}

Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt) const
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");

  if (bt==DRT::UTILS::block_displ_displ)
    return kc_;

  return Teuchos::null;
}

void CONTACT::CoNitscheStrategy::Setup(bool redistributed, bool init)
{
  if (isselfcontact_)
    dserror("no self contact with Nitsche yet");
  ReconnectParentElements();
  curr_state_=Teuchos::null;
  curr_state_eval_=false;
}

void CONTACT::CoNitscheStrategy::UpdateTraceIneqEtimates()
{
  INPAR::CONTACT::NitscheWeighting NitWgt =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(Params(),"NITSCHE_WEIGHTING");
  for (int i=0;i<(int)interface_.size();++i)
    for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
    {
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(
          interface_[i]->Discret().gElement(
              interface_[i]->Discret().ElementColMap()->GID(e)));
      if (NitWgt==INPAR::CONTACT::NitWgt_slave && mele->IsSlave()==false)
        continue;
      if (NitWgt==INPAR::CONTACT::NitWgt_master && mele->IsSlave()==true)
        continue;
      mele->EstimateNitscheTraceMaxEigenvalueCombined();
    }
  return;
}

void CONTACT::CoNitscheStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  if (DRT::INPUT::IntegralValue<int>(Params(),"NITSCHE_PENALTY_ADAPTIVE"))
    UpdateTraceIneqEtimates();
  if (friction_)
  {
    StoreToOld(MORTAR::StrategyBase::n_old);
    SetState(MORTAR::state_old_displacement,*dis);
  }

  return;
}

void CONTACT::CoNitscheStrategy::EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
{
  SetState(MORTAR::state_new_displacement,*dis);

  if (friction_)
  {
    for (int i=0;i<(int)interface_.size();++i)
    {
      interface_[i]->EvaluateNodalNormals();
      interface_[i]->ExportNodalNormals();
    }
    StoreToOld(MORTAR::StrategyBase::n_old);
  }

  UpdateTraceIneqEtimates();
  return;
}


/*----------------------------------------------------------------------------------------------*
 |  Reconnect Contact Element -- Parent Element Pointers (required for restart)       ager 04/16|
 *---------------------------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategy::ReconnectParentElements()
{
  Teuchos::RCP<DRT::Discretization> voldis=DRT::Problem::Instance()->GetDis("structure");

  for (int intidx = 0; intidx < (int) ContactInterfaces().size(); ++intidx)
  {
    const Epetra_Map* elecolmap = voldis->ElementColMap();

    const Epetra_Map* ielecolmap = ContactInterfaces()[intidx]->Discret().ElementColMap();

    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = ContactInterfaces()[intidx]->Discret().gElement(gid);
      if (!ele)
        dserror("ERROR: Cannot find element with gid %", gid);
      DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);

      int volgid = faceele->ParentElementId();
      if (elecolmap->LID(volgid) == -1) //Volume Discretization has not Element
        dserror("CoManager::ReconnectParentElements: Element %d does not exist on this Proc!",volgid);

      DRT::Element* vele = voldis->gElement(volgid);
      if (!vele)
        dserror("ERROR: Cannot find element with gid %", volgid);

      faceele->SetParentMasterElement(vele,faceele->FaceParentNumber());

      DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>* vele_plast=
          dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>*>(vele);
      if (vele_plast)
        vele_plast->SetIsNitscheContactEle(true);

    }
  }
}
