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
      mele->GetNitscheContainer().Assemble(mele,fc,kc);
      mele->GetNitscheContainer().Clear();
    }
  }
  if(fc->GlobalAssemble(Add,false)!=0) dserror("GlobalAssemble failed");
  if (f->Update(1.,*fc,1.)) dserror("update went wrong");
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

  return;
}

/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                         seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategy::SetParentState(const std::string& statename,
    const Teuchos::RCP<Epetra_Vector> vec,
    const Teuchos::RCP<DRT::Discretization> dis)
{
  if (statename == "displacement")
  {
    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(),true));
    LINALG::Export(*vec,*global);

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

