/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_strategy_tsi.cpp

\brief Nitsche contact solving strategy

\level 3

\maintainer Alexander Seitz

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



void CONTACT::CoNitscheStrategyTsi::SetState(const enum MORTAR::StateType& statename,
    const Epetra_Vector& vec)
{
  if (statename==MORTAR::state_temperature)
  {
    double inf_delta=0.;
    if (curr_state_temp_==Teuchos::null)
    {
      curr_state_temp_=Teuchos::rcp(new Epetra_Vector(vec));
      inf_delta=1.e12;
    }
    else
    {
      Epetra_Vector delta(vec);
      delta.Update(-1.,*curr_state_temp_,1.);
      delta.NormInf(&inf_delta);
    }
    if (inf_delta<1.e-16)
      return;
    else
    {
      curr_state_eval_=false;
      (*curr_state_temp_)=vec;
      SetParentState(statename,vec);
    }
  }
  else
    CONTACT::CoNitscheStrategy::SetState(statename,vec);

  return;
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategyTsi::SetParentState(const enum MORTAR::StateType& statename,
    const Epetra_Vector& vec)
{
  Teuchos::RCP<DRT::Discretization> dis = DRT::Problem::Instance()->GetDis("structure");
  if (dis==Teuchos::null)
    dserror("didn't get my discretization");

  if (statename==MORTAR::state_temperature)
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

        DRT::Element* e = idiscret_.gElement(gid);
        if (e==NULL)
          dserror("basic element not found");

        MORTAR::MortarElement* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));

        if (ele==NULL)
          dserror("element not found");

        if (ele->ParentElement()==NULL)
          dserror("parent element not set");

        if (ele->ParentElement()->NumNode()==0)
          dserror("parent element nodes not set");

        // since contact can't handle multiple dofsets, temperatures are associated with first dof of each node
        std::vector<int> lm(ele->ParentElement()->NumNode());
        for (int n=0;n<ele->ParentElement()->NumNode();++n)
          lm.at(n)=dis->Dof(0,ele->ParentElement()->Nodes()[n],0);

        std::vector<double> myval;
        DRT::UTILS::ExtractMyValues(*global,myval,lm);

        ele->MoData().ParentTemp() = myval;
      }
    }
  }
  else
    CONTACT::CoNitscheStrategy::SetParentState(statename,vec);

  return;
}

void CONTACT::CoNitscheStrategyTsi::Setup(bool redistributed, bool init)
{
  CONTACT::CoNitscheStrategy::Setup(redistributed,init);

  curr_state_temp_=Teuchos::null;
  return;
}

void CONTACT::CoNitscheStrategyTsi::UpdateTraceIneqEtimates()
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

  // add estimate for thermo-penalty
//  dserror("not yet implemented");
  return;
}


void CONTACT::CoNitscheStrategyTsi::SetAlphafThermo(const Teuchos::ParameterList& tdyn)
{
  INPAR::THR::DynamicType dyn_type = DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn,"DYNAMICTYP");
  switch (dyn_type)
  {
  case INPAR::THR::dyna_genalpha:
    thermo_alpha_ = tdyn.sublist("GENALPHA").get<double>("ALPHA_F");
    break;
  case INPAR::THR::dyna_onesteptheta:
    thermo_alpha_ = tdyn.sublist("ONESTEPTHETA").get<double>("THETA");
    break;
  case INPAR::THR::dyna_statics:
    thermo_alpha_ = 1.;
    break;
  default:
    dserror("unknown thermal time integration type");
  }
  return;
}
Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategyTsi::GetRhsBlockPtr(
     const enum DRT::UTILS::VecBlockType& bt) const
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");

  switch(bt)
  {
  case DRT::UTILS::block_temp: return Teuchos::rcp(new Epetra_Vector(Copy,*(ft_),0));
  default: return CONTACT::CoNitscheStrategy::GetRhsBlockPtr(bt);
  }
}


Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoNitscheStrategyTsi::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt,
    const CONTACT::ParamsInterface* cparams) const
{
  if (curr_state_eval_==false)
    dserror("you didn't evaluate this contact state first");

  switch(bt)
  {
  case DRT::UTILS::block_temp_temp:   return ktt_; break;
  case DRT::UTILS::block_temp_displ:  return ktd_; break;
  case DRT::UTILS::block_displ_temp:  return kdt_; break;
  default: return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bt,cparams); break;
  }
}


void CONTACT::CoNitscheStrategyTsi::Integrate(CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  if (coupST_==Teuchos::null)
  {
#ifdef DEBUG
    std::cout<<"\nWARNING: we are skipping the assembly of TSI coupling matrices in CONTACT::CoNitscheStrategyTsi::Integrate\n"
        "We do this, since the coupling object is not set. This may happen during constructor phases.\n"
        "If this warning appears long after the construction, you're in big trouble.\n" << std::endl;
#endif
    return;
  }

  Teuchos::RCP<Epetra_FEVector> ft_s = CreateRhsBlockPtr(DRT::UTILS::block_temp);
  Teuchos::RCP<Epetra_FEVector> ft_s_exp = Teuchos::rcp(new Epetra_FEVector(*(coupST_->MasterDofMap())));
  LINALG::Export(*ft_s,*ft_s_exp);
  ft_=coupST_->MasterToSlave(ft_s_exp);

  Teuchos::RCP<LINALG::SparseMatrix> tmp=CreateMatrixBlockPtr(DRT::UTILS::block_temp_temp);
  ktt_=Teuchos::rcp(new LINALG::SparseMatrix(*coupST_->SlaveDofMap(),81,true,false));
  FSI::UTILS::MatrixRowColTransform()(
      *tmp,1.,
      ADAPTER::CouplingMasterConverter(*coupST_),
      ADAPTER::CouplingMasterConverter(*coupST_),
      *ktt_,true,false);
  ktt_->Complete();

  tmp=CreateMatrixBlockPtr(DRT::UTILS::block_temp_displ);
  ktd_=Teuchos::rcp(new LINALG::SparseMatrix(*coupST_->SlaveDofMap(),81,true,false));
  FSI::UTILS::MatrixRowTransform()(
      *tmp,1.,
      ADAPTER::CouplingMasterConverter(*coupST_),
      *ktd_,true);
  ktd_->Complete(*DRT::Problem::Instance()->GetDis("structure")->DofColMap(),*coupST_->SlaveDofMap());

  tmp=CreateMatrixBlockPtr(DRT::UTILS::block_displ_temp);
  kdt_=Teuchos::rcp(new LINALG::SparseMatrix(tmp->RowMap(),81,true,false));
  FSI::UTILS::MatrixColTransform()(
      tmp->RowMap(),
      tmp->ColMap(),
      *tmp,1.,
      ADAPTER::CouplingMasterConverter(*coupST_),
      *kdt_,true,false);
  kdt_->Complete(*coupST_->SlaveDofMap(),*DRT::Problem::Instance()->GetDis("structure")->DofRowMap());

  return;
}
