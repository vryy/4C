/*!----------------------------------------------------------------------
\file adapter_coupling_ehl_mortar.cpp
\brief mortar coupling terms of ehl

\level 3
\maintainer Alexander Seitz

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "adapter_coupling_ehl_mortar.H"

#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

ADAPTER::CouplingEhlMortar::CouplingEhlMortar() : CouplingNonLinMortar()
{
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(DRT::Problem::Instance()->MortarCouplingParams(),"PARALLEL_REDIST")
      != INPAR::MORTAR::parredist_none)
    dserror("EHL does not support parallel redistribution. Set \"PARALLEL_REDIST none\" in section \"MORTAR COUPLING\"");
}

/*----------------------------------------------------------------------*
 |  read mortar condition                                               |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::ReadMortarCondition(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<DRT::Discretization>   slavedis,
    std::vector<int>                    coupleddof,
    const std::string&                  couplingcond,
    Teuchos::ParameterList&             input,
    std::map<int, DRT::Node*>& mastergnodes,
    std::map<int, DRT::Node*>& slavegnodes,
    std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
    std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements
    )
{
  ADAPTER::CouplingNonLinMortar::ReadMortarCondition(
      masterdis,slavedis,coupleddof,couplingcond,input,
      mastergnodes,slavegnodes,masterelements,slaveelements);

  input.set<int>("PROBTYPE", INPAR::CONTACT::ehl);
}
/*----------------------------------------------------------------------*
 |  perform interface integration and assembly                          |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::Integrate(
    Teuchos::RCP<const Epetra_Vector> disp,
    const double dt)
{
  // safety check
  CheckSetup();

  // return if this state has already been evaluated
  if (AlreadyEvaluated(disp))
    return;

  // set current displ state
  interface_->SetState(MORTAR::state_new_displacement,*disp);

  // init internal data
  interface_->Initialize();
  interface_->SetElementAreas();
  // call interface evaluate (d,m,gap...)
  interface_->Evaluate();

  // some first assemblies, that don't require any additional states
  D_   = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));
  M_   = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));
  interface_->AssembleDM(*D_,*M_);
  D_->Complete();
  M_->Complete(*masterdofrowmap_,*slavedofrowmap_);
  N_->Complete(*smdofrowmap_,*slavedofrowmap_);
  AssembleRealGap();
  AssembleRealGapDeriv();
  AssembleNormals();
  AssembleNormalsDeriv();
  AssembleSurfGrad();
  AssembleInterfaceVelocities(dt);

  // save that state as the last evaluated one
  evaluated_state_=Teuchos::rcp(new Epetra_Vector(*disp));

  // all done
  return;
}

bool ADAPTER::CouplingEhlMortar::AlreadyEvaluated(
    Teuchos::RCP<const Epetra_Vector> disp)
{
  if (evaluated_state_.is_null())
    return false;
  Teuchos::RCP<Epetra_Vector> diff=Teuchos::rcp(new Epetra_Vector(*disp));
  if(diff->Update(-1.,*evaluated_state_,1.))dserror("update failed");
  double inf_diff=-1.;
  if (diff->NormInf(&inf_diff)) dserror("NormInf failed");
  if (inf_diff<1.e-13)
    return true;

  return false;
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinD(const Teuchos::RCP<Epetra_Vector> x //slave dof vector
)
{
  Teuchos::RCP<LINALG::SparseMatrix> DLinEHL = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,true,false,LINALG::SparseMatrix::FE_MATRIX));
  DLinEHL->Zero();
  DLinEHL->UnComplete();

  interface_->AssembleCoupLinD(*DLinEHL,x);

  DLinEHL->Complete(*smdofrowmap_, *slavedofrowmap_);

  return DLinEHL;
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinM(const Teuchos::RCP<Epetra_Vector> x //slave dof vector
    )
{
  Teuchos::RCP<LINALG::SparseMatrix> MLinEHL =  Teuchos::rcp(new LINALG::SparseMatrix(*masterdofrowmap_,81,true,false,LINALG::SparseMatrix::FE_MATRIX));
  MLinEHL->Zero();
  MLinEHL->UnComplete();

  interface_->AssembleCoupLinM(*MLinEHL,x);

  MLinEHL->Complete(*smdofrowmap_, *masterdofrowmap_);

  return MLinEHL;
}

void ADAPTER::CouplingEhlMortar::AssembleNormals()
{
  normals_=Teuchos::rcp(new Epetra_Vector(*SlaveDofMap(),true));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    for (int d=0;d<interface_->Dim();++d)
      normals_->ReplaceGlobalValue(cnode->Dofs()[d],0,cnode->MoData().n()[d]);
  }
}


void ADAPTER::CouplingEhlMortar::AssembleNormalsDeriv()
{
  Nderiv_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));
  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    for (int d=0;d<Interface()->Dim();++d)
      for(auto p=cnode->CoData().GetDerivN()[d].begin();p!=cnode->CoData().GetDerivN()[d].end();++p)
        Nderiv_->Assemble(p->second,cnode->Dofs()[d],p->first);
  }
  Nderiv_->Complete();
}

void ADAPTER::CouplingEhlMortar::AssembleRealGap()
{
  nodal_gap_ = Teuchos::rcp(new Epetra_Vector(*slavenoderowmap_,true));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");
    double real_gap = cnode->CoData().Getg();
    switch (cnode->MoData().GetD().size())
    {
    case 0: break;
    case 1:
      if (cnode->MoData().GetD().begin()->first!=cnode->Id())
        dserror("something is wrong. Here should by my own Id");
      real_gap/=cnode->MoData().GetD().at(cnode->Id());
      break;
    default: dserror("GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using duals?");
    }
    nodal_gap_->ReplaceGlobalValue(cnode->Id(),0,real_gap);
  }
}

void ADAPTER::CouplingEhlMortar::AssembleRealGapDeriv()
{
  deriv_nodal_gap_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    if (cnode->CoData().GetDerivD().size()!=cnode->MoData().GetD().size())
      dserror("size inconsistency");

    const double w_gap = cnode->CoData().Getg();
    double d=-1.;
    switch (cnode->CoData().GetDerivD().size())
    {
    case 0: break;
    case 1:
      if (cnode->CoData().GetDerivD().begin()->first!=cnode->Id())
        dserror("something is wrong. Here should by my own Id");
      d=cnode->MoData().GetD().at(cnode->Id());
      break;
    default: dserror("GetDerivD should be of size 0 (unprojectable) or 1 (projectable). Are you not using duals?");
    }

    if (cnode->CoData().GetDerivD().size())
      for (auto p=cnode->CoData().GetDerivD().at(cnode->Id()).begin();p!=cnode->CoData().GetDerivD().at(cnode->Id()).end();++p)
      {
        const double val = -w_gap/(d*d)*p->second;
        for (int d=0;d<interface_->Dim();++d)
          deriv_nodal_gap_->Assemble(val,cnode->Dofs()[d],p->first);
      }

    if (d==-1 && cnode->CoData().GetDerivG().size()!=0)
      dserror("inconsistency");

    if (cnode->CoData().GetDerivG().size())
      for (auto p=cnode->CoData().GetDerivG().begin();p!=cnode->CoData().GetDerivG().end();++p)
      {
        const double val = p->second/d;
        for (int d=0;d<interface_->Dim();++d)
          deriv_nodal_gap_->Assemble(val,cnode->Dofs()[d],p->first);
      }
  }
  deriv_nodal_gap_->Complete(*smdofrowmap_,*slavedofrowmap_);
}

void ADAPTER::CouplingEhlMortar::AssembleInterfaceVelocities(const double dt)
{
  relTangVel_=Teuchos::rcp(new Epetra_Vector(*slavedofrowmap_));
  avTangVel_ =Teuchos::rcp(new Epetra_Vector(*slavedofrowmap_));
  relTangVel_deriv_= Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));
  avTangVel_deriv_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");


    double d_val =0.;
    switch (cnode->MoData().GetD().size())
    {
    case 0: break;
    case 1:
      if (cnode->MoData().GetD().begin()->first!=cnode->Id())
        dserror("something is wrong. Here should by my own Id");
      d_val=cnode->MoData().GetD().at(cnode->Id());
      break;
    default: dserror("GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using duals?");
    }

    if (d_val==0.)
      continue;

    for (int d=0;d<Interface()->Dim();++d)
    {
      relTangVel_->ReplaceGlobalValue(cnode->Dofs()[d],0,cnode->CoEhlData().GetWeightedRelTangVel()(d)/d_val);
      avTangVel_ ->ReplaceGlobalValue(cnode->Dofs()[d],0,cnode->CoEhlData().GetWeightedAvTangVel ()(d)/d_val);
    }

    for (auto p=cnode->CoData().GetDerivD().at(cnode->Id()).begin();p!=cnode->CoData().GetDerivD().at(cnode->Id()).end();++p)
    {
      const int col=p->first;
      for (int d=0;d<Interface()->Dim();++d)
      {
        const int row=cnode->Dofs()[d];
        const double rel_val = -cnode->CoEhlData().GetWeightedRelTangVel()(d)/(d_val*d_val)*p->second;
        const double  av_val = -cnode->CoEhlData().GetWeightedAvTangVel ()(d)/(d_val*d_val)*p->second;
        relTangVel_deriv_->Assemble(rel_val,row,col);
        avTangVel_deriv_ ->Assemble( av_val,row,col);
      }
    }
    for (auto p=cnode->CoEhlData().GetWeightedAvTangVelDeriv().begin();p!=cnode->CoEhlData().GetWeightedAvTangVelDeriv().end();++p)
    {
      const int col=p->first;
      for (int d=0;d<Interface()->Dim();++d)
      {
        const int row=cnode->Dofs()[d];
        const double val = p->second(d)/d_val;
        avTangVel_deriv_ ->Assemble(val,row,col);
      }
    }
    for (auto p=cnode->CoEhlData().GetWeightedRelTangVelDeriv().begin();p!=cnode->CoEhlData().GetWeightedRelTangVelDeriv().end();++p)
    {
      const int col=p->first;
      for (int d=0;d<Interface()->Dim();++d)
      {
        const int row=cnode->Dofs()[d];
        const double val = p->second(d)/d_val;
        relTangVel_deriv_ ->Assemble(val,row,col);
      }
    }
  }

  relTangVel_->Scale(1./dt);
  avTangVel_->Scale(1./dt);
  relTangVel_deriv_->Complete(*smdofrowmap_,*slavedofrowmap_);
  avTangVel_deriv_ ->Complete(*smdofrowmap_,*slavedofrowmap_);
  relTangVel_deriv_->Scale(1./dt);
  avTangVel_deriv_->Scale(1./dt);


}

void ADAPTER::CouplingEhlMortar::AssembleSurfGrad()
{
  SurfGrad_= Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false,LINALG::SparseMatrix::FE_MATRIX));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("ERROR: Cannot find node");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode)
      dserror("this is not a contact node");

    double dval=1.;
    switch (cnode->MoData().GetD().size())
    {
    case 0: dval=1.e32; break; // large number so no tangential gradient
    case 1:
      if (cnode->MoData().GetD().begin()->first!=cnode->Id())
        dserror("something is wrong. Here should by my own Id");
      dval=cnode->MoData().GetD().at(cnode->Id());
      break;
    default: dserror("GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using duals?");
    }

    for (auto p=cnode->CoEhlData().GetSurfGrad().begin();p!=cnode->CoEhlData().GetSurfGrad().end();++p)
      for (int d=0;d<Interface()->Dim();++d)
        SurfGrad_->Assemble(p->second(d)/dval,cnode->Dofs()[d],p->first);
  }

  SurfGrad_->Complete();
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleSurfGradDeriv(const Teuchos::RCP<const Epetra_Vector> x)
{
  Teuchos::RCP<LINALG::SparseMatrix> SurfGradDeriv = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_,81,false,false,LINALG::SparseMatrix::FE_MATRIX));

  for (int i=0;i<interface_->SlaveRowNodes()->NumMyElements();++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("ERROR: Cannot find node");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode)
      dserror("this is not a contact node");

    double dval=1.;
    switch (cnode->MoData().GetD().size())
    {
    case 0: dval=1.e32; break; // large number so no tangential gradient
    case 1:
      if (cnode->MoData().GetD().begin()->first!=cnode->Id())
        dserror("something is wrong. Here should by my own Id");
      dval=cnode->MoData().GetD().at(cnode->Id());
      break;
    default: dserror("GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using duals?");
    }

    for(auto p=cnode->CoEhlData().GetSurfGradDeriv().begin();p!=cnode->CoEhlData().GetSurfGradDeriv().end();++p)
    {
      const int col = p->first;
      for (auto q=p->second.begin();q!=p->second.end();++q)
      {
        const int lid = x->Map().LID(q->first);
        if (lid<0) dserror("not my gid");
        const double x_val = x->operator [](lid);
        for (int d=0;d<Interface()->Dim();++d)
        {
          const double val = x_val*q->second(d)/dval;
          SurfGradDeriv->Assemble(val,cnode->Dofs()[d],col);
        }
      }
    }

    if (cnode->CoData().GetDerivD().size())
      for (auto p=cnode->CoData().GetDerivD().at(cnode->Id()).begin();p!=cnode->CoData().GetDerivD().at(cnode->Id()).end();++p)
      {
        const int col=p->first;

        for (auto q=cnode->CoEhlData().GetSurfGrad().begin();q!=cnode->CoEhlData().GetSurfGrad().end();++q)
          for (int d=0;d<Interface()->Dim();++d)
          {
            const int row =cnode->Dofs()[d];
            const int x_gid = q->first;
            const int x_lid = x->Map().LID(x_gid);
            if (x_lid<0)
              dserror("not my gid");
            double x_val = x->operator [](x_lid);
            const double val = -x_val*q->second(d)/(dval*dval)*p->second;
            SurfGradDeriv->Assemble(val,row,col);
          }
      }
  }
  SurfGradDeriv->Complete(*smdofrowmap_,*slavedofrowmap_);
  return SurfGradDeriv;
}
