/*!----------------------------------------------------------------------
\file drt_contact_penalty_strategy.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_penalty_strategy.H"
#include "drt_contact_abstract_strategy.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::PenaltyStrategy::PenaltyStrategy(RCP<Epetra_Map> problemrowmap,
                                          Teuchos::ParameterList params,
                                          vector<RCP<CONTACT::Interface> > interface, int dim,
                                          RCP<Epetra_Comm> comm, double alphaf) :
AbstractStrategy(problemrowmap, params, interface, dim, comm, alphaf)
{
  // check if prerequisites for penalty strategies are met
#ifdef CONTACTBOUNDMOD
  dserror("Contact Boundary Modification not implemented for Penalty Methods.");
#endif

  // initialize constraint norm and initial penalty
  constrnorm_ = 0.0;
  constrnormtan_ = 0.0;
  initialpenalty_ = Params().get<double>("PENALTYPARAM");
  initialpenaltytan_ = Params().get<double>("PENALTYPARAMTAN");
}

/*----------------------------------------------------------------------*
 |  save the gap-scaling kappa from reference config          popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::SaveReferenceState(const RCP<Epetra_Vector> dis)
{
  // initialize the displacement field
  SetState("displacement", dis);

  // kappa will be the shape function integral on the slave sides
  // (1) build the nodal information
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // interface needs to be complete
    if (!interface_[i]->Filled() && Comm().MyPID()==0)
      dserror("ERROR: FillComplete() not called on interface %", i);

    // get out of here if not participating in interface
    if (!interface_[i]->lComm())
      continue;

    // do the computation of nodal shape function integral
    // (for convenience, the results will be stored in nodal gap)

    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int j=0; j<interface_[i]->SlaveColElements()->NumMyElements(); ++j)
    {
      int gid1 = interface_[i]->SlaveColElements()->GID(j);
      DRT::Element* ele1 = interface_[i]->Discret().gElement(gid1);
      if (!ele1)
        dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);

      interface_[i]->IntegrateKappaPenalty(*selement);
    }

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);

      // get nodal weighted gap
      // (this is where we stored the shape function integrals)
      double gap = cnode->Getg();

      // store kappa as the inverse of gap
      // (this removes the scaling introduced by weighting the gap!!!)
      cnode->Kappa() = 1.0/gap;

      //cout << "S-NODE #" << gid << " kappa=" << cnode->Kappa() << endl;
    }
  }
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::Initialize()
{
  // (re)setup global matrices containing lagrange multiplier derivatives
  linzmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));

  z_ = LINALG::CreateVector(*gsdofrowmap_, true);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate contact (incl. friction) and create linear system popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::EvaluateContact(RCP<LINALG::SparseMatrix> kteff,
                                               RCP<Epetra_Vector> feff)
{

  // uncomment this if you want to do finite difference checks
  // when enabled, kc1 and kc2 as part of kteff and fc as part of
  // feff is stored separately so they can be checked on their own

//#define FDCHECKS

  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...

	// in the beginning of this function, the regularized contact forces
	// in normal and tangential direction are evaluated from geometric
	// measures (gap and relative tangential velocity). Here, also active and
	// slip nodes are detected. Then, the insertion of the according stiffness
	// blocks takes place.

  bool isincontact = false;
  bool activesetchange = false;

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    bool localisincontact = false;
    bool localactivesetchange = false;

    // evaluate lagrange multipliers (regularized forces) in normal direction
    // and nodal derivz matrix values, store them in nodes
    interface_[i]->AssembleRegNormalForces(localisincontact, localactivesetchange);

    // evaluate lagrange multipliers (regularized forces) in tangential direction
    INPAR::CONTACT::ContactFrictionType ftype =
      Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");
    INPAR::CONTACT::SolvingStrategy soltype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");

    if((ftype==INPAR::CONTACT::friction_coulomb or ftype==INPAR::CONTACT::friction_stick)
    	  and soltype==INPAR::CONTACT::solution_penalty)
      interface_[i]->AssembleRegTangentForcesPenalty();

    if((ftype==INPAR::CONTACT::friction_coulomb or ftype==INPAR::CONTACT::friction_stick)
    	  and soltype==INPAR::CONTACT::solution_auglag)
      interface_[i]->AssembleRegTangentForcesAugmented();

    isincontact = isincontact || localisincontact;
    activesetchange = activesetchange || localactivesetchange;
  }

  // broadcast contact status & active set change
  int globalcontact, globalchange = 0;
  int localcontact = isincontact;
  int localchange = activesetchange;

  Comm().SumAll(&localcontact, &globalcontact, 1);
  Comm().SumAll(&localchange, &globalchange, 1);

  if (globalcontact>=1)
    isincontact_=true;
  else
    isincontact_ = false;

  if( (Comm().MyPID()==0) && (globalchange>=1) )
    cout << "ACTIVE SET HAS CHANGED..." << endl;

  // (re)setup active global Epetra_Maps
  // the map of global active nodes is needed for the penalty case, too.
  // this is due to the fact that we want to monitor the constraint norm
  // of the active nodes
  gactivenodes_ = null;
  gslipnodes_ = null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
  }

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if( ! IsInContact() )
    return;

  // since we will modify kteff by adding additional stiffness entries,
  // we have to uncomplete it
  // TODO move the call of EvaluateContact in ContactStruGenAlpha before the completion of kteff to avoid this
  kteff->UnComplete();

  // we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_, gmdofrowmap_, false);

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble global lagrangian multiplier vector
    interface_[i]->AssembleLM(*z_);
    // assemble global derivatives of lagrangian multipliers
    interface_[i]->AssembleLinZ(*linzmatrix_);
    // assemble global derivatives of mortar D and M matrices
    interface_[i]->AssembleLinDM(*lindmatrix_, *linmmatrix_);
  }

  // FillComplete() global matrices LinD, LinM, LinZ
  lindmatrix_->Complete(*gsmdofs, *gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofs, *gmdofrowmap_);
  linzmatrix_->Complete(*gsmdofs, *gsdofrowmap_);

#ifdef CONTACTFDPENALTYDERIVTRAC
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

  // check derivatives of penalty traction
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if( IsInContact() )
    {
    	if(ftype==INPAR::CONTACT::friction_coulomb )
      {
    		cout << "LINZMATRIX" << *linzmatrix_ << endl;
    	  interface_[i]->FDCheckPenaltyTracFric();
      }
      else if (ftype==INPAR::CONTACT::friction_none)
      {
        cout << "-- CONTACTFDDERIVZ --------------------" << endl;
        interface_[i]->FDCheckPenaltyTracNor();
        cout << "-- CONTACTFDDERIVZ --------------------" << endl;
      }
      else
      	dserror("Error: FD Check for this friction type not implemented!");
    }
  }
#endif

#ifdef FDCHECKS
  // for debugging purposes kc is stored separately
  RCP<LINALG::SparseMatrix> kc1 = rcp(new LINALG::SparseMatrix(*gsmdofs,81,true,true));
  RCP<LINALG::SparseMatrix> kc2 = rcp(new LINALG::SparseMatrix(*gsmdofs,81,true,true));
#endif

  // **********************************************************************
  // Build Contact Stiffness #1
  // **********************************************************************
  // involving contributions of derivatives of D and M:
  //  Kc,1 = delta[ 0 -M(transpose) D] * LM

  // these calculations are readily done and incorporated into linD and linM !
#ifndef FDCHECKS
  kteff->Add(*lindmatrix_, false, 1.0-alphaf_, 1.0);
  kteff->Add(*linmmatrix_, false, 1.0-alphaf_, 1.0);
#else
  kc1->Add(*lindmatrix_, false, 1.0-alphaf_, 1.0);
  kc1->Add(*linmmatrix_, false, 1.0-alphaf_, 1.0);
#endif

  // **********************************************************************
  // Build Contact Stiffness #2
  // **********************************************************************
  // involving contributions of derivatives of lagrange multipliers:
  //  Kc,2= [ 0 -M(transpose) D] * deltaLM

  // Multiply Mortar matrices D and M with LinZ
  RCP<LINALG::SparseMatrix> dtilde = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  RCP<LINALG::SparseMatrix> mtilde = rcp(new LINALG::SparseMatrix(*gmdofrowmap_));

  // do the multiplication dtile = D * LinZ
  dtilde = LINALG::Multiply(*dmatrix_, false, *linzmatrix_, false);
  // do the multiplication mtilde = -M(transpose) * LinZ
  mtilde = LINALG::Multiply(*mmatrix_, true, *linzmatrix_, false);

  // add to kteff
#ifndef FDCHECKS
  kteff->Add(*dtilde, false, 1.0-alphaf_, 1.0);
  kteff->Add(*mtilde, false, -(1.0-alphaf_), 1.0);
#else
  kc2->Add(*dtilde, false, 1.0-alphaf_, 1.0);
  kc2->Add(*mtilde, false, -(1.0-alphaf_), 1.0);
#endif

  // finalize the stiffness matrix
#ifdef FDCHECKS
  kc1->Complete();
  kc2->Complete();
  kteff->Add(*kc1, false, 1.0, 1.0);
  kteff->Add(*kc2, false, 1.0, 1.0);
#endif

  kteff->Complete();

  // **********************************************************************
  // Build RHS
  // **********************************************************************
  // feff += -alphaf * fc,n - (1-alphaf) * fc,n+1,k

#ifndef FDCHECKS
  {
    // we initialize fcmdold with dold-rowmap instead of gsdofrowmap
    // (this way, possible self contact is automatically included)

    RCP<Epetra_Vector> fcmdold = rcp(new Epetra_Vector(dold_->RowMap()));
    dold_->Multiply(false, *zold_, *fcmdold);
    RCP<Epetra_Vector> fcmdoldtemp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fcmdold, *fcmdoldtemp);
    feff->Update(-alphaf_, *fcmdoldtemp, 1.0);
  }

  {
    // we initialize fcmmold with mold-domainmap instead of gmdofrowmap
    // (this way, possible self contact is automatically included)

    RCP<Epetra_Vector> fcmmold = rcp(new Epetra_Vector(mold_->DomainMap()));
    mold_->Multiply(true, *zold_, *fcmmold);
    RCP<Epetra_Vector> fcmmoldtemp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fcmmold, *fcmmoldtemp);
    feff->Update(alphaf_, *fcmmoldtemp, 1.0);
  }

  {
    RCP<Epetra_Vector> fcmd = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(false, *z_, *fcmd);
    RCP<Epetra_Vector> fcmdtemp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fcmd, *fcmdtemp);
    feff->Update(-(1-alphaf_), *fcmdtemp, 1.0);
  }

  {
    RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
    mmatrix_->Multiply(true, *z_, *fcmm);
    RCP<Epetra_Vector> fcmmtemp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fcmm, *fcmmtemp);
    feff->Update(1-alphaf_, *fcmmtemp, 1.0);
  }

#else
  // for debugging purposes fc is stored separately
  RCP<Epetra_Vector> fc = rcp( new Epetra_Vector(*gsmdofs));

  {
    RCP<Epetra_Vector> fcmdold = rcp(new Epetra_Vector(*gsdofrowmap_));
    dold_->Multiply(false, *zold_, *fcmdold);
    RCP<Epetra_Vector> fcmdoldtemp = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*fcmdold, *fcmdoldtemp);
    fc->Update(-alphaf_, *fcmdoldtemp, 1.0);
  }

  {
    RCP<Epetra_Vector> fcmmold = LINALG::CreateVector(*gmdofrowmap_, true);
    mold_->Multiply(true, *zold_, *fcmmold);
    RCP<Epetra_Vector> fcmmoldtemp = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*fcmmold, *fcmmoldtemp);
    fc->Update(alphaf_, *fcmmoldtemp, 1.0);
  }

  {
    RCP<Epetra_Vector> fcmd = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(false, *z_, *fcmd);
    RCP<Epetra_Vector> fcmdtemp = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*fcmd, *fcmdtemp);
    fc->Update(-(1-alphaf_), *fcmdtemp, 1.0);
  }

  {
    RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
    mmatrix_->Multiply(true, *z_, *fcmm);
    RCP<Epetra_Vector> fcmmtemp = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*fcmm, *fcmmtemp);
    fc->Update(1-alphaf_, *fcmmtemp, 1.0);
  }

  RCP<Epetra_Vector> fctemp = LINALG::CreateVector(*problemrowmap_,true);
  LINALG::Export(*fc, *fctemp);
  feff->Update(1.0, *fctemp, 1.0);
#endif

#ifdef CONTACTFDGAP
   // FD check of weighted gap g derivatives (non-penetr. condition)

  cout << "-- CONTACTFDGAP -----------------------------" << endl;
  //cout << *smatrix_ << endl;
  interface_[0]->FDCheckGapDeriv();
  cout << "-- CONTACTFDGAP -----------------------------" << endl;

#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDPENALTYKTEFF
  // check assembled contact tangent stiffness matrix

  RCP<LINALG::SparseMatrix> kc = rcp(new LINALG::SparseMatrix(*gsmdofs,81,true,true));
  kc->Add(*kc1, false, 1.0, 1.0);
  kc->Add(*kc2, false, 1.0, 1.0);
  kc->Complete();

  cout << "-- CONTACTFDPENALTYKTEFF --------------------" << endl;
  //cout << "KC:" << endl;
  //cout << *kc  << endl;
  interface_[0]->FDCheckPenaltyKTeff(kc, fc, dold_, mold_, zold_, gsmdofs, alphaf_);
  cout << "-- CONTACTFDPENALTYKTEFF --------------------" << endl;

#endif

#ifdef CONTACTFDPENALTYKC1
  // check assembled contact tangent stiffness matrix, 1st part (LinD, LinM)

  bool lind = true;
  bool linm = false;

  cout << "-- CONTACTFDPENALTYKC1 ----------------------" << endl;

  cout << "zvector: " << *z_ << endl;
  cout << "linzmatrix: " << *linzmatrix_ << endl;

  if( lind )
  {
    cout << "lindmatrix: " << *lindmatrix_ << endl;
    RCP<Epetra_Vector> reffcd = LINALG::CreateVector(*gsdofrowmap_,true);
    dmatrix_->Multiply(false,*z_,*reffcd);
    cout << "reffcd: " << *reffcd << endl;
    interface_[0]->FDCheckPenaltyLinD(lindmatrix_, reffcd);
  }

  if( linm )
  {
    cout << "linmmatrix: " << *linmmatrix_ << endl;
    RCP<Epetra_Vector> reffcm = LINALG::CreateVector(*gmdofrowmap_,true);
    mmatrix_->Multiply(true,*z_,*reffcm);
    reffcm->Scale(-1.0);
    cout << "reffcm: " << *reffcm << endl;
    interface_[0]->FDCheckPenaltyLinM(linmmatrix_, reffcm);
  }

  cout << "-- CONTACTFDPENALTYKC1 ----------------------" << endl;

#endif

  return;
}

/*----------------------------------------------------------------------*
 | evaluate frictional contact and create linear system gitterle   10/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::EvaluateFriction(RCP<LINALG::SparseMatrix> kteff,
                                               RCP<Epetra_Vector> feff)
{
  // this is almost the same as in the frictionless contact
	// whereas we chose the EvaluateContact routine with
	// one difference

  // check if friction should be applied
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

	// coulomb friction case
	if(ftype == INPAR::CONTACT::friction_coulomb or
		 ftype == INPAR::CONTACT::friction_stick)
	{
		EvaluateContact(kteff,feff);
  }
  else if (ftype == INPAR::CONTACT::friction_tresca)
	{
	  dserror("Error in AbstractStrategy::Evaluate: Penalty Strategy for"
	         " Tresca friction not yet implemented");
	}
	else
	  dserror("Error in AbstractStrategy::Evaluate: Unknown friction type");
}

/*----------------------------------------------------------------------*
 | reset penalty parameter to intial value                    popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::ResetPenalty()
{
  // reset penalty parameter in strategy
  Params().set<double>("PENALTYPARAM",InitialPenalty());
  Params().set<double>("PENALTYPARAMTAN",InitialPenaltyTan());


  // reset penalty parameter in all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->IParams().set<double>("PENALTYPARAM",InitialPenalty());
    interface_[i]->IParams().set<double>("PENALTYPARAMTAN",InitialPenaltyTan());
  }
  return;
}

/*----------------------------------------------------------------------*
 | intialize second, third,... Uzawa step                     popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::InitializeUzawa(RCP<LINALG::SparseMatrix> kteff,
                                               RCP<Epetra_Vector> feff)
{
  // remove old stiffness terms
  // (FIXME: redundant code to EvaluateContact(), expect for minus sign)
  
  kteff->UnComplete();
  
  kteff->Add(*lindmatrix_, false, -(1.0-alphaf_), 1.0);
  kteff->Add(*linmmatrix_, false, -(1.0-alphaf_), 1.0);
    
  RCP<LINALG::SparseMatrix> dtilde = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  RCP<LINALG::SparseMatrix> mtilde = rcp(new LINALG::SparseMatrix(*gmdofrowmap_));
  dtilde = LINALG::Multiply(*dmatrix_, false, *linzmatrix_, false);
  mtilde = LINALG::Multiply(*mmatrix_, true, *linzmatrix_, false);
  kteff->Add(*dtilde, false, -(1.0-alphaf_), 1.0);
  kteff->Add(*mtilde, false, (1.0-alphaf_), 1.0);
  
  kteff->Complete();
  
  // remove old force terms
  // (FIXME: redundant code to EvaluateContact(), expect for minus sign)
  
  RCP<Epetra_Vector> fcmdold = rcp(new Epetra_Vector(dold_->RowMap()));
  dold_->Multiply(false, *zold_, *fcmdold);
  RCP<Epetra_Vector> fcmdoldtemp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcmdold, *fcmdoldtemp);
  feff->Update(alphaf_, *fcmdoldtemp, 1.0);

  RCP<Epetra_Vector> fcmmold = rcp(new Epetra_Vector(mold_->DomainMap()));
  mold_->Multiply(true, *zold_, *fcmmold);
  RCP<Epetra_Vector> fcmmoldtemp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcmmold, *fcmmoldtemp);
  feff->Update(-alphaf_, *fcmmoldtemp, 1.0);
    
  RCP<Epetra_Vector> fcmd = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false, *z_, *fcmd);
  RCP<Epetra_Vector> fcmdtemp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcmd, *fcmdtemp);
  feff->Update(1-alphaf_, *fcmdtemp, 1.0);
  
  RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
  mmatrix_->Multiply(true, *z_, *fcmm);
  RCP<Epetra_Vector> fcmmtemp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcmm, *fcmmtemp);
  feff->Update(-(1-alphaf_), *fcmmtemp, 1.0);
  
  // reset some matrices
  lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  
  // reset nodal derivZ values
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    for (int j=0;j<interface_[i]->Discret().NumMyColNodes();++j)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(interface_[i]->Discret().lColNode(j));
      for (int k=0; k<(int)((node->GetDerivZ()).size()); ++k)
        (node->GetDerivZ())[k].clear();
      (node->GetDerivZ()).resize(0);
    }
  }
  
  // now redo Initialize()
  Initialize();
  
  // and finally redo Evaluate()
  Evaluate(kteff,feff);
    
  return;
}

/*----------------------------------------------------------------------*
 | evaluate L2-norm of active constraints                     popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::UpdateConstraintNorm(int uzawaiter)
{
  // initialize parameters
  double cnorm = 0.0;
  double cnormtan = 0.0;
  bool updatepenalty = false;
  bool updatepenaltytan = false;
  double ppcurr = Params().get<double>("PENALTYPARAM");
  double ppcurrtan = Params().get<double>("PENALTYPARAMTAN");
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(Params(),"CONTACT");

  // gactivenodes_ is undefined
  if (gactivenodes_==Teuchos::null)
  {
    ConstraintNorm()=0;
    ConstraintNormTan()=0;
  }

  // gactivenodes_ has no elements
  else if (gactivenodes_->NumGlobalElements()==0)
  {
    ConstraintNorm()=0;
    ConstraintNormTan()=0;
  }

  // gactivenodes_ has at least one element
  else
  {
    // export weighted gap vector to gactiveN-map
    RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gact->GlobalLength()) LINALG::Export(*g_,*gact);

    // compute constraint norm
    gact->Norm2(&cnorm);

    // Evaluate norm in tangential direction for frictional contact
    if (ctype!=INPAR::CONTACT::contact_normal)
    {
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->EvaluateTangentNorm(cnormtan);

      cnormtan = sqrt(cnormtan);
    }

    //********************************************************************
    // adaptive update of penalty parameter
    // (only for Augmented Lagrange strategy)
    //********************************************************************
    INPAR::CONTACT::SolvingStrategy soltype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");

    if (soltype==INPAR::CONTACT::solution_auglag)
    {
      // check convergence of cnorm and update penalty parameter
      // only do this for second, third, ... Uzawa iteration
      // cf. Wriggers, Computational Contact Mechanics, 2nd edition (2006), p. 340
      if ((uzawaiter >= 2) && (cnorm > 0.25 * ConstraintNorm()))
      {
        updatepenalty = true;

        // update penalty parameter in strategy
        Params().set<double>("PENALTYPARAM",10*ppcurr);

        // update penalty parameter in all interfaces
        for (int i=0; i<(int)interface_.size(); ++i)
        {
          double ippcurr = interface_[i]->IParams().get<double>("PENALTYPARAM");
          if (ippcurr != ppcurr) dserror("Something wrong with penalty parameter");
          interface_[i]->IParams().set<double>("PENALTYPARAM",10*ippcurr);
        }
      }
    }

    if (soltype==INPAR::CONTACT::solution_auglag and ctype!=INPAR::CONTACT::contact_normal)
    {
      if ((uzawaiter >= 2) && (cnormtan > 0.25 * ConstraintNorm()))
      {
        updatepenaltytan = true;

        // update penalty parameter in strategy
        Params().set<double>("PENALTYPARAMTAN",10*ppcurrtan);

        // update penalty parameter in all interfaces
        for (int i=0; i<(int)interface_.size(); ++i)
        {
          double ippcurrtan = interface_[i]->IParams().get<double>("PENALTYPARAMTAN");
          if (ippcurrtan != ppcurrtan) dserror("Something wrong with penalty parameter");
          interface_[i]->IParams().set<double>("PENALTYPARAMTAN",10*ippcurrtan);
        }
      }
    }

    //********************************************************************
    // update constraint norm
    ConstraintNorm() = cnorm;
    ConstraintNormTan() = cnormtan;
  }

  // output to screen
  if (Comm().MyPID()==0)
  {
    cout << "********************************************\n";
    cout << "Normal Constraint Norm: " << cnorm << "\n";
    if (ctype != INPAR::CONTACT::contact_normal)
      cout << "Tangential Constraint Norm: " << cnormtan << "\n";
    if (updatepenalty)
      cout << "Updated normal penalty parameter: " << ppcurr << " -> " << Params().get<double>("PENALTYPARAM") << "\n";
    if (updatepenaltytan == true and ctype != INPAR::CONTACT::contact_normal)
     cout << "Updated tangential penalty parameter: " << ppcurrtan << " -> " << Params().get<double>("PENALTYPARAMTAN") << "\n";
    cout << "********************************************\n\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 | store Lagrange multipliers for next Uzawa step             popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::PenaltyStrategy::UpdateAugmentedLagrange()
{
  // store current LM into Uzawa LM
  zuzawa_ = rcp(new Epetra_Vector(*z_));

  return;
}

#endif // CCADISCRET
