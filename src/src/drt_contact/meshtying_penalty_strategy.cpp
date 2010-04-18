/*!----------------------------------------------------------------------
\file meshtying_penalty_strategy.cpp

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

#include "meshtying_penalty_strategy.H"
#include "meshtying_defines.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::MtPenaltyStrategy::MtPenaltyStrategy(DRT::Discretization& discret, RCP<Epetra_Map> problemrowmap,
                                              Teuchos::ParameterList params,
                                              vector<RCP<MORTAR::MortarInterface> > interface,
                                              int dim, RCP<Epetra_Comm> comm, double alphaf) :
MtAbstractStrategy(discret, problemrowmap, params, interface, dim, comm, alphaf)
{
  // initialize constraint norm and initial penalty
  constrnorm_ = 0.0;
  initialpenalty_ = Params().get<double>("PENALTYPARAM");
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::MortarCoupling(const RCP<Epetra_Vector> dis)
{
  // print message
  if(Comm().MyPID()==0)
  {
    cout << "Performing mortar coupling...............";
    fflush(stdout);
  }
    
  // refer call to parent class
  MtAbstractStrategy::MortarCoupling(dis);
 
  // initialize mortar matrix products
  mtm_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  mtd_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  dtm_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  dtd_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  
  // build mortar matrix products
  mtm_ = LINALG::MLMultiply(*mmatrix_,true,*mmatrix_,false,false,false,true);
  mtd_ = LINALG::MLMultiply(*mmatrix_,true,*dmatrix_,false,false,false,true);
  dtm_ = LINALG::MLMultiply(*dmatrix_,true,*mmatrix_,false,false,false,true);
  dtd_ = LINALG::MLMultiply(*dmatrix_,true,*dmatrix_,false,false,false,true);
  
  // print message
  if(Comm().MyPID()==0) cout << "done!" << endl;
    
  return;
}

/*----------------------------------------------------------------------*
 |  mesh intialization for rotational invariance              popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::MeshInitialization()
{
  // print message
  if(Comm().MyPID()==0)
  {
    cout << "Performing mesh initialization...........";
    fflush(stdout);
  }
  
  //**********************************************************************
  // (1) get master positions on global level
  //**********************************************************************
  // fill Xmaster first
  RCP<Epetra_Vector> Xmaster = LINALG::CreateVector(*gmdofrowmap_,true);
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all master row nodes on the current interface
    for (int j=0; j<interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);
      
      // do assembly (overwrite duplicate nodes)
      for (int k=0;k<Dim();++k)
      {
        int dof = mtnode->Dofs()[k];
        (*Xmaster)[(Xmaster->Map()).LID(dof)] = mtnode->X()[k];
      }
    }
  }
  
  //**********************************************************************
  // (2) solve for modified slave positions on global level
  //**********************************************************************
  // create linear problem
  RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*gsdofrowmap_,true);
  RCP<Epetra_Vector> rhs = LINALG::CreateVector(*gsdofrowmap_,true);
  mmatrix_->Multiply(false,*Xmaster,*rhs);
  
  // solve with default solver
  LINALG::Solver solver(Comm());
  solver.Solve(dmatrix_->EpetraOperator(),Xslavemod,rhs,true);
      
  //**********************************************************************
  // (3) perform mesh initialization node by node
  //**********************************************************************
  // this can be done in the AbstractStrategy now
  MtAbstractStrategy::MeshInitialization(Xslavemod);
  
  // print message
  if(Comm().MyPID()==0) cout << "done!\n" << endl;
      
  return;
}

/*----------------------------------------------------------------------*
 | evaluate meshtying and create linear system                popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::EvaluateMeshtying(RCP<LINALG::SparseOperator>& kteff,
                                                   RCP<Epetra_Vector>& feff,
                                                   RCP<Epetra_Vector> dis)
{
  // since we will modify the graph of kteff by adding additional
  // meshtyong stiffness entries, we have to uncomplete it
  kteff->UnComplete();
    
  /**********************************************************************/
  /* Global setup of kteff, feff (including meshtying)                  */
  /**********************************************************************/
  double pp = Params().get<double>("PENALTYPARAM");

  // add penalty meshtying stiffness terms
  kteff->Add(*mtm_,false,pp,1.0);
  kteff->Add(*mtd_,false,-pp,1.0);
  kteff->Add(*dtm_,false,-pp,1.0);
  kteff->Add(*dtd_,false,pp,1.0);
  
  // build constraint vector
    
  // VERSION 1: constraints for u (displacements)
  //     (-) no rotational invariance
  //     (+) no initial stresses
  // VERSION 2: constraints for x (current configuration)
  //     (+) rotational invariance
  //     (-) initial stresses
  // VERSION 3: mesh initialization (default)
  //     (+) rotational invariance
  //     (+) initial stresses
  
  // As long as we perform mesh initialization, there is no difference
  // between versions 1 and 2, as the constraints are then exactly
  // fulfilled in the reference configuration X already (version 3)!

#ifdef MESHTYINGUCONSTR
  //**********************************************************************
  // VERSION 1: constraints for u (displacements)
  //**********************************************************************
  RCP<Epetra_Vector> tempvec1 = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> tempvec2 = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dis,*tempvec1);
  dmatrix_->Multiply(false,*tempvec1,*tempvec2);
  g_->Update(-1.0,*tempvec2,0.0);
  RCP<Epetra_Vector> tempvec3 = rcp(new Epetra_Vector(*gmdofrowmap_));
  RCP<Epetra_Vector> tempvec4 = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dis,*tempvec3);
  mmatrix_->Multiply(false,*tempvec3,*tempvec4);
  g_->Update(1.0,*tempvec4,1.0);
  
#else
  //**********************************************************************
  // VERSION 2: constraints for x (current configuration)
  //**********************************************************************
  // fill xmaster
  RCP<Epetra_Vector> xmaster = LINALG::CreateVector(*gmdofrowmap_, true);
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all master row nodes on the current interface
    for (int j=0; j<interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);
      
      // prepare assembly   
      Epetra_SerialDenseVector val(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      
      for (int k=0;k<Dim();++k)
      {
        val[k] = mtnode->xspatial()[k];
        lm[k] = mtnode->Dofs()[k];
        lmowner[k] = mtnode->Owner();
      }
      
      // do assembly
      LINALG::Assemble(*xmaster,val,lm,lmowner);
    }
  }
  
  // fill xslave
  RCP<Epetra_Vector> xslave = LINALG::CreateVector(*gsdofrowmap_, true);
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all master row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);
      
      // prepare assembly   
      Epetra_SerialDenseVector val(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      
      for (int k=0;k<Dim();++k)
      {
        val[k] = mtnode->xspatial()[k];
        lm[k] = mtnode->Dofs()[k];
        lmowner[k] = mtnode->Owner();
      }
      
      // do assembly
      LINALG::Assemble(*xslave,val,lm,lmowner);
    }
  }
  
  // finally build constraint vector
  RCP<Epetra_Vector> tempvec1 = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*xslave,*tempvec1);
  g_->Update(-1.0,*tempvec1,0.0);
  RCP<Epetra_Vector> tempvec2 = rcp(new Epetra_Vector(*gsdofrowmap_));
  mmatrix_->Multiply(false,*xmaster,*tempvec2);
  g_->Update(1.0,*tempvec2,1.0); 
  //**********************************************************************
#endif // #ifdef MESHTYINGUCONSTR
  
  // update LM vector
  // (in the pure penalty case, zuzawa is zero)
  z_->Update(1.0,*zuzawa_,0.0);
  z_->Update(-pp,*g_,1.0);

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
  
  // add penalty meshtying force terms
  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*z_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  feff->Update(1.0,*fmexp,1.0);
  
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true,*z_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  feff->Update(-1.0,*fsexp,1.0);
  
  // add old contact forces (t_n)
  RCP<Epetra_Vector> fsold = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true,*zold_,*fsold);
  RCP<Epetra_Vector> fsoldexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fsold,*fsoldexp);
  feff->Update(alphaf_,*fsoldexp,1.0);

  RCP<Epetra_Vector> fmold = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*zold_,*fmold);
  RCP<Epetra_Vector> fmoldexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmold,*fmoldexp);
  feff->Update(-alphaf_,*fmoldexp,1.0);
  
  return;
}

/*----------------------------------------------------------------------*
 | intialize Uzawa step 2,3...                                popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::InitializeUzawa(RCP<LINALG::SparseOperator>& kteff,
                                                 RCP<Epetra_Vector>& feff)
{
  // remove penalty meshtying force terms
  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*z_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  feff->Update(-1.0,*fmexp,1.0);
  
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*z_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  feff->Update(1.0,*fsexp,1.0);
    
  // update LM vector
  double pp = Params().get<double>("PENALTYPARAM");
  z_->Update(1.0,*zuzawa_,0.0);
  z_->Update(-pp,*g_,1.0);
  
  // add penalty meshtying force terms
  RCP<Epetra_Vector> fmnew = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*z_,*fmnew);
  RCP<Epetra_Vector> fmexpnew = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmnew,*fmexpnew);
  feff->Update(1.0,*fmexpnew,1.0);
  
  RCP<Epetra_Vector> fsnew = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*z_,*fsnew);
  RCP<Epetra_Vector> fsexpnew = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fsnew,*fsexpnew);
  feff->Update(-1.0,*fsexpnew,1.0);
  
  return;
}

/*----------------------------------------------------------------------*
 | reset penalty parameter to intial value                    popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::ResetPenalty()
{
  // reset penalty parameter in strategy
  Params().set<double>("PENALTYPARAM",InitialPenalty());
  
  
  // reset penalty parameter in all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  { 
    interface_[i]->IParams().set<double>("PENALTYPARAM",InitialPenalty());
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | evaluate L2-norm of active constraints                     popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::UpdateConstraintNorm(int uzawaiter)
{
  // initialize parameters
  double cnorm = 0.0;
  bool updatepenalty = false;
  double ppcurr = Params().get<double>("PENALTYPARAM");   

  // compute constraint norm
  g_->Norm2(&cnorm);
  
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
  //********************************************************************
  
  // update constraint norm
  constrnorm_ = cnorm;
  
  // output to screen
  if (Comm().MyPID()==0)
  {
    cout << "********************************************\n";
    cout << "Constraint Norm: " << cnorm << "\n";
    if (updatepenalty)
      cout << "Updated penalty parameter: " << ppcurr << " -> " << Params().get<double>("PENALTYPARAM") << "\n";
    cout << "********************************************\n\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 | store Lagrange multipliers for next Uzawa step             popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtPenaltyStrategy::UpdateAugmentedLagrange()
{
  // store current LM into Uzawa LM
  zuzawa_ = rcp(new Epetra_Vector(*z_));
  
  return;
}

#endif // CCADISCRET
