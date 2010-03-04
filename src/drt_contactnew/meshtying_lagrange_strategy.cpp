/*!----------------------------------------------------------------------
\file meshtying_lagrange_strategy.cpp

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

#include "Epetra_SerialComm.h"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_defines.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::MtLagrangeStrategy::MtLagrangeStrategy(DRT::Discretization& discret, RCP<Epetra_Map> problemrowmap,
                                                Teuchos::ParameterList params,
                                                vector<RCP<MORTAR::MortarInterface> > interface,
                                                int dim, RCP<Epetra_Comm> comm, double alphaf) :
MtAbstractStrategy(discret, problemrowmap, params, interface, dim, comm, alphaf)
{
  // empty constructor                                       
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::MortarCoupling(const RCP<Epetra_Vector> dis)
{
  // print message
  if(Comm().MyPID()==0)
  {
    cout << "Performing mortar coupling...............";
    fflush(stdout);
  }
       
  // refer call to parent class
  MtAbstractStrategy::MortarCoupling(dis);
  
  // from here on, everything is Lagrange specific
  // for splitting in Evaluatemeshtying, we need the combined sm rowmap
  gsmdofs_ = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  invd_ = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd_->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd_,false,*mmatrix_,false);
  
  // print message
  if(Comm().MyPID()==0) cout << "done!" << endl;
    
  return;
}

/*----------------------------------------------------------------------*
 |  mesh intialization for rotational invariance              popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::MeshInitialization()
{
  // print message
  if(Comm().MyPID()==0)
  {
    cout << "Performing mesh initialization...........";
    fflush(stdout);
  }
  
  // this is only so simple for DUAL shape functions
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN") != INPAR::MORTAR::shape_dual)
    dserror("Mesh initialization only implemented for DUAL Lagrange multiplier strategy!");
  
  //**********************************************************************
  // (1) get master positions on global level
  //**********************************************************************
  // fill Xmaster first
  RCP<Epetra_Vector> Xmaster = LINALG::CreateVector(*gmdofrowmap_, true);
  
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
  // this is trivial for dual Lagrange multipliers
  RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*gsdofrowmap_,true);
  mhatmatrix_->Multiply(false,*Xmaster,*Xslavemod);
  
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
 |  evaluate meshtying (public)                               popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::EvaluateMeshtying(RCP<LINALG::SparseOperator>& kteff,
                                                    RCP<Epetra_Vector>& feff,
                                                    RCP<Epetra_Vector> dis)
{   
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();
    
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);
  LINALG::SplitMatrix2x2(kteffmatrix,gsmdofs_,gndofrowmap_,gsmdofs_,gndofrowmap_,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  LINALG::SplitVector(*feff,*gsmdofs_,fsm,*gndofrowmap_,fn);

  // we want to split fsm into 2 groups s,m
  fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  
  // do the vector splitting sm -> s+m
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);

  // store some stuff for static condensation of LM
  fs_   = fs;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbar)*ksn
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhatmatrix_,true,*ksn,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbar)*ksm
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhatmatrix_,true,*ksm,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kms: add T(mbar)*kss
  RCP<LINALG::SparseMatrix> kmsmod = LINALG::Multiply(*mhatmatrix_,true,*kss,false,false);
  kmsmod->Add(*kms,false,1.0,1.0);
  kmsmod->Complete(kms->DomainMap(),kms->RowMap());

  // fn: nothing to do

  // fs: subtract alphaf * old interface forces (t_n)
  RCP<Epetra_Vector> tempvecs = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*zold_,*tempvecs);
  fs->Update(-alphaf_,*tempvecs,1.0);

  // fm: add alphaf * old interface forces (t_n)
  RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*zold_,*tempvecm);
  fm->Update(alphaf_,*tempvecm,1.0); 
  
  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhatmatrix_->Multiply(true,*fs,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including meshtying)            */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81,true,false,kteffmatrix->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);

  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmsmod,false,1.0,1.0);

  // add matrices D and M to kteffnew
  kteffnew->Add(*dmatrix_,false,1.0,1.0);
  kteffnew->Add(*mmatrix_,false,-1.0,1.0);
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add s subvector (constraints) to feffnew
  
  // VERSION 1: constraints for u (displacements)
  //     (+) rotational invariance
  //     (+) no initial stresses
  // VERSION 2: constraints for x (current configuration)
  //     (+) rotational invariance
  //     (+) no initial stresses
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
  // (nothing needs not be done, as the constraints for meshtying are
  // LINEAR w.r.t. the displacements and in the first step, dis is zero.
  // Thus, the right hand side of the constraint lines is ALWAYS zero!)
  
  /*
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
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*g_,*gexp);
  feffnew->Update(1.0,*gexp,1.0);
  */

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
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*g_,*gexp);
  feffnew->Update(1.0,*gexp,1.0);
#endif // #ifdef MESHTYINGUCONSTR
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  kteff = kteffnew;
  feff = feffnew;

  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::Recover(RCP<Epetra_Vector> disi)
{
  // extract slave displacements from disi
  RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
  if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);
  
  // extract other displacements from disi
  RCP<Epetra_Vector> disin = rcp(new Epetra_Vector(*gndofrowmap_));
  if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);

  /**********************************************************************/
  /* Update Lagrange multipliers z_n+1                                  */
  /**********************************************************************/
  
  // approximate update
  //invd_->Multiply(false,*fs_,*z_);

  // full update
  z_->Update(1.0,*fs_,0.0);
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gsdofrowmap_));
  kss_->Multiply(false,*disis,*mod);
  z_->Update(-1.0,*mod,1.0);
  ksm_->Multiply(false,*disim,*mod);
  z_->Update(-1.0,*mod,1.0);
  ksn_->Multiply(false,*disin,*mod);
  z_->Update(-1.0,*mod,1.0);
  dmatrix_->Multiply(false,*zold_,*mod);
  z_->Update(-alphaf_,*mod,1.0);
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  z_->Scale(1/(1-alphaf_));

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}

#endif // CCADISCRET
