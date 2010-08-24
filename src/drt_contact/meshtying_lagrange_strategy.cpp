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

#include <Teuchos_Time.hpp>
#include "Epetra_SerialComm.h"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_defines.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::MtLagrangeStrategy::MtLagrangeStrategy(DRT::Discretization& discret, RCP<Epetra_Map> problemrowmap,
                                                Teuchos::ParameterList params,
                                                vector<RCP<MORTAR::MortarInterface> > interface,
                                                int dim, RCP<Epetra_Comm> comm, double alphaf) :
MtAbstractStrategy(discret, problemrowmap, params, interface, dim, comm, alphaf)
{
	// empty constructor body
	return;
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

	// time measurement
	Comm().Barrier();
	const double t_start = Teuchos::Time::wallTime();
       
  // refer call to parent class
  MtAbstractStrategy::MortarCoupling(dis);

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
  mhatmatrix_ = LINALG::MLMultiply(*invd_,false,*mmatrix_,false,false,false,true);
  
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  // These modifications are applied once right here, thus the
  // following code (EvaluateMeshtying, Recover) remains unchanged.
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
		// modify dmatrix_, invd_ and mhatmatrix_
		RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
		RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*trafo_,false,*invd_,false,false,false,true);
		RCP<LINALG::SparseMatrix> temp3 = LINALG::MLMultiply(*trafo_,false,*mhatmatrix_,false,false,false,true);
		dmatrix_    = temp1;
		invd_       = temp2;
		mhatmatrix_ = temp3;
  }

  /**********************************************************************/
  /* Build constraint matrix (containing D and M)                       */
  /**********************************************************************/
  // case 1: saddle point system
  //    -> constraint matrix with rowmap=Problemmap, colmap=LMmap
  // case 2: two static condensations
  //    -> no explicit constraint matrix needed
  // case 3: one static condensation
  //    -> constraint matrix with rowmap=Problemmap, colmap=Slavemap
  /**********************************************************************/
  bool setup = true;
  INPAR::CONTACT::SystemType systype = Teuchos::getIntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");
  if (systype==INPAR::CONTACT::system_condensed)
  {
#ifdef MESHTYINGTWOCON
  	setup = false;
#endif // #ifdef MESHTYINGTWOCON
  }

  // build constraint matrix only if necessary
  if (setup)
  {
		// first setup
		RCP<LINALG::SparseMatrix> constrmt = rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
		constrmt->Add(*dmatrix_,true,1.0,1.0);
		constrmt->Add(*mmatrix_,true,-1.0,1.0);
		constrmt->Complete(*gsdofrowmap_,*gdisprowmap_);

		// transform constraint matrix
		if (systype==INPAR::CONTACT::system_condensed)
		{
			// transform parallel row / column distribution
		  // (only necessary in the MESHTYINGPAR case)
#ifdef MESHTYINGPAR
			conmatrix_ = MORTAR::MatrixRowColTransform(constrmt,problemrowmap_,pgsdofrowmap_);
#else
			conmatrix_ = constrmt;
#endif // #ifdef MESHTYINGPAR
		}
		else
		{
			// transform parallel row distribution
		  // (only necessary in the MESHTYINGPAR case)
#ifdef MESHTYINGPAR
			RCP<LINALG::SparseMatrix> temp = MORTAR::MatrixRowTransform(constrmt,problemrowmap_);
#else
			RCP<LINALG::SparseMatrix> temp = constrmt;
#endif // #ifdef MESHTYINGPAR

			// always transform column GIDs of constraint matrix
			conmatrix_ = MORTAR::MatrixColTransformGIDs(temp,glmdofrowmap_);
		}
  }

  // time measurement
  Comm().Barrier();
  const double t_end = Teuchos::Time::wallTime()-t_start;
  if (Comm().MyPID()==0) cout << "in...." << t_end << " secs........";

  // print message
  if(Comm().MyPID()==0) cout << "done!" << endl;

  return;
}

/*----------------------------------------------------------------------*
 |  mesh initialization for rotational invariance             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::MeshInitialization()
{
  // print message
  if(Comm().MyPID()==0)
  {
    cout << "Performing mesh initialization...........";
    fflush(stdout);
  }
  
  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

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
  // initialize modified slave positions
  RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*gsdofrowmap_,true);
    
  // shape function type
  INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
    
  // CASE A: DUAL LM SHAPE FUNCTIONS
  if (shapefcn == INPAR::MORTAR::shape_dual)
  {
    // this is trivial for dual Lagrange multipliers
    mhatmatrix_->Multiply(false,*Xmaster,*Xslavemod);
  }
  
  // CASE B: STANDARD LM SHAPE FUNCTIONS
  else if (shapefcn == INPAR::MORTAR::shape_standard)
  {
    // create linear problem
    RCP<Epetra_Vector> rhs = LINALG::CreateVector(*gsdofrowmap_,true);
    mmatrix_->Multiply(false,*Xmaster,*rhs);
    
    // solve with default solver
    LINALG::Solver solver(Comm());
    solver.Solve(dmatrix_->EpetraOperator(),Xslavemod,rhs,true);
  }
  
  //**********************************************************************
  // (3) perform mesh initialization node by node
  //**********************************************************************
  // this can be done in the AbstractStrategy now
  MtAbstractStrategy::MeshInitialization(Xslavemod);

  // time measurement
  Comm().Barrier();
  const double t_end = Teuchos::Time::wallTime()-t_start;
  if (Comm().MyPID()==0) cout << "in...." << t_end << " secs........";

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
  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  INPAR::CONTACT::SystemType systype = Teuchos::getIntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn!=INPAR::MORTAR::shape_dual) dserror("Condensation only for dual LM");
        
    //********************************************************************
    // VERSION 1: CONDENSE LAGRANGE MULTIPLIERS AND SLAVE DOFS
    //********************************************************************
#ifdef MESHTYINGTWOCON
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
#ifdef MESHTYINGPAR
    LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
#else
    LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
#endif // #ifdef MESHTYINGPAR

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
    LINALG::SplitVector(*problemrowmap_,*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);

    // we want to split fsm into 2 groups s,m
    fs = rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = rcp(new Epetra_Vector(*gmdofrowmap_));
    
    // do the vector splitting sm -> s+m
    LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

    // store some stuff for static condensation of LM
    fs_   = fs;
    ksn_  = ksn;
    ksm_  = ksm;
    kss_  = kss;

    /**********************************************************************/
    /* Build the constraint vector                                        */
    /**********************************************************************/
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
    // VERSION 1: constraints for u (displacements)
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
    */

#else
    // VERSION 2: constraints for x (current configuration)
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
  #endif // #ifdef MESHTYINGUCONSTR
    
    /**********************************************************************/
    /* Build the final K and f blocks                                     */
    /**********************************************************************/
    // knn: nothing to do

    // knm: add kns*mbar
    RCP<LINALG::SparseMatrix> knmmod = rcp(new LINALG::SparseMatrix(*gndofrowmap_,100));
    knmmod->Add(*knm,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> knmadd = LINALG::MLMultiply(*kns,false,*mhatmatrix_,false,false,false,true);
    knmmod->Add(*knmadd,false,1.0,1.0);
    knmmod->Complete(knm->DomainMap(),knm->RowMap());

    // kms: add T(mbar)*kss
    RCP<LINALG::SparseMatrix> kmsmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmsmod->Add(*kms,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmsadd = LINALG::MLMultiply(*mhatmatrix_,true,*kss,false,false,false,true);
    kmsmod->Add(*kmsadd,false,1.0,1.0);
    kmsmod->Complete(kms->DomainMap(),kms->RowMap());
          
    // kmn: add T(mbar)*ksn
    RCP<LINALG::SparseMatrix> kmnmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmnmod->Add(*kmn,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhatmatrix_,true,*ksn,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

    // kmm: add T(mbar)*ksm + kmsmod*mbar
    RCP<LINALG::SparseMatrix> kmmmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmmmod->Add(*kmm,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhatmatrix_,true,*ksm,false,false,false,true);
    kmmmod->Add(*kmmadd,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmmadd2 = LINALG::MLMultiply(*kmsmod,false,*mhatmatrix_,false,false,false,true);
    kmmmod->Add(*kmmadd2,false,1.0,1.0);
    kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

    // fn: subtract kns*inv(D)*g
#ifndef MESHTYINGUCONSTR
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gsdofrowmap_));
    invd_->Multiply(false,*g_,*tempvec);
    RCP<Epetra_Vector> fnmod = rcp(new Epetra_Vector(*gndofrowmap_));
    kns->Multiply(false,*tempvec,*fnmod);
    fnmod->Update(1.0,*fn,1.0);
#endif

    // fs: subtract alphaf * old interface forces (t_n)
    RCP<Epetra_Vector> tempvecs = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*zold_,*tempvecs);
    tempvecs->Update(1.0,*fs,-alphaf_);

    // fm: add alphaf * old interface forces (t_n)
    RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*zold_,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0); 
    
    // fm: add T(mbar)*fs
    RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
    mhatmatrix_->Multiply(true,*tempvecs,*fmmod);
    fmmod->Update(1.0,*fm,1.0);

    // fm: subtract kmsmod*inv(D)*g
#ifndef MESHTYINGUCONSTR
    RCP<Epetra_Vector> fmadd = rcp(new Epetra_Vector(*gmdofrowmap_));
    kmsmod->Multiply(false,*tempvec,*fmadd);
    fmmod->Update(1.0,*fmadd,1.0);
#endif  
    
    // build identity matrix for slave dofs
    RCP<Epetra_Vector> ones = rcp (new Epetra_Vector(*gsdofrowmap_));
    ones->PutScalar(1.0);
    RCP<LINALG::SparseMatrix> onesdiag = rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();
    
    /**********************************************************************/
    /* Global setup of kteffnew, feffnew (including meshtying)            */
    /**********************************************************************/
    RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81,true,false,kteffmatrix->GetMatrixtype()));
    RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);

#ifdef MESHTYINGPAR
   kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
   kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
   onesdiag = MORTAR::MatrixRowTransform(onesdiag,pgsdofrowmap_);
#endif // #ifdef MESHTYINGPAR

    // add n submatrices to kteffnew
    kteffnew->Add(*knn,false,1.0,1.0);
    kteffnew->Add(*knmmod,false,1.0,1.0);
    //kteffnew->Add(*kns,false,1.0,1.0);

    // add m submatrices to kteffnew
    kteffnew->Add(*kmnmod,false,1.0,1.0);
    kteffnew->Add(*kmmmod,false,1.0,1.0);
    //kteffnew->Add(*kmsmod,false,1.0,1.0);

    // add identitiy for slave increments
    kteffnew->Add(*onesdiag,false,1.0,1.0);
    
    // FillComplete kteffnew (square)
    kteffnew->Complete();
          
    // add n subvector to feffnew
    RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
#ifdef MESHTYINGUCONSTR
    LINALG::Export(*fn,*fnexp);
#else
    LINALG::Export(*fnmod,*fnexp);
#endif
    feffnew->Update(1.0,*fnexp,1.0);

    // add m subvector to feffnew
    RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fmmod,*fmmodexp);
    feffnew->Update(1.0,*fmmodexp,1.0);

    /**********************************************************************/
    /* Replace kteff and feff by kteffnew and feffnew                     */
    /**********************************************************************/
    kteff = kteffnew;
    feff = feffnew;
    
    
    //********************************************************************
    // VERSION 2: CONDENSE ONLY LAGRANGE MULTIPLIERS
    //********************************************************************
#else
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
#ifdef MESHTYINGPAR
    LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
#else
    LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
#endif // #ifdef MESHTYINGPAR

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
    LINALG::SplitVector(*problemrowmap_,*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);

    // we want to split fsm into 2 groups s,m
    fs = rcp(new Epetra_Vector(*gsdofrowmap_));
    fm = rcp(new Epetra_Vector(*gmdofrowmap_));
    
    // do the vector splitting sm -> s+m
    LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

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
    RCP<LINALG::SparseMatrix> kmnmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmnmod->Add(*kmn,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhatmatrix_,true,*ksn,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

    // kmm: add T(mbar)*ksm
    RCP<LINALG::SparseMatrix> kmmmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmmmod->Add(*kmm,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhatmatrix_,true,*ksm,false,false,false,true);
    kmmmod->Add(*kmmadd,false,1.0,1.0); 
    kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

    // kms: add T(mbar)*kss
    RCP<LINALG::SparseMatrix> kmsmod = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmsmod->Add(*kms,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmsadd = LINALG::MLMultiply(*mhatmatrix_,true,*kss,false,false,false,true);
    kmsmod->Add(*kmsadd,false,1.0,1.0);
    kmsmod->Complete(kms->DomainMap(),kms->RowMap());

    // fn: nothing to do

    // fs: subtract alphaf * old interface forces (t_n)
    RCP<Epetra_Vector> tempvecs = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*zold_,*tempvecs);
    tempvecs->Update(1.0,*fs,-alphaf_);

    // fm: add alphaf * old interface forces (t_n)
    RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*zold_,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0); 
    
    // fm: add T(mbar)*fs
    RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
    mhatmatrix_->Multiply(true,*tempvecs,*fmmod);
    fmmod->Update(1.0,*fm,1.0);

    /**********************************************************************/
    /* Global setup of kteffnew, feffnew (including meshtying)            */
    /**********************************************************************/
    RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81,true,false,kteffmatrix->GetMatrixtype()));
    RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);

#ifdef MESHTYINGPAR
    kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    kmsmod = MORTAR::MatrixRowTransform(kmsmod,pgmdofrowmap_);
#endif // #ifdef MESHTYINGPAR

    // add n submatrices to kteffnew
    kteffnew->Add(*knn,false,1.0,1.0);
    kteffnew->Add(*knm,false,1.0,1.0);
    kteffnew->Add(*kns,false,1.0,1.0);

    // add m submatrices to kteffnew
    kteffnew->Add(*kmnmod,false,1.0,1.0);
    kteffnew->Add(*kmmmod,false,1.0,1.0);
    kteffnew->Add(*kmsmod,false,1.0,1.0);

    // add matrices D and M to kteffnew
    kteffnew->Add(*conmatrix_,true,1.0,1.0);
    
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
#endif // #ifdef MESHTYINGTWOCON 
  }
  
  //**********************************************************************
  //**********************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // add meshtying force terms
    RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*z_,*fs);
    RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fs,*fsexp);
    feff->Update(-(1.0-alphaf_),*fsexp,1.0);
    
    RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*z_,*fm);
    RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fm,*fmexp);
    feff->Update(1.0-alphaf_,*fmexp,1.0);

    // add old contact forces (t_n)
    RCP<Epetra_Vector> fsold = rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*zold_,*fsold);
    RCP<Epetra_Vector> fsoldexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fsold,*fsoldexp);
    feff->Update(-alphaf_,*fsoldexp,1.0);

    RCP<Epetra_Vector> fmold = rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*zold_,*fmold);
    RCP<Epetra_Vector> fmoldexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*fmold,*fmoldexp);
    feff->Update(alphaf_,*fmoldexp,1.0);
  } 

  return;
}

/*----------------------------------------------------------------------*
 | Solve linear system of saddle point type                   popp 03/10|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::SaddlePointSolve(LINALG::Solver& solver,
                  RCP<LINALG::SparseOperator> kdd,  RCP<Epetra_Vector> fd,
                  RCP<Epetra_Vector>  sold, RCP<Epetra_Vector> dirichtoggle,
                  int numiter)
{
  //**********************************************************************
  // prepare saddle point system
  //**********************************************************************
  // get system type
  INPAR::CONTACT::SystemType systype = Teuchos::getIntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");
  
  // some pointers and variables
  RCP<LINALG::SparseMatrix> stiffmt   = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kdd);
  RCP<Epetra_Map>           dispmap   = problemrowmap_;
  RCP<Epetra_Map>           slavemap  = gsdofrowmap_;
  RCP<Epetra_Map>           mastermap = gmdofrowmap_;
  RCP<Epetra_Map>           lmmap     = glmdofrowmap_;
      
  // initialize merged system (matrix, rhs, sol)
  RCP<Epetra_Map>           mergedmap   = LINALG::MergeMap(dispmap,lmmap,false); 
  RCP<LINALG::SparseMatrix> mergedmt    = rcp(new LINALG::SparseMatrix(*mergedmap,100,false,true));
  RCP<Epetra_Vector>        mergedrhs   = LINALG::CreateVector(*mergedmap);
  RCP<Epetra_Vector>        mergedsol   = LINALG::CreateVector(*mergedmap);
  RCP<Epetra_Vector>        mergedzeros = LINALG::CreateVector(*mergedmap);

  //**********************************************************************
  // finalize matrix and vector blocks
  //**********************************************************************
  // get constraint matrix
  RCP<LINALG::SparseMatrix> constrmt = conmatrix_;

  // remove meshtying force terms again
  // (solve directly for z_ and not for increment of z_)
  RCP<Epetra_Vector> fs = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true,*z_,*fs);
  RCP<Epetra_Vector> fsexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fs,*fsexp);
  fd->Update((1.0-alphaf_),*fsexp,1.0);
  
  RCP<Epetra_Vector> fm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*z_,*fm);
  RCP<Epetra_Vector> fmexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fm,*fmexp);
  fd->Update(-(1.0-alphaf_),*fmexp,1.0);
  
  // build constraint rhs (=empty)
  RCP<Epetra_Vector> constrrhs = rcp(new Epetra_Vector(*lmmap));
#ifndef MESHTYINGUCONSTR
  dserror("ERROR: Meshtying saddle point system only implemented for MESHTYINGUCONSTR");
#endif // #ifndef MESHTYINGUCONSTR
  
  //**********************************************************************
  // Build and solve saddle point system
  // (A) Standard coupled version
  //**********************************************************************
  if (systype==INPAR::CONTACT::system_spcoupled)
  {
    // build merged matrix    
    mergedmt->Add(*stiffmt,false,1.0,1.0);
    mergedmt->Add(*constrmt,false,1.0-alphaf_,1.0);
    mergedmt->Add(*constrmt,true,1.0,1.0);
    mergedmt->Complete();    
       
    // build merged rhs
    RCP<Epetra_Vector> fresmexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*fd,*fresmexp);
    mergedrhs->Update(1.0,*fresmexp,1.0);
    RCP<Epetra_Vector> constrexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*constrrhs,*constrexp);
    mergedrhs->Update(1.0,*constrexp,1.0);
    
    // adapt dirichtoggle vector and apply DBC
    RCP<Epetra_Vector> dirichtoggleexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*dirichtoggle,*dirichtoggleexp);
    LINALG::ApplyDirichlettoSystem(mergedmt,mergedsol,mergedrhs,mergedzeros,dirichtoggleexp);
    
    // standard solver call
    solver.Solve(mergedmt->EpetraMatrix(),mergedsol,mergedrhs,true,numiter==0);
  }
  
  //**********************************************************************
  // Build and solve saddle point system
  // (B) SIMPLER preconditioner version
  //**********************************************************************
  else if (systype==INPAR::CONTACT::system_spsimpler)
  {
    // build transposed constraint matrix
    RCP<LINALG::SparseMatrix> trconstrmt = rcp(new LINALG::SparseMatrix(*lmmap,100,false,true));
    trconstrmt->Add(*constrmt,true,1.0,0.0);
    trconstrmt->Complete(*dispmap,*lmmap);
    
    // apply Dirichlet conditions to (0,1) block
    RCP<Epetra_Vector> zeros   = rcp(new Epetra_Vector(*dispmap,true));
    RCP<Epetra_Vector> rhscopy = rcp(new Epetra_Vector(*fd));
    LINALG::ApplyDirichlettoSystem(stiffmt,sold,rhscopy,zeros,dirichtoggle);
    constrmt->ApplyDirichlet(dirichtoggle,false);
    
    // row map (equals domain map) extractor
    LINALG::MapExtractor rowmapext(*mergedmap,lmmap,dispmap);
    LINALG::MapExtractor dommapext(*mergedmap,lmmap,dispmap);

    // make solver SIMPLER-ready
    solver.PutSolverParamsToSubParams("SIMPLER", DRT::Problem::Instance()->FluidPressureSolverParams());
    
    // build block matrix for SIMPLER
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat =
      rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));
    mat->Assign(0,0,View,*stiffmt);
    mat->Assign(0,1,View,*constrmt);
    mat->Assign(1,0,View,*trconstrmt);
    mat->Complete();
    
    // we also need merged rhs here
    RCP<Epetra_Vector> fresmexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*fd,*fresmexp);
    mergedrhs->Update(1.0,*fresmexp,1.0);
    RCP<Epetra_Vector> constrexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*constrrhs,*constrexp);
    mergedrhs->Update(1.0,*constrexp,1.0);
    
    // we need a dummy merged matrix here in order to be able
    // to apply Dirichlet B.C. to mergedrhs and mergedsol
    mergedmt->Complete();    

    // adapt dirichtoggle vector and apply DBC
    RCP<Epetra_Vector> dirichtoggleexp = rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*dirichtoggle,*dirichtoggleexp);
    LINALG::ApplyDirichlettoSystem(mergedmt,mergedsol,mergedrhs,mergedzeros,dirichtoggleexp);
       
    // SIMPLER preconditioning solver call
    solver.Solve(mat->EpetraOperator(),mergedsol,mergedrhs,true,numiter==0);   
  }
  
  //**********************************************************************
  // invalid system types
  //**********************************************************************
  else dserror("ERROR: Invalid system type in SaddlePontSolve");
  
  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  RCP<Epetra_Vector> sollm = rcp(new Epetra_Vector(*lmmap));
  LINALG::MapExtractor mapext(*mergedmap,dispmap,lmmap);
  mapext.ExtractCondVector(mergedsol,sold);
  mapext.ExtractOtherVector(mergedsol,sollm);
  sollm->ReplaceMap(*slavemap);
  z_->Update(1.0,*sollm,0.0);
  
  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                            popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtLagrangeStrategy::Recover(RCP<Epetra_Vector> disi)
{
  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  INPAR::CONTACT::SystemType systype = Teuchos::getIntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");
  
  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed)
  {
    // double-check if this is a dual LM system
    if (shapefcn!=INPAR::MORTAR::shape_dual) dserror("Condensation only for dual LM");
    
    // extract slave displacements from disi
    RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

    // extract master displacements from disi
    RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);
    
    // extract other displacements from disi
    RCP<Epetra_Vector> disin = rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);
    
#ifdef MESHTYINGTWOCON
    /**********************************************************************/
    /* Update slave increment \Delta d_s                                  */
    /**********************************************************************/
    mhatmatrix_->Multiply(false,*disim,*disis);

    // if constraint vector non-zero, we need an additonal term
#ifndef MESHTYINGUCONSTR
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gsdofrowmap_));
    invd_->Multiply(false,*g_,*tempvec);
    disis->Update(1.0,*tempvec,1.0);
#endif // #ifndef MESHTYINGUNCONSTR

    RCP<Epetra_Vector> disisexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*disis,*disisexp);
    disi->Update(1.0,*disisexp,1.0);
#endif // #ifdef MESHTYINGTWOCON
          
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
		dmatrix_->Multiply(true,*zold_,*mod);
		z_->Update(-alphaf_,*mod,1.0);
		RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
		invd_->Multiply(true,*zcopy,*z_);
		z_->Scale(1/(1-alphaf_));
  }
  
  //**********************************************************************
  //**********************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // do nothing (z_ was part of solution already)
  }
  

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}

#endif // CCADISCRET
