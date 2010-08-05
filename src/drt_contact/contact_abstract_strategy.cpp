/*!----------------------------------------------------------------------
\file contact_abstract_strategy.cpp

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
#include "contact_abstract_strategy.H"
#include "contact_defines.H"
#include "contact_interface.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy::CoAbstractStrategy(RCP<Epetra_Map> problemrowmap, Teuchos::ParameterList params,
                                                vector<RCP<CONTACT::CoInterface> > interface,
                                                int dim, RCP<Epetra_Comm> comm, double alphaf) :
MORTAR::StrategyBase(problemrowmap,params,dim,comm,alphaf),
interface_(interface),
isincontact_(false),
isselfcontact_(false),
friction_(false),
dualquadslave3d_(false)
{
  // set potential global self contact status
  // (this is TRUE if at least one contact interface is a self contact interface)
  bool selfcontact = 0;
  for (int i=0;i<(int)interface_.size();++i)
    if (interface_[i]->SelfContact()) ++selfcontact;
  
  if (selfcontact) isselfcontact_=true;
  
  // check for infeasible self contact combinations
  INPAR::CONTACT::FrictionType ftype = Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(params,"FRICTION");
  if (isselfcontact_ && ftype != INPAR::CONTACT::friction_none)
    dserror("ERROR: Self contact only implemented for frictionless contact!");
  
  // set frictional contact status
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;
  
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

	// make numbering of LM dofs consecutive and unique across N interfaces
	int offset_if = 0;

  // merge interface maps to global maps
  for (int i=0; i<(int)interface_.size(); ++i)
  {
  	// build Lagrange multiplier dof map
		interface_[i]->UpdateLagMultSets(offset_if);

		// merge interface Lagrange multiplier dof maps to global LM dof map
		glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface_[i]->LagMultDofs());
		offset_if = glmdofrowmap_->NumGlobalElements();
		if (offset_if < 0) offset_if = 0;
        
    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_, interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_, interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, interface_[i]->MasterRowDofs());
  
    // merge active sets and slip sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    interface_[i]->InitializeActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map)
  gndofrowmap_ = LINALG::SplitMap(*problemrowmap_, *gsdofrowmap_);
  gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);

  // initialize flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  
  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // setup Lagrange muliplier vectors
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  // setup global Mortar matrices Dold and Mold
  dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  dold_->Zero();
  dold_->Complete();
  mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  mold_->Zero();
  mold_->Complete(*gmdofrowmap_, *gsdofrowmap_);

  // output contact stress vectors
  stressnormal_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  //----------------------------------------------------------------------
	// CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
	//----------------------------------------------------------------------
	// These matrices need to be applied to the slave displacements
	// in the cases of dual LM interpolation for tet10/hex20 meshes
	// in 3D. Here, the displacement basis functions have been modified
	// in order to assure positivity of the D matrix entries and at
	// the same time biorthogonality. Thus, to scale back the modified
	// discrete displacements \hat{d} to the nodal discrete displacements
	// {d}, we have to apply the transformation matrix T and vice versa
	// with the transformation matrix T^(-1).
	//----------------------------------------------------------------------
	INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
	if (shapefcn == INPAR::MORTAR::shape_dual)
		for (int i=0; i<(int)interface_.size(); ++i)
			dualquadslave3d_ += interface_[i]->Quadslave3d();

	//----------------------------------------------------------------------
	// IF SO, COMPUTE TRAFO MATRIX AND ITS INVERSE
	//----------------------------------------------------------------------
	if (Dualquadslave3d())
	{
		trafo_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
		invtrafo_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));

		// set of already processed nodes
		// (in order to avoid double-assembly for N interfaces)
		set<int> donebefore;

		// for all interfaces
		for (int i=0; i<(int)interface_.size(); ++i)
			interface_[i]->AssembleTrafo(*trafo_,*invtrafo_,donebefore);

		// FillComplete() transformation matrices
		trafo_->Complete();
		invtrafo_->Complete();
	}

	return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CoAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 | set current and old deformation state                     popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
  if( (statename=="displacement") || statename=="olddisplacement")
  {
    // set state on interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->SetState(statename, vec);
  }

  return;
}

/*----------------------------------------------------------------------*
 | update global master and slave sets (public)               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::UpdateMasterSlaveSetsGlobal()
{
  // reset global slave / master Epetra Maps
  gsnoderowmap_  = rcp(new Epetra_Map(0,0,Comm()));
  gsdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  gmdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  glmdofrowmap_  = rcp(new Epetra_Map(0,0,Comm()));

	// make numbering of LM dofs consecutive and unique across N interfaces
	int offset_if = 0;

  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i=0;i<(int)interface_.size();++i)
  {
  	// build Lagrange multiplier dof map
		interface_[i]->UpdateLagMultSets(offset_if);

		// merge interface Lagrange multiplier dof maps to global LM dof map
		glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface_[i]->LagMultDofs());
		offset_if = glmdofrowmap_->NumGlobalElements();
		if (offset_if < 0) offset_if = 0;

		// merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,interface_[i]->MasterRowDofs());
  }

  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate interface for next Newton step       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitEvalInterface()
{
  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();

    // evaluate interfaces
    interface_[i]->Evaluate();
  }

  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate mortar stuff for next Newton step    popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitEvalMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  if (IsSelfContact()) UpdateMasterSlaveSetsGlobal();

  // initialize Dold and Mold if not done already
  if (dold_==null)
  {
    dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    dold_->Zero();
    dold_->Complete();
  }
  if (mold_==null)
  {
    mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  }

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_       = LINALG::CreateVector(*gsnoderowmap_, true);
  if (friction_)
    jump_    = rcp(new Epetra_Vector(*gsdofrowmap_));

  // in the case of dual quad 3D, also the modified D matrices are setup  
  if (Dualquadslave3d())
   {  
    // initialize Dold and Mold if not done already
    if (doldmod_==null)
    {
      doldmod_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
   } 
  
  // (re)setup global matrices containing fc derivatives
  lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));

  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDM(*dmatrix_,*mmatrix_);
    interface_[i]->AssembleG(*g_);

#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    cout << " -- CONTACTFDNORMAL- -----------------------------------" << endl;
    interface_[i]->FDCheckNormalDeriv();
    cout << " -- CONTACTFDNORMAL- -----------------------------------" << endl;
#endif // #ifdef CONTACTFDNORMAL
  
#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    cout << " -- CONTACTFDMORTARD -----------------------------------" << endl;
    dmatrix_->Complete();
    if( dmatrix_->NormOne() )
      interface_[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    cout << " -- CONTACTFDMORTARD -----------------------------------" << endl;
#endif // #ifdef CONTACTFDMORTARD
  
#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    cout << " -- CONTACTFDMORTARM -----------------------------------" << endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if( mmatrix_->NormOne() )
        interface_[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    cout << " -- CONTACTFDMORTARM -----------------------------------" << endl;
#endif // #ifdef CONTACTFDMORTARM
  }

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate reference state                                gitterle 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateReferenceState(const RCP<Epetra_Vector> vec)
{
#ifndef CONTACTFORCEREFCONFIG 
  
	// only do something for frictional case
  if (!friction_) return;
  
  // set state and do mortar calculation
  SetState("displacement",vec);
  InitEvalInterface();
  InitEvalMortar();
  
  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // store D and M to old ones
  StoreDM("old");

  // store nodal entries from D and M to old ones
  StoreToOld(MORTAR::StrategyBase::dm);
  
  // transform dold_ in the case of dual quadratic 3d
  if (Dualquadslave3d())
  {
    RCP<LINALG::SparseMatrix> tempold = LINALG::MLMultiply(*dold_,false,*invtrafo_,false,false,false,true);
    doldmod_ = tempold;
  }  

  // evaluate relative movement
  EvaluateRelMov();

#else  
  
  // set state and do mortar calculation
  SetState("displacement",vec);
  InitEvalInterface();
  InitEvalMortar();
  
  //----------------------------------------------------------------------
	// CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
	//----------------------------------------------------------------------
	// Concretely, we apply the following transformations:
	// D         ---->   D * T^(-1)
	//----------------------------------------------------------------------
	if (Dualquadslave3d())
	{
		// modify dmatrix_
		RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
		dmatrix_ = temp;
	}

  // dump mortar Matrix D
  // this is always the matrix from the reference configuration,
  // also in the restart case
  dmatrix_->Dump("DREF");
  
  // only continue for frictional case
  if (!friction_) return;
  
  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // store D and M to old ones
  StoreDM("old");

  // store nodal entries from D and M to old ones
  StoreToOld(MORTAR::StrategyBase::dm);

  // evaluate relative movement
  EvaluateRelMov();
  
#endif
  
  return;
}

/*----------------------------------------------------------------------*
 | evaluate relative movement of contact bodies            gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateRelMov()
{
  
  if (friction_ == false)
    return;
  
  // transformation of slave displacement dofs
  // Dmod       ---->   D * T^(-1)
  if (Dualquadslave3d())
  {
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrixmod_ = temp;
  }
  
  // vector of slave coordinates xs
  RCP<Epetra_Vector> xsmod = rcp(new Epetra_Vector(*gsdofrowmap_));
 
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->AssembleSlaveCoord(xsmod);
 
  // in case of 3D dual quadratic case, slave coordinates xs are modified
  if (Dualquadslave3d())
    invtrafo_->Multiply(false,*xsmod,*xsmod);
  
  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->EvaluateRelMov(xsmod,dmatrixmod_,doldmod_);
    interface_[i]->AssembleRelMov(*jump_);
  }
  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Evaluate(RCP<LINALG::SparseOperator>& kteff,
                                           RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis)
{
  // treat frictional and frictionless cases differently
  if (friction_) EvaluateFriction(kteff,feff);
  else           EvaluateContact(kteff,feff);

  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n interfaces!");

    // get global quantity to be stored in nodes
    RCP<Epetra_Vector> vectorglobal = null;
    switch (type)
    {
      case MORTAR::StrategyBase::lmcurrent:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case MORTAR::StrategyBase::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      case MORTAR::StrategyBase::activeold:
      {
        break;
      }
      case MORTAR::StrategyBase::jump:
      {
        vectorglobal = Jump();
        break;
      }
      default:
        dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
    } // switch

    // export global quantity to current interface slave dof row map
    RCP<Epetra_Map> sdofrowmap = interface_[i]->SlaveRowDofs();
    RCP<Epetra_Vector> vectorinterface = rcp(new Epetra_Vector(*sdofrowmap));

    if (vectorglobal != null) // necessary for case "activeold"
      LINALG::Export(*vectorglobal, *vectorinterface);

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      vector<int> locindex(dim);

      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(cnode->Dofs()[dof]);
        if (locindex[dof]<0) dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch(type)
        {
        case MORTAR::StrategyBase::lmcurrent:
        {
          cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmold:
        {
          cnode->MoData().lmold()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmuzawa:
        {
          cnode->MoData().lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmupdate:
        {
#ifndef CONTACTPSEUDO2D
          // throw a dserror if node is Active and DBC
          if (cnode->IsDbc() && cnode->Active())
            dserror("ERROR: Slave node %i is active AND carries D.B.C.s!", cnode->Id());

          // explicity set global Lag. Mult. to zero for D.B.C nodes
          if (cnode->IsDbc()) (*vectorinterface)[locindex[dof]] = 0.0;
#endif // #ifndef CONTACTPSEUDO2D

          // explicity set global Lag. Mult. to zero for inactive nodes
          // (this is what we wanted to enforce anyway before condensation)
          if (cnode->Active()==false)
            (*vectorinterface)[locindex[dof]] = 0.0;

          // store updated LM into node
          cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::activeold:
        {
          if (!friction_)
            dserror("ERROR: This should not be called for contact without friction");
          CONTACT::FriNode* frinode = static_cast<FriNode*>(cnode);
          frinode->Data().ActiveOld() = frinode->Active();
          break;
        }
        case MORTAR::StrategyBase::jump:
        {
          if(!friction_)
            dserror("ERROR: This should not be called for contact without friction");
          FriNode* frinode = static_cast<FriNode*>(cnode);
          frinode->Data().jump()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        default:
          dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
        } // switch
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses        gitterle 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::OutputStresses ()
{
  // reset contact stress class variables
  stressnormal_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: OutputStresses: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      double nn[3];
      double nt1[3];
      double nt2[3];
      double lmn = 0.0;
      double lmt1 = 0.0;
      double lmt2 = 0.0;

      for (int j=0;j<3;++j)
      {
        nn[j]=cnode->MoData().n()[j];
        nt1[j]=cnode->CoData().txi()[j];
        nt2[j]=cnode->CoData().teta()[j];
        lmn +=  nn[j]* cnode->MoData().lm()[j];
        lmt1 += nt1[j]* cnode->MoData().lm()[j];
        lmt2 += nt2[j]* cnode->MoData().lm()[j];
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      vector<int> locindex(dim);

      // normal stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->Dofs()[dof]);
        (*stressnormal_)[locindex[dof]] = -lmn*nn[dof];
      }

      // tangential stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->Dofs()[dof]);
        (*stresstangential_)[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k=0;k<cnode->NumDof();++k)
      {
        int currdof = cnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid>=0)
        {
          cnode->SetDbc() = true;
          break;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDM(const string& state)
{
  //store Dold and Mold matrix in D and M
  if (state=="current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  }

  // store D and M matrix in Dold and Mold
  else if (state=="old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
    doldmod_=dmatrixmod_;
  }

  // unknown conversion
  else
    dserror("ERROR: StoreDM: Unknown conversion requested!");

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal quant. to old ones   (last conv. time step)gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreToOld(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreDMToNodes: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      switch (type)
      {
        case MORTAR::StrategyBase::dm:
        {
          // store D and M entries
          cnode->StoreDMOld();
          break;
        }
        case MORTAR::StrategyBase::pentrac:
        {
          // store penalty tractions to old ones
          cnode->StoreTracOld();
          break;
        }
        default:
          dserror("ERROR: StoreDMToNodes: Unknown state string variable!");
      } // switch
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Update(int istep, RCP<Epetra_Vector> dis)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (IsSelfContact()) zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  zold_->Update(1.0,*z_,0.0);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);
  StoreDM("old");

  // old displacements in nodes
  // (this is NOT only needed for friction but also for calculating
  // the auxiliary positions in binarytree contact search)
  SetState("olddisplacement",dis);
  
#ifdef MORTARGMSH1  
  VisualizeGmsh(istep);
#endif // #ifdef MORTARGMSH1

  // double-check if active set is really converged
  // (necessary e.g. for monolithic FSI with Lagrange multiplier contact,
  // because usually active set convergence check has been integrated into
  // structure Newton scheme, but now the monolithic FSI Newton scheme decides)
  if (!ActiveSetConverged() || !ActiveSetSemiSmoothConverged())
    dserror("ERROR: Active set not fully converged!");
  
  // reset active set status for next time step
  ResetActiveSet();

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // information and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information / quantities. 
  if (friction_)
  {
    // store contact state to contact nodes (active or inactive)
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
  
    // store nodal entries of D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // store nodal entries form penalty contact tractions to old ones
    StoreToOld(MORTAR::StrategyBase::pentrac);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoWriteRestart(RCP<Epetra_Vector>& activetoggle,
                                                 RCP<Epetra_Vector>& sliptoggle)
{
  // initalize
  activetoggle = rcp(new Epetra_Vector(*gsnoderowmap_));
  if (friction_)
    sliptoggle = rcp(new Epetra_Vector(*gsnoderowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoWriteRestart: Double active node check needed for n interfaces!");
        
    // loop over all slave nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = static_cast<CoNode*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      // set value active / inactive in toggle vector
      if (cnode->Active()) (*activetoggle)[dof]=1;
      if (friction_)
      {
        if (static_cast<CONTACT::FriNode*>(cnode)->Data().Slip()) (*sliptoggle)[dof]=1;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoReadRestart(IO::DiscretizationReader& reader,
                                                RCP<Epetra_Vector> dis)
{
  // set restart displacement state
  SetState("displacement", dis);
  SetState("olddisplacement", dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitEvalInterface();
  InitEvalMortar();

  //----------------------------------------------------------------------
	// CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
	//----------------------------------------------------------------------
	// Concretely, we apply the following transformations:
	// D         ---->   D * T^(-1)
	//----------------------------------------------------------------------
	if (Dualquadslave3d())
	{
		// modify dmatrix_
		RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
		dmatrix_ = temp;
	}

  // read restart information on actice set and slip set
  RCP<Epetra_Vector> activetoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
  reader.ReadVector(activetoggle,"activetoggle");
  
  // friction
  RCP<Epetra_Vector> sliptoggle;
  if(friction_)
  {  
    sliptoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
    reader.ReadVector(sliptoggle,"sliptoggle");
  }
  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoReadRestart: Double active node check needed for n interfaces!");
        
    // loop over all slave nodes on the current interface
    for (int j=0; j<(interface_[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof]==1)
      {
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = static_cast<CoNode*>(node);

        // set value active / inactive in cnode
        cnode->Active()=true;
        
        if (friction_)
        {
          // set value stick / slip in cnode
          if ((*sliptoggle)[dof]==1) static_cast<CONTACT::FriNode*>(cnode)->Data().Slip()=true;
        }
      }
    }
  }

  // read restart information on Lagrange multipliers
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  reader.ReadVector(LagrMult(),"lagrmultold");
  reader.ReadVector(LagrMultOld(),"lagrmultold");

  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // only for Augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  INPAR::CONTACT::SolvingStrategy st = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
  if (st == INPAR::CONTACT::solution_auglag)
  {
    zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    reader.ReadVector(LagrMultUzawa(),"lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  StoreDM("old");
  
  if (friction_)
  {
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
    StoreToOld(MORTAR::StrategyBase::dm);
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }  
  }

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements()) IsInContact()=true;

  return;
}


/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InterfaceForces(bool output)
{
  // Note that we ALWAYS use a TR-like approach to compute the contact
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}

  // check chosen output option
  INPAR::CONTACT::EmOutputType emtype =
    Teuchos::getIntegralValue<INPAR::CONTACT::EmOutputType>(Params(),"EMOUTPUT");
  
  // get out of here if no output wanted
  if (emtype==INPAR::CONTACT::output_none) return;
  
  // compute two subvectors of fc each via Lagrange multipliers z_n+1, z_n
  RCP<Epetra_Vector> fcslavetemp = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  RCP<Epetra_Vector> fcslavetempend = rcp(new Epetra_Vector(dold_->RowMap()));
  RCP<Epetra_Vector> fcmastertempend = rcp(new Epetra_Vector(mold_->DomainMap()));

  // for self contact, slave and master sets may have changed,
  // thus we have to export z and zold to new D and M dimensions
  if (IsSelfContact())
  {
    RCP<Epetra_Vector> zexp = rcp(new Epetra_Vector(dmatrix_->RowMap()));
    RCP<Epetra_Vector> zoldexp = rcp(new Epetra_Vector(dold_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    dmatrix_->Multiply(true,*zexp,*fcslavetemp);
    mmatrix_->Multiply(true,*zexp,*fcmastertemp);
    dold_->Multiply(true,*zoldexp,*fcslavetempend);
    mold_->Multiply(true,*zoldexp,*fcmastertempend);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->Multiply(true, *z_, *fcslavetemp);
    mmatrix_->Multiply(true, *z_, *fcmastertemp);
    dold_->Multiply(true, *zold_, *fcslavetempend);
    mold_->Multiply(true, *zold_, *fcmastertempend);
  }

  // export the contact forces to full dof layout
  RCP<Epetra_Vector> fcslave = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcslaveend = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmasterend = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);
  LINALG::Export(*fcslavetempend, *fcslaveend);
  LINALG::Export(*fcmastertempend, *fcmasterend);

  // contact forces and moments
  vector<double> gfcs(3);
  vector<double> ggfcs(3);
  vector<double> gfcm(3);
  vector<double> ggfcm(3);
  vector<double> gmcs(3);
  vector<double> ggmcs(3);
  vector<double> gmcm(3);
  vector<double> ggmcm(3);

  vector<double> gmcsnew(3);
  vector<double> ggmcsnew(3);
  vector<double> gmcmnew(3);
  vector<double> ggmcmnew(3);

  // weighted gap vector
  RCP<Epetra_Vector> gapslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapslave,posnode,lm,lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j=0;j<interface_[i]->MasterRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcmastertemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapmaster,posnode,lm,lmowner);
    }
  }

  // weighted gap
  RCP<Epetra_Vector> gapslavefinal  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmasterfinal = rcp(new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false,*gapslave,*gapslavefinal);
  mmatrix_->Multiply(false,*gapmaster,*gapmasterfinal);
  RCP<Epetra_Vector> gapfinal = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0,*gapslavefinal,0.0);
  gapfinal->Update(-1.0,*gapmasterfinal,1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> lm(3);
      vector<double> nodegaps(3);
      vector<double> nodegapm(3);

      // LMs and gaps
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = cnode->MoData().lm()[d];
      }

      // moments
      vector<double> nodemoments(3);
      vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1]*lm[2]-nodegaps[2]*lm[1];
      nodemoments[1] = nodegaps[2]*lm[0]-nodegaps[0]*lm[2];
      nodemoments[2] = nodegaps[0]*lm[1]-nodegaps[1]*lm[0];
      nodemomentm[0] = nodegapm[1]*lm[2]-nodegapm[2]*lm[1];
      nodemomentm[1] = nodegapm[2]*lm[0]-nodegapm[0]*lm[2];
      nodemomentm[2] = nodegapm[0]*lm[1]-nodegapm[1]*lm[0];
      for (int d=0;d<3;++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i=0;i<3;++i)
  {
    Comm().SumAll(&gfcs[i],&ggfcs[i],1);
    Comm().SumAll(&gfcm[i],&ggfcm[i],1);
    Comm().SumAll(&gmcs[i],&ggmcs[i],1);
    Comm().SumAll(&gmcm[i],&ggmcm[i],1);
    Comm().SumAll(&gmcsnew[i],&ggmcsnew[i],1);
    Comm().SumAll(&gmcmnew[i],&ggmcmnew[i],1);
  }

  // print interface results to file
  if (emtype == INPAR::CONTACT::output_file ||
      emtype == INPAR::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && Comm().MyPID()==0)
    {
      FILE* MyFile = NULL;
      MyFile = fopen("o/scilab_output/OutputInterface.txt", "at+");
      
      if (MyFile)
      {
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggmcm[i]);
        //for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsfgh[i]);
        //for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsmgh[i]);
        fprintf(MyFile, "\n");
        fclose(MyFile);
      }
      else
        dserror("ERROR: File for writing meshtying forces could not be opened.");
    } 
  }
  
  // print interface results to screen
  if (emtype == INPAR::CONTACT::output_screen ||
      emtype == INPAR::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && Comm().MyPID()==0)
    {
    	double snorm = sqrt(ggfcs[0]*ggfcs[0]+ggfcs[1]*ggfcs[1]+ggfcs[2]*ggfcs[2]);
    	double mnorm = sqrt(ggfcm[0]*ggfcm[0]+ggfcm[1]*ggfcm[1]+ggfcm[2]*ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n",ggfcs[0],ggfcs[1],ggfcs[2], snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n",ggfcm[0],ggfcm[1],ggfcm[2], mnorm);
      printf("Slave Contact Moment:  % e  % e  % e\n",ggmcs[0],ggmcs[1],ggmcs[2]);
      //printf("Slave Contact Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Contact Moment: % e  % e  % e\n",ggmcm[0],ggmcm[1],ggmcm[2]);
      //printf("Master Contact Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | evaluate contact forces w.r.t. reference configuration     mgit 07/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ForceRefConfig()
{
  // multiply current D matrix with current LM
  RCP<Epetra_Vector> forcecurr = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*z_,*forcecurr);
  
  // read mortar matrix D of reference configuration 
  RCP<LINALG::SparseMatrix> dref = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  std::string inputname = "DREF";
  dref->Load(Comm(),inputname);
  dref->Complete();
  
  // LM in reference / current configuration
  RCP<Epetra_Vector> zref  = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> zcurr  = rcp(new Epetra_Vector(*z_));

  // solve with default solver
	LINALG::Solver solver(Comm());
	solver.Solve(dref->EpetraOperator(),zref,forcecurr,true);

  // store reference LM into global vector and nodes
	z_ = zref;
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
  
  // print message
  if (Comm().MyPID()==0)
    cout << "\n**** First output is w.r.t. the reference configuration! ****" << endl;

  // print active set
  PrintActiveSet();

  // restore current LM into global vector and nodes
  z_ = zcurr;
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------- CONTACT::CoAbstractStrategy\n"
       << "Contact interfaces: " << (int)interface_.size() << endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    cout << *(interface_[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::PrintActiveSet()
{
  // output message
	Comm().Barrier();
  if (Comm().MyPID()==0)
  {
    printf("\nActive contact set--------------------------------------------------------------\n");
	  fflush(stdout);
  }

	//**********************************************************************
	// detailed active set output
	//**********************************************************************

#ifdef CONTACTASOUTPUT
  // loop over all interfaces
	for (int i=0; i<(int)interface_.size(); ++i)
	{
		//if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

		// loop over all slave nodes on the current interface
		for (int j=0;j<interface_[i]->SlaveFullNodes()->NumMyElements();++j)
		{
			// gid of current node
			int gid = interface_[i]->SlaveFullNodes()->GID(j);
			DRT::Node* node = interface_[i]->Discret().gNode(gid);
			if (!node) dserror("ERROR: Cannot find node with gid %",gid);

			// introduce local integer variable status
			// (0=inactive, 1=active, 2=slip, 3=stick)
			// (this is necessary as all data will be written by proc 0, but
			// the knowledge of the above status ONLY exists on the owner
			// processor of the respective node. Thus this information must
			// also be communicated to proc 0 in addition to the actual data!)
			int status = 0;

			//--------------------------------------------------------------------
			// FRICTIONLESS CASE
			//--------------------------------------------------------------------
			if (!friction_)
			{
				// cast to CoNode
				CoNode* cnode = static_cast<CoNode*>(node);

				// initialize output variables
				double wgap = 0.0;
				double nz = 0.0;

				// do processing only for local owner proc
				if (interface_[i]->lComm()->MyPID()==interface_[i]->Procmap()[cnode->Owner()])
				{
					// compute weighted gap
					wgap = (*g_)[g_->Map().LID(gid)];

					// compute normal part of Lagrange multiplier
					for (int k=0;k<3;++k)
						nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

					// compute status
					if (cnode->Active()) status = 1;
				}

				// communicate (locally on interface)
				interface_[i]->lComm()->Broadcast(&wgap,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&nz,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&status,1,interface_[i]->Procmap()[cnode->Owner()]);

				// output is done by local proc 0
				if (interface_[i]->lComm()->MyPID()==0)
				{
					// print nodes of inactive set *************************************
					if (status==0)
					{
						printf("INACTIVE: %d \t wgap: % e \t lm: % e \n",gid,wgap,nz);
						fflush(stdout);
					}

					// print nodes of active set ***************************************
					else
					{
						printf("ACTIVE:   %d \t wgap: % e \t lm: % e \n",gid,wgap,nz);
						fflush(stdout);
					}
				}
			}

			//--------------------------------------------------------------------
			// FRICTIONAL CASE
			//--------------------------------------------------------------------
			else
			{
				// cast to CoNode and FriNode
				CoNode* cnode = static_cast<CoNode*>(node);
				FriNode* frinode = static_cast<FriNode*>(cnode);

				// initialize output variables
				double wgap = 0.0;
				double nz = 0.0;
				double tz = 0.0;
				double jumptxi = 0.0;
				double jumpteta = 0.0;

				// do processing only for local owner proc
				if (interface_[i]->lComm()->MyPID()==interface_[i]->Procmap()[cnode->Owner()])
				{
					// compute weighted gap
					wgap = (*g_)[g_->Map().LID(gid)];

					// compute normal part of Lagrange multiplier
					for (int k=0;k<3;++k)
						nz += frinode->MoData().n()[k] * frinode->MoData().lm()[k];

					// additional output quantities for friction
					double txiz = 0.0;
					double tetaz = 0.0;

					// compute tangential parts of Lagrange multiplier and jumps
					for (int k=0;k<Dim();++k)
					{
						txiz += frinode->CoData().txi()[k] * frinode->MoData().lm()[k];
						tetaz += frinode->CoData().teta()[k] * frinode->MoData().lm()[k];
						jumptxi += frinode->CoData().txi()[k] * frinode->Data().jump()[k];
						jumpteta += frinode->CoData().teta()[k] * frinode->Data().jump()[k];
					}

					// total tangential component
					tz = sqrt(txiz*txiz+tetaz*tetaz);

					// check for dimensions
					if (Dim()==2 and abs(jumpteta)>0.0001)
						dserror("Error: Jumpteta should be zero for 2D");

					// compute status
					if (cnode->Active())
					{
						if (frinode->Data().Slip()) status = 2;
						else                        status = 3;
					}
				}

				// communicate (locally on interface)
				interface_[i]->lComm()->Broadcast(&wgap,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&nz,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&tz,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&jumptxi,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&jumpteta,1,interface_[i]->Procmap()[cnode->Owner()]);
				interface_[i]->lComm()->Broadcast(&status,1,interface_[i]->Procmap()[cnode->Owner()]);

				// output is now done by local proc 0
				if (interface_[i]->lComm()->MyPID()==0)
				{
					// print nodes of slip set **************************************
					if (status == 2)
					{
						printf("SLIP:  %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \n",gid,nz,tz,jumptxi,jumpteta);
						fflush(stdout);
					}
					// print nodes of stick set *************************************
					else if (status == 3)
					{
						printf("STICK: %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \n",gid,nz,tz,jumptxi,jumpteta);
						fflush(stdout);
					}
				}
		  }
		}
	}

#else
	//**********************************************************************
  // reduced active set output
	//**********************************************************************

	// counters
	int activenodes = 0;
	int gactivenodes = 0;
	int inactivenodes = 0;
	int ginactivenodes = 0;
	int slipnodes = 0;
	int gslipnodes = 0;

	// loop over all interfaces
	for (int i=0; i<(int)interface_.size(); ++i)
	{
		//if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

		// loop over all slave nodes on the current interface
		for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
		{
			int gid = interface_[i]->SlaveRowNodes()->GID(j);
			DRT::Node* node = interface_[i]->Discret().gNode(gid);
			if (!node) dserror("ERROR: Cannot find node with gid %",gid);

			// increase active counters
			CoNode* cnode = static_cast<CoNode*>(node);
			if (cnode->Active()) activenodes   += 1;
			else                 inactivenodes += 1;

			// increase friction counters
			if (friction_)
			{
				FriNode* frinode = static_cast<FriNode*>(cnode);
				if (frinode->Data().Slip()) slipnodes += 1;
			}
		}
	}

	// sum among all processors
	Comm().SumAll(&activenodes,&gactivenodes,1);
	Comm().SumAll(&inactivenodes,&ginactivenodes,1);
	Comm().SumAll(&slipnodes,&gslipnodes,1);

	// print active set information
	if (Comm().MyPID()==0)
	{
		if (friction_)
		{
			cout << BLUE2_LIGHT  << "Total     SLIP nodes:\t" << gslipnodes << END_COLOR << endl;
			cout << BLUE2_LIGHT  << "Total    STICK nodes:\t" << gactivenodes-gslipnodes << END_COLOR << endl;
			cout << YELLOW_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR << endl;
		}
		else
		{
			cout << BLUE2_LIGHT <<  "Total   ACTIVE nodes:\t" << gactivenodes << END_COLOR << endl;
			cout << YELLOW_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR << endl;
		}
	}
#endif // #ifdef CONTACTASOUTPUT

  // output line
  Comm().Barrier();
	if (Comm().MyPID()==0)
	{
		printf("--------------------------------------------------------------------------------\n\n");
		fflush(stdout);
	}

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::VisualizeGmsh(const int step, const int iter)
{
  // visualization with gmsh
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->VisualizeGmsh(step, iter);
}

#endif // CCADISCRET
