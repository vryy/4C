/*!----------------------------------------------------------------------
\file drt_contact_abstract_strategy.cpp

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

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialComm.h"
#include "../drt_lib/linalg_utils.H"
#include "drt_contact_abstract_strategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/linalg_ana.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::AbstractStrategy::AbstractStrategy(RCP<Epetra_Map> problemrowmap, Teuchos::ParameterList params,
                                            vector<RCP<CONTACT::Interface> > interface, int dim, RCP<Epetra_Comm> comm,
                                            double alphaf) :
interface_(interface),
scontact_(params),
dim_(dim),
comm_(comm),
alphaf_(alphaf),
problemrowmap_(problemrowmap),
activesetconv_(false),
activesetsteps_(1),
isincontact_(false),
isselfcontact_(false)
{
  // set potential global self contact status
  // (this is TRUE if at least one contact interface is a self contact interface)
  bool selfcontact = 0;
  for (int i=0;i<(int)interface_.size();++i)
    if (interface_[i]->SelfContact()) ++selfcontact;
  
  if (selfcontact) isselfcontact_=true;
  
  // check for infeasible self contact combinations
  INPAR::CONTACT::SolvingStrategy stype =
        Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(params,"STRATEGY");
  
  if (selfcontact && stype != INPAR::CONTACT::solution_lagmult)
    dserror("ERROR: Self contact only implemented for Lagrange multiplier strategy!");
  
  INPAR::CONTACT::ContactType ctype =
        Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(params,"CONTACT");
  
  if (selfcontact > 0 && ctype != INPAR::CONTACT::contact_normal)
    dserror("ERROR: Self contact only implemented for frictionless contact!");
  
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------                     

  // merge interface maps to global maps
  for (int i=0; i<(int)interface_.size(); ++i)
  {
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
    gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
    gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map) 
  gndofrowmap_ = LINALG::SplitMap(*problemrowmap_, *gsdofrowmap_);
  gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);

  
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

  // friction: setup vector of displacement jumps (slave dof) 
  jump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  
  // output contact stress vectors 
  stressnormal_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = rcp(new Epetra_Vector(*gsdofrowmap_));
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::AbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}
                                            
/*----------------------------------------------------------------------*
 | set current and old deformation state                     popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
  if( (statename=="displacement") || (statename=="olddisplacement") )
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
void CONTACT::AbstractStrategy::UpdateMasterSlaveSetsGlobal()
{
  // reset global slave / master Epetra Maps
  gsnoderowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  gsdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  gmdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  
  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i=0;i<(int)interface_.size();++i)
  {
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,interface_[i]->MasterRowDofs());
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate interface for next Newton step       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::InitEvalInterface()
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
void CONTACT::AbstractStrategy::InitEvalMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them befor initialiting D,M etc.
  if (IsSelfContact()) UpdateMasterSlaveSetsGlobal();
  
  // intitialize Dold and Mold if not done already
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
  dmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_          = LINALG::CreateVector(*gsnoderowmap_,true);
   
  // (re)setup global matrices containing fc derivatives
  lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  
  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
    
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
 | evaluate relative movement of contact bodies            gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::EvaluateRelMov(RCP<Epetra_Vector> disi)
{
#ifdef CONTACTRELVELMATERIAL
  
  // extract slave displacements from disi
  RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*disis);

  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);

  // recover incremental jump (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  mhatmatrix_->Multiply(false,*disim,*incrjump_);
  incrjump_->Update(1.0,*disis,-1.0);

  // sum up incremental jumps from active set nodes
  jump_->Update(1.0,*incrjump_,1.0);
  
  // friction
  // store updaded jumps to nodes
  StoreNodalQuantities(AbstractStrategy::jump);

#else

  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->EvaluateRelMov();
    interface_[i]->AssembleRelMov(*jump_);
  }
  
#endif

  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::Evaluate(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff)
{
  // check if friction should be applied
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

   	if (ftype == INPAR::CONTACT::friction_tresca ||
        ftype == INPAR::CONTACT::friction_coulomb ||
  	    ftype == INPAR::CONTACT::friction_stick )
  	  EvaluateFriction(kteff,feff);
 	  // Frictionless contact case
  	else
  	  EvaluateContact(kteff,feff);

  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::StoreNodalQuantities(AbstractStrategy::QuantityType type)
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
      case AbstractStrategy::lmcurrent:
      {
        vectorglobal = LagrMult();
        break;
      }
      case AbstractStrategy::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case AbstractStrategy::activeold:
      {
        break;
      }
      case AbstractStrategy::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case AbstractStrategy::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      case AbstractStrategy::jump:
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
      CNode* cnode = static_cast<CNode*>(node);

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
        case AbstractStrategy::lmcurrent:
        {
          cnode->lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case AbstractStrategy::lmold:
        {
          cnode->lmold()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case AbstractStrategy::lmuzawa:
        {
          cnode->lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case AbstractStrategy::activeold:
        {
          cnode->ActiveOld() = cnode->Active();
          break;
        }
        case AbstractStrategy::lmupdate:
        {
          // print a warning if a non-DBC inactive dof has a non-zero value
          // (only in semi-smooth Newton case, of course!)
          // bool semismooth = Teuchos::getIntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
          //if (semismooth && !cnode->Dbc()[dof] && !cnode->Active() && abs((*vectorinterface)[locindex[dof]])>1.0e-8)
          //  cout << "***WARNING***: Non-D.B.C. inactive node " << cnode->Id() << " has non-zero Lag. Mult.: dof "
          //       << cnode->Dofs()[dof] << " lm " << (*vectorinterface)[locindex[dof]] << endl;

#ifndef CONTACTPSEUDO2D
          // throw a dserror if node is Active and DBC
          if (cnode->Dbc()[dof] && cnode->Active())
            dserror("ERROR: Slave Node %i is active and at the same time carries D.B.C.s!", cnode->Id());

          // explicity set global Lag. Mult. to zero for D.B.C nodes
          if (cnode->IsDbc())
            (*vectorinterface)[locindex[dof]] = 0.0;
#endif // #ifndef CONTACTPSEUDO2D

          // explicity set global Lag. Mult. to zero for inactive nodes
          // (this is what we wanted to enforce anyway before condensation)
          if (cnode->Active()==false)
            (*vectorinterface)[locindex[dof]] = 0.0;
                    
          // store updated LM into node
          cnode->lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case AbstractStrategy::jump:
        {
          cnode->jump()[dof] = (*vectorinterface)[locindex[dof]];
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
void CONTACT::AbstractStrategy::OutputStresses ()
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
      CNode* cnode = static_cast<CNode*>(node);

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
        nn[j]=cnode->n()[j];
        nt1[j]=cnode->txi()[j];
        nt2[j]=cnode->teta()[j];
        lmn +=  nn[j]* cnode->lm()[j];
        lmt1 += nt1[j]* cnode->lm()[j];
        lmt2 += nt2[j]* cnode->lm()[j];
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
void CONTACT::AbstractStrategy::StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps)
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
      CNode* cnode = static_cast<CNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k=0;k<cnode->NumDof();++k)
      {
        int currdof = cnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid>=0) cnode->Dbc()[k] = true;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store DM to Nodes (in vecor of last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::StoreDMToNodes(AbstractStrategy::QuantityType type)
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
      CNode* cnode = static_cast<CNode*>(node);

      switch (type)
      {
        case AbstractStrategy::dm:
        {
          // store D and M entries
          cnode->StoreDMOld();
          break;
        }
        case AbstractStrategy::pentrac:
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
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::StoreDM(const string& state)
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
  }

  // unknown conversion
  else
    dserror("ERROR: StoreDM: Unknown conversion requested!");

  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::Update(int istep)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (IsSelfContact()) zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  
  zold_->Update(1.0,*z_,0.0);
  StoreNodalQuantities(AbstractStrategy::lmold);
  StoreDM("old");

#ifdef CONTACTGMSH1
  VisualizeGmsh(istep);
#endif // #ifdef CONTACTGMSH1

  // reset active set status for next time step
  ActiveSetConverged() = false;
  ActiveSetSteps() = 1;
  
  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::DoWriteRestart(RCP<Epetra_Vector>& activetoggle, RCP<Epetra_Vector>& sliptoggle)
{
  activetoggle = rcp(new Epetra_Vector(*gsnoderowmap_));
  sliptoggle = rcp(new Epetra_Vector(*gsnoderowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CNode* cnode = static_cast<CNode*>(node);

      // set value active / inactive in toggle vector
      if (cnode->Active())
        (*activetoggle)[j]=1;
      if (cnode->Slip())
        (*sliptoggle)[j]=1;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::DoReadRestart(IO::DiscretizationReader& reader,
                                              RCP<Epetra_Vector> dis)
{
  // set restart displacement state
  SetState("displacement", dis);
  
  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitEvalInterface();
  InitEvalMortar(); 
  
  // read restart information on actice set and slip set
  RCP<Epetra_Vector> activetoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
  reader.ReadVector(activetoggle,"activetoggle");
  RCP<Epetra_Vector> sliptoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
  reader.ReadVector(sliptoggle,"sliptoggle");

  // store restart information on active set and slip set  
  // into nodes, therefore first loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0; j<(interface_[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);
      
      if ((*activetoggle)[dof]==1)
      {
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);
        CNode* cnode = static_cast<CNode*>(node);
        
        // set value active / inactive in cnode
        cnode->Active()=true;

        // set value stick / slip in cnode
        if ((*sliptoggle)[dof]==1)
          cnode->Slip()=true;
      }
    }
  }
 
  // read restart information on Lagrange multipliers
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  reader.ReadVector(LagrMultOld(),"lagrmultold");
  reader.ReadVector(LagrMult(),"lagrmultold");
  
  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(AbstractStrategy::lmold);
  StoreNodalQuantities(AbstractStrategy::lmcurrent);
    
  // store restart Mortar quantities into nodes
  StoreDM("old");
  StoreNodalQuantities(AbstractStrategy::activeold);
  StoreDMToNodes(AbstractStrategy::dm);
    
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
    gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
  }

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements()) IsInContact()=true;
  
  return;
}


/*----------------------------------------------------------------------*
 |  Compute contact forces (for debugging only)               popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::ContactForces(RCP<Epetra_Vector> fresm)
{
  // Note that we ALWAYS use a TR-like approach to compute the contact
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}

  // FIXME: fresm is only here for debugging purposes!
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
    LINALG::Export(*z_,*zexp);
    LINALG::Export(*zold_,*zoldexp);
    dmatrix_->Multiply(false,*zexp,*fcslavetemp);
    mmatrix_->Multiply(true,*zexp,*fcmastertemp);
    dold_->Multiply(false,*zoldexp,*fcslavetempend);
    mold_->Multiply(true,*zoldexp,*fcmastertempend);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->Multiply(false, *z_, *fcslavetemp);
    mmatrix_->Multiply(true, *z_, *fcmastertemp);
    dold_->Multiply(false, *zold_, *fcslavetempend);
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

  // build total contact force vector (TR-like!!!)
  fc_=fcslave;
  fc_->Update(-(1.0-alphaf_), *fcmaster, 1.0-alphaf_);
  fc_->Update(alphaf_, *fcslaveend, 1.0);
  fc_->Update(-alphaf_, *fcmasterend, 1.0);

  /*
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
  RCP<Epetra_Vector> fresmslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fresmmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  LINALG::Export(*fresm,*fresmslave);
  LINALG::Export(*fresm,*fresmmaster);

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
      CNode* cnode = static_cast<CNode*>(node);

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
      CNode* cnode = static_cast<CNode*>(node);

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
      CNode* cnode = static_cast<CNode*>(node);

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
        lm[d] = cnode->lm()[d];
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

      cout << "NORMAL: " << cnode->n()[0] << " " << cnode->n()[1] << " " << cnode->n()[2] << endl;
      cout << "LM:     " << lm[0] << " " << lm[1] << " " << lm[2] << endl;
      cout << "GAP:    " << nodegaps[0]-nodegapm[0] << " " << nodegaps[1]-nodegapm[1] << " " << nodegaps[2]-nodegapm[2] << endl;
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

  // output
  double slavenorm = 0.0;
  fcslavetemp->Norm2(&slavenorm);
  double slavenormend = 0.0;
  fcslavetempend->Norm2(&slavenormend);
  double fresmslavenorm = 0.0;
  fresmslave->Norm2(&fresmslavenorm);
  if (Comm().MyPID()==0)
  {
    cout << "Slave Contact Force Norm (n+1):  " << slavenorm << endl;
    cout << "Slave Contact Force Norm (n):  " << slavenormend << endl;
    cout << "Slave Residual Force Norm: " << fresmslavenorm << endl;
    cout << "Slave Contact Force Vector: " << ggfcs[0] << " " << ggfcs[1] << " " << ggfcs[2] << endl;
    cout << "Slave Contact Moment Vector: " << ggmcs[0] << " " << ggmcs[1] << " " << ggmcs[2] << endl;
    cout << "Slave Contact Moment Vector (2nd version): " << ggmcsnew[0] << " " << ggmcsnew[1] << " " << ggmcsnew[2] << endl;
  }
  double masternorm = 0.0;
  fcmastertemp->Norm2(&masternorm);
  double masternormend = 0.0;
  fcmastertempend->Norm2(&masternormend);
  double fresmmasternorm = 0.0;
  fresmmaster->Norm2(&fresmmasternorm);
  if (Comm().MyPID()==0)
  {
    cout << "Master Contact Force Norm (n+1): " << masternorm << endl;
    cout << "Master Contact Force Norm (n): " << masternormend << endl;
    cout << "Master Residual Force Norm " << fresmmasternorm << endl;
    cout << "Master Contact Force Vector: " << ggfcm[0] << " " << ggfcm[1] << " " << ggfcm[2] << endl;
    cout << "Master Contact Moment Vector: " << ggmcm[0] << " " << ggmcm[1] << " " << ggmcm[2] << endl;
    cout << "Master Contact Moment Vector (2nd version): " << ggmcmnew[0] << " " << ggmcmnew[1] << " " << ggmcmnew[2] << endl;
  }
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
  */

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "----------------------------------- CONTACT::AbstractStrategy\n"
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
void CONTACT::AbstractStrategy::PrintActiveSet()
{
  // get input parameter ctype
  INPAR::CONTACT::ContactType ctype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(Params(),"CONTACT");
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");

  if (Comm().MyPID()==0)
    cout << "Active contact set--------------------------------------------------------------\n";
  Comm().Barrier();

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
      CNode* cnode = static_cast<CNode*>(node);

      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }

      // friction
      double zt = 0.0;
      double ztxi = 0.0;
      double zteta = 0.0;
      double jumptxi = 0.0;
      double jumpteta = 0.0;

      if(ftype == INPAR::CONTACT::friction_tresca ||
         ftype == INPAR::CONTACT::friction_coulomb ||
         ftype == INPAR::CONTACT::friction_stick)
      {
        // compute tangential parts of Lagrange multiplier and jumps
        for (int k=0;k<Dim();++k)
        {
          ztxi += cnode->txi()[k] * cnode->lm()[k];
          zteta += cnode->teta()[k] * cnode->lm()[k];
          jumptxi += cnode->txi()[k] * cnode->jump()[k];
          jumpteta += cnode->teta()[k] * cnode->jump()[k];
        }
      
        zt = sqrt(ztxi*ztxi+zteta*zteta);
        
        // check for dimensions        
        if (Dim()==2 and abs(jumpteta)>0.0001)
        {
        	dserror("Error: Jumpteta should be zero for 2D");
        }
      }

      if (ctype == INPAR::CONTACT::contact_normal)
      {
        // print nodes of inactive set *************************************
        if (cnode->Active()==false)
        {
          printf("INACTIVE: %d \t wgap: %e \t lm: %e \n",gid,wgap,nz);
          fflush(stdout);
        }

        // print nodes of active set ***************************************
        else
        {
          printf("ACTIVE:   %d \t wgap: %e \t lm: %e \n",gid,wgap,nz);
          fflush(stdout);
        }
        
      }
      else
      if(cnode->Active())
      {
        if(cnode->Slip())
          cout << "SLIP " << gid << " Normal " << nz << " Tangential " << zt << " Movement x " << cnode->xspatial()[0] - cnode->X()[0] << endl;
        else
         cout << "STICK " << gid << " Normal " << nz << " Tangential " << zt << " Movement x " << cnode->xspatial()[0] - cnode->X()[0] << endl;
      }
    }
  }

  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::VisualizeGmsh(const int step, const int iter)
{
  // check for frictional contact
  bool fric = false;
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(Params(),"FRICTION");
  if (ftype != INPAR::CONTACT::friction_none) fric=true;

  // visualization with gmsh
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->VisualizeGmsh(interface_[i]->CSegs(), step, iter, fric);
}

#endif // CCADISCRET
