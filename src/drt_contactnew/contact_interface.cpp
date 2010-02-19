/*!----------------------------------------------------------------------
\file contact_interface.cpp
\brief One contact interface

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

#ifndef PARALLEL
#include "Epetra_SerialComm.h"
#endif
#include "contact_interface.H"
#include "contact_integrator.H"
#include "contact_coupling2d.H"
#include "contact_coupling3d.H"
#include "contact_defines.H"
#include "../drt_mortar/mortar_binarytree.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoInterface::CoInterface(const int id, const Epetra_Comm& comm,
                                  const int dim,
                                  const Teuchos::ParameterList& icontact,
                                  bool selfcontact) :
MORTAR::MortarInterface(id,comm,dim,icontact),
selfcontact_(selfcontact),
friction_(false)
{
  // set frictional contact status
  INPAR::CONTACT::FrictionType ftype = Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(icontact,"FRICTION");
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CoInterface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
    os << "Contact ";
  MORTAR::MortarInterface::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::CreateSearchTree()
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  
  // warning
#ifdef MORTARGMSHCTN
  if (Dim()==3 && Comm().MyPID()==0)
  {
    cout << "\n******************************************************************\n";
    cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    cout << "******************************************************************\n";
  }
#endif // #ifdef MORTARGMSHCTN

  // binary tree search
  if (SearchAlg()==INPAR::MORTAR::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (SelfContact())
    {
      // set state in interface to intialize all kinds of quantities
      RCP<Epetra_Vector> zero =rcp(new Epetra_Vector(*idiscret_->DofRowMap()));
      SetState("displacement",zero);

      // create fully overlapping map of all contact elements
      RCP<Epetra_Map> elefullmap = rcp(new Epetra_Map(*idiscret_->ElementColMap()));

      // create binary tree object for self contact search
      binarytreeself_ = rcp(new CONTACT::SelfBinaryTree(Discret(),elefullmap,Dim(),SearchParam()));

    }
    //*****TWO BODY CONTACT*****
    else
    {
      // create binary tree object for contact search and setup tree
      binarytree_ = rcp(new MORTAR::BinaryTree(Discret(),selecolmap_,melefullmap_,Dim(),SearchParam()));

      // initialize active contact nodes via binarytree
      // binarytree_->SearchContactInit(binarytree_->Sroot(), binarytree_->Mroot());
    }
  }
  
  // no binary tree search
  else
  {
    if (SelfContact()) dserror("ERROR: Binarytree search needed for self contact");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CoNode* node = static_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));

    // reset feasible projection status
    node->HasProj() = false;
   }
  
  // loop over procs modes nodes to reset (column map)
  for (int i=0;i<OldColNodes()->NumMyElements();++i)
  {
    int gid = OldColNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    
    //reset nodal normal and tangents and jumps
    for (int j=0;j<3;++j)
    {
      cnode->MoData().n()[j]=0.0;
      cnode->CoData().txi()[j]=0.0;
      cnode->CoData().teta()[j]=0.0;
    }

    // reset nodal Mortar maps
    for (int j=0;j<(int)((cnode->MoData().GetD()).size());++j)
      (cnode->MoData().GetD())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetM()).size());++j)
      (cnode->MoData().GetM())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetMmod()).size());++j)
      (cnode->MoData().GetMmod())[j].clear();
    
    (cnode->MoData().GetD()).resize(0);
    (cnode->MoData().GetM()).resize(0);
    (cnode->MoData().GetMmod()).resize(0);
  
    
    // reset derivative maps of normal vector
    for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
      (cnode->CoData().GetDerivN())[j].clear();
    (cnode->CoData().GetDerivN()).resize(0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
      (cnode->CoData().GetDerivTxi())[j].clear();
    (cnode->CoData().GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
      (cnode->CoData().GetDerivTeta())[j].clear();
    (cnode->CoData().GetDerivTeta()).resize(0);
    
    // reset derivative map of Mortar matrices
    (cnode->CoData().GetDerivD()).clear();
    (cnode->CoData().GetDerivM()).clear();
    
    // reset nodal weighted gap and derivative
    cnode->CoData().Getg() = 1.0e12;
    (cnode->CoData().GetDerivG()).clear();
    
    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((cnode->CoData().GetDerivZ()).size()); ++j)
      (cnode->CoData().GetDerivZ())[j].clear();
    (cnode->CoData().GetDerivZ()).resize(0);

    if (friction_)
    {  
      FriNode* frinode = static_cast<FriNode*>(cnode);

      // reset nodal Mortar maps (Petrov-Galerkin approach)
      for (int j=0;j<(int)((frinode->Data().GetDPG()).size());++j)
        (frinode->Data().GetDPG())[j].clear();
      for (int j=0;j<(int)((frinode->Data().GetMPG()).size());++j)
        (frinode->Data().GetMPG())[j].clear();
      
      (frinode->Data().GetDPG()).resize(0);
      (frinode->Data().GetMPG()).resize(0);
      
      // reset SNodes and Mnodes
      frinode->Data().GetSNodes().clear();
      frinode->Data().GetMNodes().clear();
      
      // reset derivative map of jump
      for (int j=0; j<(int)((frinode->Data().GetDerivJump()).size()); ++j)
        (frinode->Data().GetDerivJump())[j].clear();
      (frinode->Data().GetDerivJump()).resize(0);
      
      
      // reset derivative map of Mortar matrices
      (frinode->Data().GetDerivDPG()).clear();
      (frinode->Data().GetDerivMPG()).clear();
    }
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CoElement* element = static_cast<CONTACT::CoElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::EvaluateSearchBinarytree()
{
  // get out of here if not participating in interface
  if (!lComm()) return true;

  // *********************************************************************
  // Possible versions for self contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchContactCombined(), which
  //    does both bottom-up update and search. Then the dynamics master/slave
  //    assignment routine UpdateMasterSlaveSets() is called.
  //
  // *********************************************************************
  if (SelfContact())
  {
    // calculate minimal element length
    binarytreeself_->SetEnlarge(false);

    // search for contact with an combined algorithm
    binarytreeself_->SearchContactCombined();

    // update master/slave sets of interface
    UpdateMasterSlaveSets();
  }

  // *********************************************************************
  // Possible versions for 2-body contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchCOntactCombined(), which
  //    does buth top-down update (where necessary) and search.
  //
  // 2) Separate Update and Contact Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchContactSeparate(). Of course, the bottom-up
  //    update makes more sense here!
  //
  // *********************************************************************
  else
  {
    // calculate minimal element length
    binarytree_->SetEnlarge(false);

    // update tree in a top down way
    //binarytree_->UpdateTreeTopDown();

    // update tree in a bottom up way
    //binarytree_->UpdateTreeBottomUp();

  #ifdef MORTARGMSHCTN
    for (int i=0;i<(int)(binarytree_->ContactMap().size());i++)
      binarytree_->ContactMap()[i].clear();
    binarytree_->ContactMap().clear();
    binarytree_->ContactMap().resize(2);
  #endif // #ifdef MORTARGMSHCTN

    // search for contact with a separate algorithm
    //binarytree_->SearchSeparate();

    // search for contact with an combined algorithm
    binarytree_->SearchCombined();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate Mortar matrix D on slave element (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateSlave(MORTAR::MortarElement& sele)
{
  // create a CONTACT integrator instance with correct NumGP and Dim
  CONTACT::CoIntegrator integrator(shapefcn_,sele.Shape());

  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // do the element integration (integrate and linearize D)
  int nrow = sele.NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  integrator.IntegrateDerivSlave2D3D(sele,sxia,sxib,dseg);

  // do the assembly into the slave nodes
  integrator.AssembleD(Comm(),sele,*dseg);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateCoupling(MORTAR::MortarElement& sele,
                                             MORTAR::MortarElement& mele)
{
  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim()==2)
  {
    // create instance of coupling class
    CONTACT::CoCoupling2d coup(shapefcn_,Discret(),Dim(),sele,mele);
    
    // do coupling
    coup.EvaluateCoupling();
  }
  // ************************************************************** 3D ***
  else if (Dim()==3)
  {
    bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

    // ************************************************** quadratic 3D ***
    if (sele.IsQuad() || mele.IsQuad())
    {
      // only for auxiliary plane 3D version
      if (!auxplane) dserror("ERROR: Quadratic 3D coupling only for AuxPlane case!");

      // build linear integration elements from quadratic MortarElements
      vector<RCP<MORTAR::IntElement> > sauxelements(0);
      vector<RCP<MORTAR::IntElement> > mauxelements(0);
      SplitIntElements(sele,sauxelements);
      SplitIntElements(mele,mauxelements);

      // loop over all IntElement pairs for coupling
      for (int i=0;i<(int)sauxelements.size();++i)
      {
        for (int j=0;j<(int)mauxelements.size();++j)
        {
          // create instance of coupling class
          CONTACT::CoCoupling3dQuad coup(shapefcn_,Discret(),Dim(),true,auxplane,
                        sele,mele,*sauxelements[i],*mauxelements[j]);
          // do coupling
          coup.EvaluateCoupling();
        }
      }
    }

    // ***************************************************** linear 3D ***
    else
    {
      // create instance of coupling class
      CONTACT::CoCoupling3d coup(shapefcn_,Discret(),Dim(),false,auxplane,
                                 sele,mele);
      // do coupling
      coup.EvaluateCoupling();
    }
  }
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateKappaPenalty(CONTACT::CoElement& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // check for quadratic 3D case
  bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

  // ************************************************** quadratic 3D ***
  if (Dim()==3 && sele.IsQuad())
  {
    // only for auxiliary plane 3D version
    if (!auxplane) dserror("ERROR: Quadratic 3D contact only for AuxPlane case!");

    // build linear integration elements from quadratic CElements
    vector<RCP<MORTAR::IntElement> > sauxelements(0);
    SplitIntElements(sele,sauxelements);

    for (int i=0;i<(int)sauxelements.size();++i)
    {
      // do the int element integration of kappa and store into gap
#ifdef MORTARPETROVGALERKIN
      int nrow = sauxelements[i]->NumNode();
#else
      int nrow = sele.NumNode();
#endif // #ifdef MORTARPETROVGALERKIN
      RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

      // create a CONTACT integrator instance with correct NumGP and Dim
      CONTACT::CoIntegrator integrator(shapefcn_,sauxelements[i]->Shape());
      integrator.IntegrateKappaPenalty(sele,*(sauxelements[i]),sxia,sxib,gseg);

      // do the assembly into the slave nodes
#ifdef MORTARPETROVGALERKIN
      integrator.AssembleG(Comm(),*(sauxelements[i]),*gseg);
#else
      integrator.AssembleG(Comm(),sele,*gseg);
#endif //#ifdef MORTARPETROVGALERKIN
    }
  }

  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.NumNode();
    RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

    // create a CONTACT integrator instance with correct NumGP and Dim
    CONTACT::CoIntegrator integrator(shapefcn_,sele.Shape());
    integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

    // do the assembly into the slave nodes
    integrator.AssembleG(Comm(),sele,*gseg);
  }
  
  // Check if PETROVGALERKIN-approach is switched on for tri6/hex20
  // in the frictional case
#ifdef MORTARPETROVGALERKIN
#ifndef CONTACTPETROVGALERKINFRIC
  if (friction_ && (sele.Shape()==DRT::Element::tri6 || sele.Shape()==DRT::Element::quad8) )
    dserror("Frictional contact needs flag: PETROVGALERKINFRIC for tri6/quad8");
#endif
#endif

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateRelMov()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
  
  if (friction_ == false)
    dserror("Error in CoInterface::EvaluateRelMov(): Only evaluated for frictional contact");

  // parameters
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k=0;k<3;++k)
      nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

    vector <double> jump(dim);
    for(int dim=0;dim<Dim();dim++)
      jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    double kappa = cnode->Kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_penalty)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             Teuchos::getIntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")!=1)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             Teuchos::getIntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")==1)
    {
      if(nz - cn*gap > 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_auglag)
    {
      if(lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
      dserror("Error in Interface::EvaluateRelMov(): Solution strategy not known!");

    if(activeinfuture==true)
    {
#ifdef CONTACTPETROVGALERKINFRIC
      vector<map<int,double> > dmap = cnode->Data().GetDPG();
      vector<map<int,double> > dmapold = cnode->Data().GetDOldPG();
#else
      vector<map<int,double> > dmap = cnode->MoData().GetD();
      vector<map<int,double> > dmapold = cnode->Data().GetDOld();
#endif

      set <int> snodes = cnode->Data().GetSNodes();

      // check if there are entries in the old D map
      if(dmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old D-Map!");

      map<int,double>::iterator colcurr;
      set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = static_cast<CoNode*>(snode);
        const int* sdofs = csnode->Dofs();

        double dik = (dmap[0])[sdofs[0]];
        double dikold = (dmapold[0])[sdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<csnode->NumDof();++dim)
        {
            jump[dim]-=(dik-dikold)*(csnode->xspatial()[dim]);
        }
      } //  loop over adjacent slave nodes

#ifdef CONTACTPETROVGALERKINFRIC
      vector<map<int,double> > mmap = cnode->Data().GetMPG();
      vector<map<int,double> > mmapold = cnode->Data().GetMOldPG();
#else
      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->Data().GetMOld();
#endif

      set <int> mnodescurrent = cnode->Data().GetMNodes();
      set <int> mnodesold = cnode->Data().GetMNodesOld();

      // check if there are entries in the old M map
      if(mmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old M-Map!");

      set <int> mnodes;
      set<int>::iterator mcurr;

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
           jump[dim]+= (mik-mikold)*(cmnode->xspatial()[dim]);
        }
      } //  loop over master nodes

      // write it to nodes
      for(int dim=0;dim<Dim();dim++)
        cnode->Data().jump()[dim] = jump[dim];
 
      // linearization of jump vector
      /*** 01  **********************************************************/
      // loop over all slave nodes (find adjacent ones to this stick node)
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = static_cast<CoNode*>(snode);
        const int* sdofs = csnode->Dofs();

        double dik = (dmap[0])[sdofs[0]];
        double dikold=(dmapold[0])[sdofs[0]];

        for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
        {
          int col = csnode->Dofs()[dimrow];
          double val = -(dik-dikold);
          cnode->AddDerivJumpValue(dimrow,col,val);
        }
      }

      /*** 02  **********************************************************/
      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold=(mmapold[0])[mdofs[0]];

        for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
        {
            int col = cmnode->Dofs()[dimrow];
            double val = (mik-mikold);
            cnode->AddDerivJumpValue(dimrow,col,val);
        }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
#ifdef CONTACTPETROVGALERKINFRIC
      map<int,map<int,double> >& ddmap = cnode->Data().GetDerivDPG();
#else
      map<int,map<int,double> >& ddmap = cnode->CoData().GetDerivD();
#endif
      map<int,map<int,double> >::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr=ddmap.begin();dscurr!=ddmap.end();++dscurr)
      {
        int gid = dscurr->first;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = static_cast<CoNode*>(snode);
        double* dxi = csnode->xspatial();

        // compute entry of the current stick node / slave node pair
#ifdef CONTACTPETROVGALERKINFRIC
        map<int,double>& thisdmmap = cnode->Data().GetDerivDPG(gid);
#else
        map<int,double>& thisdmmap = cnode->CoData().GetDerivD(gid);
#endif

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            double val =-colcurr->second*dxi[dimrow];
            cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
#ifdef CONTACTPETROVGALERKINFRIC
      map<int,map<int,double> >& dmmap = cnode->Data().GetDerivMPG();
#else
      map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
#endif

      map<int,map<int,double> >::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
      {
        int gid = dmcurr->first;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
#ifdef CONTACTPETROVGALERKINFRIC
        map<int,double>& thisdmmap = cnode->Data().GetDerivMPG(gid);
#else
        map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);
#endif

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            double val =colcurr->second*mxi[dimrow];
            cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }
    } // active nodes
  } // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble relative movement / jump (global)             gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRelMov(Epetra_Vector& jumpglobal)
{
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    int dim = cnode->NumDof();
    double* jump = cnode->Data().jump();

    Epetra_SerialDenseVector jumpnode(dim);
    vector<int> jumpdof(dim);
    vector<int> jumpowner(dim);

    for( int k=0; k<dim; ++k )
    {
      jumpnode(k) = jump[k];
      jumpdof[k] = cnode->Dofs()[k];
      jumpowner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(jumpglobal, jumpnode, jumpdof, jumpowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateTangentNorm(double& cnormtan)
{
  cnormtan=0;

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some information form node
    double* n = cnode->MoData().n();
    int dim = cnode->NumDof();

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    Epetra_SerialDenseMatrix jumpvec(dim,1);
    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->Data().jump()[i];

    // jump vector
    Epetra_SerialDenseMatrix forcevec(dim,1);
    for (int i=0;i<dim;i++)
      forcevec(i,0) = cnode->MoData().lm()[i];

    // evaluate jump in tangential direction
    Epetra_SerialDenseMatrix jumptan(dim,1);
    jumptan.Multiply('N','N',1,tanplane,jumpvec,0.0);

    // norm of tangential jumps for stick nodes
    if (cnode->Active()== true and cnode->Data().Slip()==false)
    {
      for( int j=0;j<cnode->NumDof();++j)
        cnormtan+=jumptan(j,0)*jumptan(j,0);
    }
    else if (cnode->Active()== true and cnode->Data().Slip()==true)
    {
      double jumptxi = 0;
      double jumpteta = 0;
      double forcen = 0;
      double forcetxi = 0;
      double forceteta = 0;

      for (int i=0;i<dim;i++)
      {
        jumptxi+=cnode->CoData().txi()[i]*cnode->Data().jump()[i];
        jumpteta+=cnode->CoData().teta()[i]*cnode->Data().jump()[i];

        forcen+=cnode->MoData().n()[i]*cnode->MoData().lm()[i];
        forcetxi+=cnode->CoData().txi()[i]*cnode->MoData().lm()[i];
        forceteta+=cnode->CoData().teta()[i]*cnode->MoData().lm()[i];
      }

      //cout << "FACTOR-Direction " << (jumptxi/jumpteta)/(forcetxi/forceteta) << endl;
      //cout << "FACTOR-Magnitude" << (frcoeff*forcen)/(sqrt(forcetxi*forcetxi+forceteta*forceteta))<< endl;
    }
  } // loop over slave nodes

  // get cnorm from all procs
  double sumcnormtanallprocs=0.0;
  Comm().SumAll(&cnormtan,&sumcnormtanallprocs,1);
  cnormtan=sumcnormtanallprocs;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                 popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegNormalForces(bool& localisincontact,
                                                 bool& localactivesetchange)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
  
  // penalty parameter
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->CoData().Getg();

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    //double modgap = cnode->CoData().Getg();
    //if (abs(modgap) < 1.0e-10) modgap=0.0;

    double kappa = cnode->Kappa();

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitely to constant
    // and corresponding derivatives to zero

    for( int j=0;j<dim;++j)
      cnode->MoData().lm()[j] = i*j;

    cnode->CoData().GetDerivZ().clear();

    continue;
#endif

    //********************************************************************
    // Decision on active /  inactive nodes (regularization)
    //
    // CASE 1: Penalty approach
    // A node is activated if its weighted gap is negative or deactivated
    // if its gap is equal zero or positive.
    // -> the regularization reads: lambda_n = kappa * pp * < -gap >
    //
    // CASE 2: Augmented Lagrange approach
    // A node is activated if its Lagrange multiplier, stemming from the
    // last Uzawa Lagrange multiplier AND the current regularization is
    // negative or deactivated if its LM is equal zero or positive.
    // -> the regularization reads: lambda_n = < lmuzawa_n - kappa * pp * gap >
    //
    // As the Uzawa Lagrange multipliers are zero in the penalty approach,
    // the two cases can formally be treted identically, see below.
    // We do not need an explicit separation of cases!
    //
    //********************************************************************

    // Activate/Deactivate node and notice any change
    if( (cnode->Active() == false) && (lmuzawan - kappa * pp * gap >= 0) )
    {
        cnode->Active() = true;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now active (";
        //for( int j=0; j<dim; j++)
        //  cout << " " << cnode->Dofs()[j] << " ";
        //cout << ") gap=" << gap << endl;
    }

    else if( (cnode->Active() == true) && (lmuzawan - kappa * pp * gap < 0) )
    {
        cnode->Active() = false;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now inactive, gap=" << gap << endl;
    }
    //********************************************************************

    // Compute derivZ-entries with the Macauley-Bracket
    // of course, this is only done for active constraints in order
    // for linearization and r.h.s to match!
    if( cnode->Active()==true )
    {

//      cout << "GID " << gid << endl;
//      cout << "LMUZAWAN " << lmuzawan << endl;
//      cout << "GAP " << gap << endl;

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      // contribution of derivative of normal
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();
      map<int,double>::iterator ncurr;

      for( int j=0;j<dim;++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
          cnode->AddDerivZValue(j, gcurr->first, - kappa * pp * (gcurr->second) * normal[j]);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, - kappa * pp * gap * ncurr->second);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, + lmuzawan * ncurr->second);
      }
    }

    // be sure to remove all LM-related stuff from inactive nodes
    else
    {
      // clear lagrange multipliers
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;

      // clear derivz
      cnode->CoData().GetDerivZ().clear();

    } // Macauley-Bracket
  } // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                 gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesPenalty()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);
       tanplane(0,2)=  -(n[0]*n[2]);
       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
       tanplane(1,2)=  -(n[1]*n[2]);

       tanplane(2,0)=  -(n[2]*n[0]);
       tanplane(2,1)=  -(n[2]*n[1]);
       tanplane(2,2)= 1-(n[2]*n[2]);
     }
     else if (dim==2)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);

       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
     }
     else
       dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->Data().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // fill vector tractionold
    vector<double> tractionold(dim);
    for (int i=0;i<dim;i++)
      tractionold[i] = cnode->Data().tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=tractionold[i]+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))
    {
      //cout << "Node " << gid << " is stick" << endl;
      cnode->Data().Slip() = false;

      // in the stick case, traction is trailtraction
      for (int i=0;i<dim;i++)
        cnode->Data().traction()[i]=trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
      //cout << "Node " << gid << " is slip" << endl;
      cnode->Data().Slip() = true;

      // in the slip case, traction is evaluated with a return map algorithm
      for (int i=0;i<dim;i++)
        cnode->Data().traction()[i]=maxtantrac/magnitude*trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap)+maxtantrac/magnitude*trailtraction[j];
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->Data().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->Data().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /**************************************** deriv(tanplane).jump  ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->Data().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->Data().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }
    // slip nodes
    else if (cnode->Active() == true && cnode->Data().Slip()== true)
    {
      /******************** tanplane.deriv(jump).maxtantrac/magnidude ***/

      vector<map<int,double> >& derivjump = cnode->Data().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** deriv(tanplane).jump.maxtantrac/magnitude ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->Data().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->Data().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j, gcurr->first, - frcoeff*kappa * ppnor * (gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
      for (int dim=0;dim<cnode->NumDof();++dim)
        temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->Data().jump()[dim]*kappa*pptan;

        traction+= tractionold[dimout];

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->Data().jump()[dim]*kappa*pptan;

        traction+=tractionold[dimout];

        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
        {
          int col = colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-colcurr->second*n[dim]*cnode->Data().jump()[dim]*traction/magnitude*pptan*kappa;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->Data().jump()[dim]*kappa*pptan;

          traction += tractionold[dimout];

          for (int dim=0;dim<cnode->NumDof();++dim)
          {

          // loop over all entries of the current derivative map
          for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
          {
            int col = colcurr->first;

            double val =-colcurr->second*n[dimout]*cnode->Data().jump()[dim]*traction/magnitude*pptan*kappa;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
               double val1 = val*temp[dimrow];
               cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces (Aug. Lagr.)    gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesAugmented()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->Data().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // Evaluate trailtraction
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=lmuzawatan(i,0)+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))    {
      //cout << "Node " << gid << " is stick" << endl;
      cnode->Data().Slip() = false;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j];
    }
    else
    {
      //cout << "Node " << gid << " is slip" << endl;
      cnode->Data().Slip() = true;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j]*maxtantrac/magnitude;
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->Data().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->Data().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*(cnode->Data().jump()[dim]);
            val = val - (colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*(cnode->Data().jump()[dim]);
            val = val-(colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }

    // slip nodes
    else if (cnode->Active() == true && cnode->Data().Slip()== true)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->Data().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->Data().jump()[dim];
            val = (val - (colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->Data().jump()[dim];
            val = (val-(colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j,gcurr->first,- frcoeff*kappa*ppnor*(gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( colcurr = (derivn[j]).begin(); colcurr != (derivn[j]).end(); ++colcurr )
        {
          for( int k=0;k<cnode->NumDof();++k)
          {
            double val = frcoeff*(colcurr->second)*lmuzawa(j,0)*trailtraction[k]/magnitude;
            cnode->AddDerivZValue(k,colcurr->first,val);
          }
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
        for (int dim=0;dim<cnode->NumDof();++dim)
          temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->Data().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->Data().jump()[dim]*kappa*pptan);

        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
        {
          int col = colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-colcurr->second*n[dim]*(lmuzawa(dim,0)+cnode->Data().jump()[dim]*pptan*kappa)*traction/magnitude;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->Data().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =-colcurr->second*n[dimout]*(lmuzawa(dim,0)+cnode->Data().jump()[dim]*pptan*kappa)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble derivatives of lagrange multipliers              popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinZ(LINALG::SparseMatrix& linzglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinZ: Node ownership inconsistency!");

    // derivz is the vector<map> we want to assemble
    vector<map<int,double> >& derivz = cnode->CoData().GetDerivZ();

    if ( (int) derivz.size()>0 )
    {
      int rowsize = cnode->NumDof();
      int colsize = (int) derivz[0].size();

      // consistency check
      for (int j=0; j<rowsize-1; ++j)
        if ((int)derivz[j].size() != (int)derivz[j+1].size())
          dserror("ERROR: AssembleLinZ: Column dim. of nodal derivz-map is inconsistent!");

      map<int,double>::iterator colcurr;

      // loop over dofs
      for ( int k=0; k<rowsize; ++k )
      {
        int row = cnode->Dofs()[k]; // row index equals global dof index of this #i node's dof k
        int l = 0;

        // loop over all directional derivative entries using the map iterator
        for( colcurr = derivz[k].begin(); colcurr != derivz[k].end(); ++colcurr )
        {
          int col = colcurr->first; // col index equals global id of directional derivative component ,l
          double val = colcurr->second;
          linzglobal.Assemble(val, row, col);
          l++;
        }

        if( l != colsize )
          dserror("ERROR: AssembleLinZ: l = %i but colsize = %i",k,colsize);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices with nodal normals / tangents           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNT(LINALG::SparseMatrix& nglobal,
                                     LINALG::SparseMatrix& tglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleNT: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(1);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(1);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(i);
      lmrowownerT[0] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->MoData().n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(1,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->CoData().txi()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }

    else if (Dim()==3)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(2);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(2);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(2*i);
      lmrowT[1] = activet_->GID(2*i+1);
      lmrowownerT[0] = cnode->Owner();
      lmrowownerT[1] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->MoData().n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(2,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->CoData().txi()[j];
        Tnode(1,j) = cnode->CoData().teta()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }
    else
      dserror("ERROR: Dim() must be either 2D or 3D");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives           popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    map<int,double>& dgmap = cnode->CoData().GetDerivG();
    map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;
      //cout << "Assemble S: " << row << " " << col << " " << val << endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix P containing tangent derivatives          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleP(LINALG::SparseMatrix& pglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleP: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      vector<map<int,double> >& dtmap = cnode->CoData().GetDerivTxi();
      map<int,double>::iterator colcurr;
      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();
      int row = activet_->GID(i);

      if (mapsize==3) mapsize=2;

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleP: Column dim. of nodal DerivT-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions)
      for (int j=0;j<mapsize;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << row << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,row,col);
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsize);
      }
    }
    else if (Dim()==3)
    {
      // prepare assembly
      vector<map<int,double> >& dtximap = cnode->CoData().GetDerivTxi();
      vector<map<int,double> >& dtetamap = cnode->CoData().GetDerivTeta();
      map<int,double>::iterator colcurr;
      int colsizexi = (int)dtximap[0].size();
      int colsizeeta = (int)dtetamap[0].size();
      int mapsizexi = (int)dtximap.size();
      int mapsizeeta = (int)dtetamap.size();
      int rowxi = activet_->GID(2*i);
      int roweta = activet_->GID(2*i+1);

      for (int j=0;j<mapsizexi-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTXi-map is inconsistent!");

      for (int j=0;j<mapsizeeta-1;++j)
        if ((int)dtetamap[j].size() != (int)dtetamap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTEta-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions) for TXi
      for (int j=0;j<mapsizexi;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << rowxi << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,rowxi,col);
          ++k;
        }

        if (k!=colsizexi)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizexi);
      }

      // loop over all derivative maps (=dimensions) for TEta
      for (int j=0;j<mapsizeeta;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << roweta << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,roweta,col);
          ++k;
        }

        if (k!=colsizeeta)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizeeta);
      }
    }
    else
      dserror("ERROR: Dim() must be either 2 or 3!");

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices LinDM containing fc derivatives         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                                       LINALG::SparseMatrix& linmglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
  
  /********************************************** LinDMatrix **********/
  // This is easy and can be done without communication, as the global
  // matrix lind has the same row map as the storage of the derivatives
  // of the D matrix.
  /**********************************************************************/

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i) // jtilde
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    int dim = cnode->NumDof();

    // Mortar matrix D derivatives
    map<int,map<int,double> >& dderiv = cnode->CoData().GetDerivD();

    map<int,double>::iterator colcurr;
    map<int,map<int,double> >::iterator dcurr;

    for( int j=0; j<dim; ++j ) // j
    {
      int row = cnode->Dofs()[j];

      // loop over all slave nodes dof again using the map iterator
      for ( dcurr = dderiv.begin(); dcurr != dderiv.end(); ++dcurr  ) // ktilde
      {
        // get the corresponding node
        DRT::Node* snode = idiscret_->gNode(dcurr->first);
        if (!snode)
          dserror("ERROR: Cannot find node with gid %",dcurr->first);
        CoNode* csnode = static_cast<CoNode*>(snode);

        double* lm = csnode->MoData().lm();

        // get dderiv row
        map<int,double>& dderivrow = dderiv[dcurr->first];

        // loop over all directional derivative entries using the map iterator
        for( colcurr = dderivrow.begin(); colcurr != dderivrow.end(); ++colcurr ) // l
        {
          int col = colcurr->first;
          double val = lm[j] * (colcurr->second);
          lindglobal.Assemble(val, row, col);
        }
      }
    }

    /** old version, using dual shape functions **
    // Mortar matrix D derivatives
    map<int,double>& dderiv = cnode->CoData().GetDerivD();
    map<int,double>::iterator colcurr;

    // current Lagrange multipliers
    double* lm = cnode->MoData().lm();

    // loop over all slave node dofs for assembly
    for (int j=0; j<numdof; ++j)
    {
      int row = cnode->Dofs()[j]; // row index equals global dof index of this #i slave node's dof j

      cout << "  j=" << j  << " row=" << row ;

      // loop over all directional derivative entries
      for (colcurr=dderiv.begin(); colcurr!=dderiv.end(); ++colcurr)
      {
        int col = colcurr->first; // col index equals global id of directional derivative component ,l

        double val = lm[j] * (colcurr->second);

        cout << " | col=" << col << ", val=" << val << " | ";

        if (abs(val)>1.0e-12)
          lindglobal.Assemble(val, row, col);
      }

      cout << " " << endl;
    }
    */
  }

  /********************************************** LinMMatrix ************/
  // This is a quite complex task and we have to do a lot of parallel
  // communication here!!! The reason for this is that the directional
  // derivative of M is stored slave-node-wise meaning that we can only
  // address the derivative of M via the slave node row map! But now
  // we have to assemble into a matrix linm which is based on the
  // master node row map!!! This is not straight forward as we CANNOT
  // do this without communication: when we have collected one slave
  // node's contribution to linm we have to assemble this portion to
  // the respective master nodes in parallel! Only the processor owning
  // the master node can do this! Therefore, we communicate the derivatives
  // of M to all procs and then explicitly call the current master proc
  // to do the assembly into the sparse matrix linm.
  /**********************************************************************/

  // loop over all slave nodes (full map)
  for (int i=0;i<snodefullmap_->NumMyElements();++i)
  {
    int gid = snodefullmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    int dim = cnode->NumDof();
    
    // create variables
     map<int,double>::iterator colcurr;
     map<int,map<int,double> >::iterator mcurr;
     double* lm = NULL;
     int mastersize = 0;
     map<int,map<int,double> > tempmap;
     map<int,map<int,double> >& mderiv = tempmap;

     // inter-proc. communication
     // (we want to keep all procs around, although at the moment only
     // the cnode owning proc does relevant work!)
     if (Comm().MyPID()==cnode->Owner())
     {
       // Mortar matrix M derivatives
       mderiv = cnode->CoData().GetDerivM();

       // current Lagrange multipliers
       lm = cnode->MoData().lm();

       // get sizes and iterator start
       mastersize = (int)mderiv.size();
       mcurr = mderiv.begin();
     }

    //********************************************************************
    // important notes:
    // 1) we use lComm here, as we only communicate among participating
    // procs of this interface (all others are already out of here)
    // 2) the broadcasting proc has to be given in lComm numbering, not
    // Comm numbering, which is achieved by calling the map procmap_
    //********************************************************************
    lComm()->Broadcast(&mastersize,1,procmap_[cnode->Owner()]);
    
    // loop over all master nodes in the DerivM-map of the current slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = 0;
      if (Comm().MyPID()==cnode->Owner())
      {
        mgid = mcurr->first;
        ++mcurr;
      }
      lComm()->Broadcast(&mgid,1,procmap_[cnode->Owner()]);

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
      CoNode* cmnode = static_cast<CoNode*>(mnode);
      
      // create variables
       int mapsize = 0;
       map<int,double> tempmap2;
       map<int,double>& thismderiv = tempmap2;

       // Mortar matrix M derivatives
       if (Comm().MyPID()==cnode->Owner())
       {
         thismderiv = cnode->CoData().GetDerivM()[mgid];
         mapsize = (int)(thismderiv.size());
       }
       lComm()->Broadcast(&mapsize,1,procmap_[cnode->Owner()]);
 
      // inner product M_{lj,c} * z_j for index j
      for (int j=0;j<dim;++j)
      {
        int row = cmnode->Dofs()[j];

        if (Comm().MyPID()==cnode->Owner())
          colcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = 0;
          double val = 0.0;
          if (Comm().MyPID()==cnode->Owner())
          {
            col = colcurr->first;
            val = lm[j] * (colcurr->second);
            ++colcurr;
          }
          // barrier here!
          lComm()->Barrier();
          lComm()->Broadcast(&col,1,procmap_[cnode->Owner()]);
          lComm()->Broadcast(&val,1,procmap_[cnode->Owner()]);

          // owner of master node has to do the assembly!!!
          if (Comm().MyPID()==cmnode->Owner())
          {
            //cout << "Assemble LinM: " << row << " " << col << " " << val << endl;
            if (abs(val)>1.0e-12) linmglobal.Assemble(-val,row,col);
          }
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (Comm().MyPID()==cnode->Owner() && colcurr!=thismderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
      }
    }

    // check for completeness of DerivM-Master-iteration
    if (Comm().MyPID()==cnode->Owner() && mcurr!=mderiv.end())
      dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble normal weighted gap                              popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleG(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDMG: Node ownership inconsistency!");

    /**************************************************** g-vector ******/
    if (cnode->CoData().Getg()!=0.0)
    {
      double gap = cnode->CoData().Getg();

      // cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj()
      //      << " IsActive: " << cnode->Active() << " Gap: " << gap << endl;

      // check if this active node has a feasible projection
      // else, one would at first think that a dserror has to be thrown!
      // (but this is not true in general, as there might indeed be an
      // active node which nervertheless has no feasible projection,
      // e.g. a slave node which is just over the edge of the master surface)
      // -> thus this check has been removed (popp 03/09)

      //if (!cnode->HasProj() && cnode->Active())
      //  dserror("ERROR: Active node ID: %i without feasible projection", cnode->Id());

      // check if this inactive node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      // when applying a Petrov-Galerkin scheme with standard shape functions
      // for the weighted gap interpolation this problem does not exist
      // FIXME: Only the linear case considered here (popp 04/09)
#ifndef MORTARPETROVGALERKIN
      if (!cnode->HasProj() && !cnode->Active())
      {
        gap = 1.0e12;
        cnode->CoData().Getg()=gap;
      }
#endif // #ifndef CONTACTPETROVGALERKIN

      Epetra_SerialDenseVector gnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      gnode(0) = gap;
      lm[0] = cnode->Id();
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(gglobal,gnode,lm,lmowner);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                          LINALG::SparseMatrix& linstickDISglobal,
                                          Epetra_Vector& linstickRHSglobal)
{
  // FIXGIT: Assemble LinStick is containing a matrix for the de-
  // rivatives of the Lagrange multipliers. This is according to Hüeber.
  // Because of worse convergence, this is not implemented, but the
  // code is commented after the algorithm.

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // create map of stick nodes
  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  // loop over all stick nodes of the interface
  for (int i=0;i<sticknodes->NumMyElements();++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    // prepare assembly, get information from node
    vector<map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
    vector<map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

   for (int j=0;j<Dim()-1;++j)
      if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

    if (Dim()==3)
    {
      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
    }

    // more information from node
    double* xi = cnode->xspatial();
    double* txi = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* jump = cnode->Data().jump();
    
    // iterator for maps
    map<int,double>::iterator colcurr;

    // row number of entries
    vector<int> row (Dim()-1);
    if (Dim()==2)
    {
      row[0] = stickt->GID(i);
    }
    else if (Dim()==3)
    {
      row[0] = stickt->GID(2*i);
      row[1] = stickt->GID(2*i)+1;
    }
    else
      dserror("ERROR: AssemblelinStick: Dimension not correct");

    // evaluation of specific components of entries to assemble
    double jumptxi=0;
    double jumpteta=0;
    for (int dim = 0;dim < Dim();dim++)
    {
      jumptxi += txi[dim]*jump[dim];
      jumpteta += teta[dim]*jump[dim];
    }

     // check for dimensions
     if(Dim()==2 and (jumpteta != 0.0))
      dserror ("ERROR: AssembleLinSlip: jumpteta must be zero in 2D");

    // Entries on right hand side
    /************************************************ (-utxi, -uteta) ***/
    Epetra_SerialDenseVector rhsnode(Dim()-1);
    vector<int> lm(Dim()-1);
    vector<int> lmowner(Dim()-1);

    rhsnode(0) = -jumptxi;
    lm[0] = cnode->Dofs()[1];
    lmowner[0] = cnode->Owner();

    if (Dim()==3)
    {
      rhsnode(1) = -jumpteta;
      lm[1] = cnode->Dofs()[2];
      lmowner[1] = cnode->Owner();
    }

    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);

    // Entries from differentiation with respect to displacements
    /*** 01 **************************************** tangent.(D-Dn-1) ***/

    // we need the nodal entries of the D-matrix and the old ones
    // they must have an entry
    vector<map<int,double> > dmapold=cnode->Data().GetDOld();
    if(dmapold.size()<1)
        dserror("ERROR: AssembleLinStick: No dmap-entries form previous time step");

    double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
    double Dold= dmapold[0][cnode->Dofs()[0]];

    // loop over dimensions
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
      int col = cnode->Dofs()[dim];
      double valtxi = (-1)*txi[dim]*(D-Dold);

      // do not assemble zeros into matrix
      if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

      if(Dim()==3)
      {
        double valteta = (-1)*teta[dim]*(D-Dold);
        if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }
    }

    /*** 02 **************************************** tangent.(M-Mn-1) ***/
    // we need the nodal entries of the M-matrix and the old ones
    vector<map<int,double> > mmap = cnode->MoData().GetM();
    vector<map<int,double> > mmapold = cnode->Data().GetMOld();

    // map from previous time step must have an entry
    if(mmapold.size()<1)
       dserror("ERROR: AssembleLinStick: No mmap-entries form previous time step");

    // create a set of nodes including nodes according to M entries
    // from current and previous time step
    set <int> mnodes;

    // iterator
    set<int>::iterator mcurr;

    set <int> mnodescurrent = cnode->Data().GetMNodes();
    set <int> mnodesold = cnode->Data().GetMNodesOld();

    for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
      mnodes.insert(*mcurr);

    for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
      mnodes.insert(*mcurr);

    // loop over all master nodes (find adjacent ones to this stick node)
    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
    {
      int gid = *mcurr;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cmnode = static_cast<FriNode*>(mnode);
      const int* mdofs = cmnode->Dofs();

      double mik = (mmap[0])[mdofs[0]];
      double mikold = (mmapold[0])[mdofs[0]];

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cmnode->Dofs()[dim];
        double valtxi = txi[dim]*(mik-mikold);

        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

        if (Dim()==3)
        {
          double valteta = teta[dim]*(mik-mikold);
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }
    }

    /*** 1 ******************************** deriv(tangent).(D-Dn-1).x ***/
    // loop over dimensions
    for (int j=0;j<Dim();++j)
    {
      // loop over all entries of the current derivative map (txi)
      for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
      {
        int col = colcurr->first;
        double val = (-1)*(D-Dold)*xi[j]*colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
      }

      if(Dim()==3)
      {
        // loop over all entries of the current derivative map (teta)
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*(D-Dold)*xi[j]*colcurr->second;

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
        }
      }
    }

    /*** 2 ******************************** deriv(tangent).(M-Mn-1).x ***/
    // loop over all master nodes (find adjacent ones to this stick node)
    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
    {
      int gid = *mcurr;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cmnode = static_cast<FriNode*>(mnode);
      const int* mdofs = cmnode->Dofs();

      double mik = (mmap[0])[mdofs[0]];
      double mikold = (mmapold[0])[mdofs[0]];
      double* mxi = cmnode->xspatial();

      // loop over dimensions
      for (int j=0;j<Dim();++j)
      {
        // loop over all entries of the current derivative map (dtxi)
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = (mik-mikold)*mxi[j]*colcurr->second;
          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
        }

        if(Dim()==3)
        {
          // loop over all entries of the current derivative map (dteta)
          for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (mik-mikold)*mxi[j]*colcurr->second;
            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
          }
        }
      }
    }

    /*** 3 ************************************** tangent.Deriv(D).x  ***/
    // we need the dot product txi.x and teta.x of this node
    double txidotx = 0.0;
    double tetadotx = 0.0;
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
      txidotx += txi[dim]*xi[dim];
      tetadotx += teta[dim]*xi[dim];
    }

    // prepare assembly
    map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

    // loop over all entries of the current derivative map
    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi = (-1)*txidotx*colcurr->second;

      if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

      if(Dim()==3)
      {
        double valteta = (-1)*tetadotx*colcurr->second;
        if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }
    }

    /*** 4 *************************************** tangent.Deriv(M).x ***/
    // we need the Lin(M-matrix) entries of this node
    map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
    map<int,map<int,double> >::iterator dmcurr;

    // loop over all master nodes in the DerivM-map of the active slave node
    for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
    {
      int gid = dmcurr->first;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cmnode = static_cast<FriNode*>(mnode);
      double* mxi = cmnode->xspatial();

      // we need the dot product ns*xm of this node pair
      double txidotx = 0.0;
      double tetadotx = 0.0;
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        txidotx += txi[dim]*mxi[dim];
        tetadotx += teta[dim]*mxi[dim];

      }
      // compute matrix entry of the current active node / master node pair
      map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

      // loop over all entries of the current derivative map
      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
      {
        int col = colcurr->first;
        double val = txidotx*colcurr->second;

        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);

        if (Dim()==3)
        {
          double val = tetadotx*colcurr->second;
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
        }
      }
    }
  }

//  have a look to the beginning of the function
//  // only for coulomb friction
//  string ftype   = IParams().get<string>("friction type","none");
//  double frcoeff = IParams().get<double>("friction coefficient",0.0);
//  double cn = IParams().get<double>("semismooth cn",0.0);
//  double ct = IParams().get<double>("semismooth ct",0.0);
//
//  if (ftype == "tresca")
//    dserror ("Error: AssemblelinStick: complementary function according"
//             " to Hueber only available for Coulomb friction");
//
//  // get out of here if not participating in interface
//  if (!lComm())
//    return;
//
//  // create map of stick nodes
//  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
//  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);
//
//  // nothing to do if no stick nodes
//  if (sticknodes->NumMyElements()==0)
//    return;
//
//  // not yet implemented for 3D
//    if (Dim()==3)
//      dserror("ERROR: AssembleLinStick: 3D not yet implemented");
//
//  // loop over all stick nodes of the interface
//  for (int i=0;i<sticknodes->NumMyElements();++i)
//  {
//    int gid = sticknodes->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
//    FriNode* cnode = static_cast<FriNode*>(node);
//
//    if (cnode->Owner() != Comm().MyPID())
//      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");
//
//    // prepare assembly
//    vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();
//    map<int,double>::iterator colcurr;
//
//    // calculate DerivT from DerivN
//    // only for 2D so far, in this case calculation is very easy
//    // dty =  dnx
//    // dtx = -dny
//
//    vector <map<int,double> > dtmap(Dim());
//
//    for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
//      dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));
//
//    for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
//      dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));
//
//    int colsize = (int)dtmap[0].size();
//    int mapsize = (int)dtmap.size();
//    int row = stickt->GID(i);
//    double* xi = cnode->xspatial();
//    double* txi = cnode->CoData().txi();
//    double* jump = cnode->Data().jump();
//    double& wgap = cnode->CoData().Getg();
//
//    double utan = 0;
//    double nz = 0;
//
//    for (int dim = 0;dim < Dim();dim++)
//    {
//      utan += txi[dim]*jump[dim];
//      nz += cnode->MoData().n()[dim] * cnode->MoData().lm()[dim];
//    }
//
//    // initialization of nz if nz = wgap = 0
//    if(nz==wgap and wgap==0)
//    {
//      nz=1;
//      cout << "Warning: Initialization of nz to 1" << endl;
//    }
//
//    for (int j=0;j<mapsize-1;++j)
//      if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
//        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivT-map is inconsistent!");
//
//    // Entries on right hand side
//    /**************************************** frcoeff*cn*wgap*ct*utan ***/
//
//    Epetra_SerialDenseVector rhsnode(1);
//    vector<int> lm(1);
//    vector<int> lmowner(1);
//
//    rhsnode(0) = -frcoeff*cn*wgap*ct*utan;
//    lm[0] = cnode->Dofs()[1];
//    lmowner[0] = cnode->Owner();
//
//    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
//
//    // Entries from differentiation with respect to lagrange multipliers
//    /*******************/
//
//    // loop over the dimension
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*cnode->MoData().n()[dim]*ct*utan;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickLMglobal.Assemble(val,row,col);
//    }
//
//    // Entries from differentiation with respect to displacements
//    /************************************************** -tan.(D-Dn-1) ***/
//
//    // we need the nodal entries of the D-matrix and the old one
//    double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
//    double Dold= (cnode->Data().GetDOld()[0])[cnode->Dofs()[0]];
//
//    // loop over all derivative maps (=dimensions)
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*txi[dim]*(D-Dold);
//     //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//     // do not assemble zeros into matrix
//     if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /*************************************************** tan.(M-Mn-1) ***/
//
//    // we need the nodal entries of the M-matrix and the old one
//    vector<map<int,double> > mmap = cnode->MoData().GetM();
//    vector<map<int,double> > mmapold = cnode->Data().GetMOld();
//
//    // create a set of nodes including nodes according to M entries
//    // from current and previous time step
//    set <int> mnodes;
//
//    for (colcurr=mmap[0].begin(); colcurr!=mmap[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    if(mmapold.size()<1)
//    {
//      cout << "GID " << gid << endl;
//      dserror("vector too small");
//    }
//
//    for (colcurr=mmapold[0].begin(); colcurr!=mmapold[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    set<int>::iterator mcurr;
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int dim=0;dim<cnode->NumDof();++dim)
//      {
//        int col = cmnode->Dofs()[dim];
//        double val = -frcoeff*(nz-cn*wgap)*ct*txi[dim]*(mik-mikold);
//        //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//       // do not assemble zeros into matrix
//       if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /***************************************  -DerivT.(D-Dn-1).xs  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*(D-Dold)*xi[j]*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//
//    /***************************************  -DerivT.(M-Mn-1).xm  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      double* mxi = cmnode->xspatial();
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int j=0;j<mapsize;++j)
//      {
//        int k=0;
//
//        // loop over all entries of the current derivative map
//        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//        {
//          int col = colcurr->first;
//          double val = -frcoeff*(nz-cn*wgap)*ct*(mik-mikold)*mxi[j]*colcurr->second;
//          // do not assemble zeros into matrix
//          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//          ++k;
//        }
//
//        if (k!=colsize)
//          dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//      }
//    }
//
//    /**********************************************  -T.DerivD.x  *******/
//
//    // we need the dot product n*x of this node
//    double tdotx = 0.0;
//    for (int dim=0;dim<cnode->NumDof();++dim)
//      tdotx += txi[dim]*xi[dim];
//
//    // prepare assembly
//    map<int,double>& ddmap = cnode->CoData().GetDerivD();
//
//    // loop over all entries of the current derivative map
//    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*tdotx*colcurr->second;
//
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /***********************************************   -T.DerivM.x ******/
//
//    // we need the Lin(M-matrix) entries of this node
//    map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
//    map<int,map<int,double> >::iterator dmcurr;
//
//    // loop over all master nodes in the DerivM-map of the active slave node
//    for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
//    {
//      int gid = dmcurr->first;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      double* mxi = cmnode->xspatial();
//
//      // we need the dot product ns*xm of this node pair
//      double tdotx = 0.0;
//      for (int dim=0;dim<cnode->NumDof();++dim)
//        tdotx += txi[dim]*mxi[dim];
//
//      // compute matrix entry of the current active node / master node pair
//      map<int,double>& thisdmmap = cnode->GetDerivM(gid);
//
//      // loop over all entries of the current derivative map
//      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*tdotx*colcurr->second;
//
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /**************************************-frcoeff*cn*ct*utan*DerivG ***/
//
//    // prepare assembly
//    map<int,double>& dgmap = cnode->CoData().GetDerivG();
//
//    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = +frcoeff*cn*ct*utan*colcurr->second;
//      //cout << "Assemble LinStick: " << row << " " << col << " " << val << endl;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /************************************** -frcoeff*DerivN*z+ct*utan ***/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*cnode->MoData().lm()[j]*ct*utan*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//  }
  return;
}

/*----------------------------------------------------------------------*
|  Assemble matrix LinSlip with tangential+D+M derivatives    mgit 02/09|
*----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
                                         LINALG::SparseMatrix& linslipDISglobal,
                                         Epetra_Vector& linslipRHSglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
  double frbound = IParams().get<double>("FRBOUND");
  double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  bool fulllin = Teuchos::getIntegralValue<int>(IParams(),"FULL_LINEARIZATION");

  // Coulomb Friction
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
#ifdef CONTACTCOMPHUEBER

    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();
      vector<map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      vector<map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       if (Dim()==3)
       {
         for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
       }

      // more information from node
      double* jump = cnode->Data().jump();
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      double& wgap = cnode->CoData().Getg();

      // iterator for maps
      map<int,double>::iterator colcurr;

      // row number of entries
      vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = slipt_->GID(2*i);
        row[1] = slipt_->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinSlip: Dimension not correct");

      // boolean variable if flag "CONTACTFRICTIONLESSFIRST" AND
      // ActiveOld = true
      bool friclessandfirst = false;

      // evaluation of specific components of entries to assemble
      double znor = 0;
      double ztxi = 0;
      double zteta = 0;
      double jumptxi = 0;
      double jumpteta = 0;
      double euclidean = 0;
      for (int i=0;i<Dim();i++)
      {
        znor += n[i]*z[i];
        ztxi += txi[i]*z[i];
        zteta += teta[i]*z[i];
        jumptxi += txi[i]*jump[i];
        jumpteta += teta[i]*jump[i];
      }

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      vector<double> sum1 (Dim()-1,0);
      sum1[0] = ztxi+ct*jumptxi;
      if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
      if (Dim()==2) euclidean = abs(sum1[0]);
      if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

      // check of dimensions
      if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      // check of euclidean norm
      if (euclidean==0.0)
        dserror ("ERROR: AssemblelinSlip: Euclidean norm is zero");

#ifdef CONTACTFRICTIONLESSFIRST

      // in the case of frictionless contact for nodes just coming into
      // contact, the frictionless contact condition is applied.
      if (cnode->Data().ActiveOld()==false)
      {
        friclessandfirst=true;
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];
          double valteta = 0;
          if (Dim()==3) valteta = teta[dim];

          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);

        }
        if(fulllin)
        {
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (colcurr->second)*z[dim];
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            }

            if(Dim()==3)
            {
              for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
              {
                int col = colcurr->first;
                double valteta = (colcurr->second)*z[dim];
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
        }
      }
#endif

      // this is not evaluated if "FRICTIONLESSFIRST" is flaged on AND the node
      // is just coming into contact
      if(friclessandfirst==false)
      {
        /******************************************************************/
        // calculation of matrix entries of the linearized slip condition
        /******************************************************************/
        // 1) Entries from differentiation with respect to lagrange multipliers
        // 2) Entries on right hand side
        // 3) Entries from differentiation with respect to displacements

        // 1) Entries from differentiation with respect to lagrange multipliers
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          double valtxi = 0;
          int col = cnode->Dofs()[dim];
          double valtxi0 = euclidean*txi[dim];
          double valtxi1 = ((ztxi+ct*jumptxi)/euclidean*ztxi)*txi[dim];
          double valtxi3 = (zteta+ct*jumpteta)/euclidean*ztxi*teta[dim];
          double valtxi2 = -frcoeff*(znor-cn*wgap)*txi[dim]-frcoeff*(ztxi+ct*jumptxi)*n[dim];
          valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

          double valteta = 0;
          if (Dim()==3)
          {
            double valteta0 = euclidean*teta[dim];
            double valteta1 = ((ztxi+ct*jumptxi)/euclidean*zteta)*txi[dim];
            double valteta3 = (zteta+ct*jumpteta)/euclidean*zteta*teta[dim];
            double valteta2 = -frcoeff*(znor-cn*wgap)*teta[dim]-frcoeff*(zteta+ct*jumpteta)*n[dim];
            valteta = valteta0 + valteta1 + valteta2 + valteta3;
          }

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
        }

        // 2) Entries on right hand side
        /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

        double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
        double valuetxi2 = +euclidean*ztxi;
        double valuetxi3 = (ztxi+ct*jumptxi)/euclidean*ztxi*ztxi;
        double valuetxi4 = (zteta+ct*jumpteta)/euclidean*zteta*ztxi;
        double valuetxi5 = -(frcoeff*(znor-cn*wgap))*ztxi-(frcoeff*znor)*(ztxi+ct*jumptxi);

        Epetra_SerialDenseVector rhsnode(Dim()-1);
        vector<int> lm(Dim()-1);
        vector<int> lmowner(Dim()-1);

        rhsnode(0) = (valuetxi1+valuetxi2+valuetxi3+valuetxi4+valuetxi5);
        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        if(Dim()==3)
        {
          double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
          double valueteta2 = +euclidean*zteta;
          double valueteta3 = (ztxi+ct*jumptxi)/euclidean*ztxi*zteta;
          double valueteta4 = (zteta+ct*jumpteta)/euclidean*zteta*zteta;
          double valueteta5 = -(frcoeff*(znor-cn*wgap))*zteta-(frcoeff*znor)*(zteta+ct*jumpteta);

          rhsnode(1) = (valueteta1+valueteta2+valueteta3+valueteta4+valueteta5);
          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }

        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

        // 3) Entries from differentiation with respect to displacements
        /*** 01  *********** -Deriv(euclidean).ct.tangent.(D-Dn-1)*ztan ***/

        // we need the nodal entries of the D-matrix and the old one
        vector<map<int,double> > dmapold=cnode->Data().GetDOld();
        if(dmapold.size()<1)
            dserror("ERROR: AssembleLinSlip: No dmap-entries form previous time step");

        double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
        double Dold= dmapold[0][cnode->Dofs()[0]];

        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          double valtxi1=0;
          double valtxi2=0;
          double valteta1=0;
          double valteta2=0;

          int col = cnode->Dofs()[dim];
          valtxi1 = (ztxi+ct*jumptxi)/euclidean*(-1)*ct*txi[dim]*(D-Dold)*ztxi;

          if(Dim()==3)
          {
            valteta1 = (ztxi+ct*jumptxi)/euclidean*(-1)*ct*txi[dim]*(D-Dold)*zteta;
            valtxi2 = (zteta+ct*jumpteta)/euclidean*(-1)*ct*teta[dim]*(D-Dold)*ztxi;
            valteta2 = (zteta+ct*jumpteta)/euclidean*(-1)*ct*teta[dim]*(D-Dold)*zteta;
          }

          // do not assemble zeros into matrix
          if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
          if (Dim()==3)
          {
            if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
            if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
            if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
          }
         }

        /*** 02  *********** -Deriv(euclidean).ct.tangent.(M-Mn-1)*ztan ***/

        // we need the nodal entries of the M-matrix and the old one
        vector<map<int,double> > mmap = cnode->MoData().GetM();
        vector<map<int,double> > mmapold = cnode->Data().GetMOld();

        // map from previous time step must have an entry
        if(mmapold.size()<1)
           dserror("ERROR: AssembleLinStick: No mmap-entries form previous time step");

        // create a set of nodes including nodes according to M entries
         // from current and previous time step
         set <int> mnodes;

         // iterator
         set<int>::iterator mcurr;

         set <int> mnodescurrent = cnode->Data().GetMNodes();
        set <int> mnodesold = cnode->Data().GetMNodesOld();

        for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
          mnodes.insert(*mcurr);

        for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
          mnodes.insert(*mcurr);

        // loop over all master nodes (find adjacent ones to this slip node)
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          const int* mdofs = cmnode->Dofs();

          double mik = (mmap[0])[mdofs[0]];
          double mikold = (mmapold[0])[mdofs[0]];

          double valtxi1=0;
          double valtxi2=0;
          double valteta1=0;
          double valteta2=0;

          // compute linstick-matrix entry of the current active node / master node pair
          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            valtxi1 = (ztxi+ct*jumptxi)/euclidean*(+1)*ct*txi[dim]*(mik-mikold)*ztxi;

            if(Dim()==3)
            {
              valteta1 = (ztxi+ct*jumptxi)/euclidean*(+1)*ct*txi[dim]*(mik-mikold)*zteta;
              valtxi2 = (zteta+ct*jumpteta)/euclidean*(+1)*ct*teta[dim]*(mik-mikold)*ztxi;
              valteta2 = (zteta+ct*jumpteta)/euclidean*(+1)*ct*teta[dim]*(mik-mikold)*zteta;
            }

            // do not assemble zeros into matrix
            if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
            if (Dim()==3)
            {
              if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
              if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
              if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
            }
          }
        }

        /*** 03 ********************** frcoeff*znor*ct*tangent.(D-Dn-1) ***/
        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          double valtxi=0;
          double valteta=0;
          int col = cnode->Dofs()[dim];
          valtxi = (frcoeff*(znor-cn*wgap))*ct*txi[dim]*(D-Dold);

          if (Dim()==3)
            valteta = (frcoeff*(znor-cn*wgap))*ct*teta[dim]*(D-Dold);

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }

        /*** 04 ********************** frcoeff*znor*ct*tangent.(M-Mn-1) ***/

        // loop over all master nodes
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          double valtxi=0;
          double valteta=0;

          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          const int* mdofs = cmnode->Dofs();

          double mik = (mmap[0])[mdofs[0]];
          double mikold = (mmapold[0])[mdofs[0]];

          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            valtxi = (frcoeff*(znor-cn*wgap))*(-1)*ct*txi[dim]*(mik-mikold);

            if (Dim()==3)
              valteta = (frcoeff*(znor-cn*wgap))*(-1)*ct*teta[dim]*(mik-mikold);

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if(Dim()==3)
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
           }
         }

        // remaining terms only in case of full linearization
        if(fulllin)
        {
          /*** 1 ********************************* euclidean.deriv(T).z ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean*(colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if (Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = euclidean*(colcurr->second)*z[j];

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*zteta;

             // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (Dim()==3)
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // 3D loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*zteta;

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/

          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 4 ******************* deriv(euclidean).tan.deriv(D).ztan ***/
          // we need the dot product t*x of this node
          double txidotx = 0.0;
          double tetadotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            txidotx += txi[dim]*xi[dim];
            tetadotx += teta[dim]*xi[dim];
          }

          // prepare assembly
          map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
          for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi1 = (-1)*(ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*ztxi;
            double valteta1 = (-1)*(ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*zteta;
            double valtxi2 = (-1)*(zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*ztxi;
            double valteta2 = (-1)*(zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*zteta;

            // do not assemble zeros into matrix
            if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
            if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
            if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
            if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
          }

          /*** 5 ******************* deriv(euclidean).tan.deriv(M).ztan ***/
          // we need the Lin(M-matrix) entries of this node
          map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
          map<int,map<int,double> >::iterator dmcurr;

          // loop over all master nodes in the DerivM-map of the active slave node
          for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
          {
            int gid = dmcurr->first;
            DRT::Node* mnode = idiscret_->gNode(gid);
            if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
            FriNode* cmnode = static_cast<FriNode*>(mnode);
            double* mxi = cmnode->xspatial();

            // we need the dot product ns*xm of this node pair
            double txidotx = 0.0;
            double tetadotx = 0.0;
            for (int dim=0;dim<cnode->NumDof();++dim)
            {
              txidotx += txi[dim]*mxi[dim];
              tetadotx += teta[dim]*mxi[dim];
            }

            // compute entry of the current active node / master node pair
            map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

            // loop over all entries of the current derivative map
            for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*ztxi;
              double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*zteta;
              double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*ztxi;
              double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*zteta;

              // do not assemble zeros into matrix
              if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
              if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
              if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
              if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
            }
          }

          /*** 6 **************************** (frcoeff*znor).deriv(T).z ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

                // do not assemble zeros into matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 7 ************************* (frcoeff*znor).deriv(T).jump ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 8 ********************* (frcoeff*znor).ct.tan.deriv(D).x ***/
          // we need the dot product t*x of this node
          txidotx = 0.0;
          tetadotx = 0.0;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
           txidotx += txi[dim]*xi[dim];
           tetadotx += teta[dim]*xi[dim];
          }

          // loop over all entries of the current derivative map
          for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (-1)*(-1)*(frcoeff*(znor-cn*wgap))*ct*txidotx*colcurr->second;
            double valteta = (-1)*(-1)*(frcoeff*(znor-cn*wgap))*ct*tetadotx*colcurr->second;

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
             if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

          /*** 9 ********************* (frcoeff*znor).ct.tan.deriv(M).x ***/

          // loop over all master nodes in the DerivM-map of the active slave node
          for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
          {
            int gid = dmcurr->first;
            DRT::Node* mnode = idiscret_->gNode(gid);
            if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
            FriNode* cmnode = static_cast<FriNode*>(mnode);
            double* mxi = cmnode->xspatial();

            // we need the dot product t*xm of this node pair
            double txidotx = 0.0;
            double tetadotx = 0.0;
            for (int dim=0;dim<cnode->NumDof();++dim)
            {
              txidotx += txi[dim]*mxi[dim];
              tetadotx += teta[dim]*mxi[dim];
            }

            // compute entry of the current active node / master node pair
            map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

             // loop over all entries of the current derivative map
             for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
             {
               int col = colcurr->first;
               double valtxi = (-1)*(frcoeff*(znor-cn*wgap))*ct*txidotx*colcurr->second;
               double valteta = (-1)*(frcoeff*(znor-cn*wgap))*ct*tetadotx*colcurr->second;

               // do not assemble zeros into matrix
               if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
               if (Dim()==3)
               {
                 if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
               }
            }
          }

          /*** 10 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map
            for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (-1)*(ztxi+ct*jumptxi)*frcoeff*(colcurr->second)*z[j];
              double valteta = (-1)*(zteta+ct*jumpteta)*frcoeff*(colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

          /*** 11 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
          // prepare assembly
          map<int,double>& dgmap = cnode->CoData().GetDerivG();

          // loop over all entries of the current derivative map
          for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
            double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }
        } // if fullin
      } // if (frictionlessandfirst == false)
    } // loop over all slip nodes of the interface
#else
    
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->Data().jump();
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();
      
      cout << dtmap.size() << endl;
      cout << dtmap[0].size() << endl;
      cout << dtmap[1].size() << endl;
      cout << dtmap[2].size() << endl;
       
      exit(0);

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double znor    = n[0]*z[0] + n[1]*z[1];
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frcoff*znor).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frcoeff*znor)*txi[dim]-frcoeff*(ztan+ct*jumptan)*n[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1 = -(abs(ztan+ct*jumptan))*ztan+(frcoeff*znor)*(ztan+ct*jumptan);
      double value2 = -(frcoeff*znor)*ztan-(frcoeff*znor)*(ztan+ct*jumptan);
      double value3 = +sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      rhsnode(0) = (value1+value2+value3);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->Data().GetDOld()[0])[cnode->Dofs()[0]];

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif


       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->Data().GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->Data().GetMNodes();
      set <int> mnodesold = cnode->Data().GetMNodesOld();

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /********************************* frcoeff*znor*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (frcoeff*znor)*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -frcoeff*znor*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = (frcoeff*znor)*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->Data().ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

          // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /*********************************** -(frcoeff*znor).DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** (frcoeff*znor).ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** +(frcoeff*znor).ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /***************************** -(frcoeff*znor).ct.T.DerivM.x ***/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           FriNode* cmnode = static_cast<FriNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /***************************** -frcoeff*DerivN.z(ztan+ct*utan) ***/

       // loop over all derivative maps (=dimensions)
       for (int j=0;j<mapsize;++j)
       {
         int k=0;

         // loop over all entries of the current derivative map
         for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(ztan+ct*jumptan)*frcoeff*(colcurr->second)*z[j];
           //cout << "10 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into s matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
           ++k;
         }

         if (k!=colsize)
           dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }
      } // if fullin
    }
#endif
  } // Coulomb friction

  // Tresca Friction
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->Data().jump();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frbound)*txi[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1= -(abs(ztan+ct*jumptan))*ztan+frbound*(ztan+ct*jumptan);
      double value2= -frbound*ztan+sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);
      rhsnode(0) = (value1+value2);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->Data().GetDOld()[0])[cnode->Dofs()[0]];

      if (abs(Dold)<0.0001)
        dserror ("Error:No entry for Dold");

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->Data().GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->Data().GetMNodes();
      set <int> mnodesold = cnode->Data().GetMNodesOld();

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /************************************** frbound*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = frbound*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

        // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = frbound*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->Data().ActiveOld()==false) val = 0;
#endif
          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
            if (cnode->Data().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->Data().ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*frbound*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /********************************  -frbound.ct.T.DerivM.x ******/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           FriNode* cmnode = static_cast<FriNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*frbound*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->Data().ActiveOld()==false) val = 0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      } // if fullin
    }
  }// Tresca friction

  return;
}

/*----------------------------------------------------------------------*
 |  initialize active set (nodes / dofs)                      popp 03/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::InitializeActiveSet()
{
  // define local variables
  int countnodes = 0;
  int countdofs = 0;
  vector<int> mynodegids(snoderowmap_->NumMyElements());
  vector<int> mydofgids(sdofrowmap_->NumMyElements());
  
  // define local variables for slip maps
  vector<int> myslipnodegids(0);
  vector<int> myslipdofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (cnode->IsInitActive())
    {
      // FULL: all slave nodes are assumed to be active
      cnode->Active()=true;
      mynodegids[countnodes] = cnode->Id();
      ++countnodes;

      for (int j=0;j<numdof;++j)
      {
        mydofgids[countdofs] = cnode->Dofs()[j];
        ++countdofs;
      }
    }
  }

  // resize the temporary vectors
  mynodegids.resize(countnodes);
  mydofgids.resize(countdofs);

  // communicate countnodes and countdofs among procs
  int gcountnodes, gcountdofs;
  Comm().SumAll(&countnodes,&gcountnodes,1);
  Comm().SumAll(&countdofs,&gcountdofs,1);

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(gcountdofs,countdofs,&mydofgids[0],0,Comm()));

  // create an empty slip node map and slip dof map
  // for the first time step (t=0) all nodes are stick nodes
  if(friction_)
  {
    slipnodes_   = rcp(new Epetra_Map(0,0,Comm()));
    slipdofs_    = rcp(new Epetra_Map(0,0,Comm()));
  }
  
  // split active dofs into Ndofs, Tdofs and slipTdofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::BuildActiveSet()
{
  // define local variables
  int countnodes = 0;
  int countdofs = 0;
  vector<int> mynodegids(snoderowmap_->NumMyElements());
  vector<int> mydofgids(sdofrowmap_->NumMyElements());
  
  // FRICTION - define local variables
  int countslipnodes = 0;
  int countslipdofs = 0;
  vector<int> myslipnodegids(snoderowmap_->NumMyElements());
  vector<int> myslipdofgids(sdofrowmap_->NumMyElements());
  
  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add node / dofs to temporary map IF active
    if (cnode->Active())
    {
      mynodegids[countnodes] = cnode->Id();
      ++countnodes;

      for (int j=0;j<numdof;++j)
      {
        mydofgids[countdofs] = cnode->Dofs()[j];
        ++countdofs;
      }
    }
        
    // add node / dofs to temporary map IF slip
    if (friction_)
    {
      if (static_cast<FriNode*>(cnode)->Data().Slip())
      {
        myslipnodegids[countslipnodes] = cnode->Id();
        ++countslipnodes;

        for (int j=0;j<numdof;++j)
        {
          myslipdofgids[countslipdofs] = cnode->Dofs()[j];
          ++countslipdofs;
        }
      }  
    }
  }

  // resize the temporary vectors
  mynodegids.resize(countnodes);
  mydofgids.resize(countdofs);

  // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
  int gcountnodes, gcountdofs, gcountslipnodes,gcountslipdofs;
  Comm().SumAll(&countnodes,&gcountnodes,1);
  Comm().SumAll(&countdofs,&gcountdofs,1);

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(gcountdofs,countdofs,&mydofgids[0],0,Comm()));
  
  // the last three steps for frictional contact
  if (friction_)
  {
    myslipnodegids.resize(countslipnodes);
    myslipdofgids.resize(countslipdofs);

    Comm().SumAll(&countslipnodes,&gcountslipnodes,1);
    Comm().SumAll(&countslipdofs,&gcountslipdofs,1);

    slipnodes_ = rcp(new Epetra_Map(gcountslipnodes,countslipnodes,&myslipnodegids[0],0,Comm()));
    slipdofs_  = rcp(new Epetra_Map(gcountslipdofs,countslipdofs,&myslipdofgids[0],0,Comm()));
  }
  
  // split active dofs into Ndofs and Tdofs and slip dofs in Tslipdofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs          popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::SplitActiveDofs()
{
  // get out of here if active set is empty
  if (activenodes_==null)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements()==0)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countN=0;
  int countT=0;
  vector<int> myNgids(activenodes_->NumMyElements());
  vector<int> myTgids((Dim()-1)*activenodes_->NumMyElements());

  // dimension check
  double dimcheck =(activedofs_->NumGlobalElements())/(activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j=1;j<numdof;++j)
    {
      myTgids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  Comm().SumAll(&countN,&gcountN,1);
  Comm().SumAll(&countT,&gcountT,1);

  // check global dimensions
  if ((gcountN+gcountT)!=activedofs_->NumGlobalElements())
    dserror("ERROR: SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));
  activet_ = rcp(new Epetra_Map(gcountT,countT,&myTgids[0],0,Comm()));
  
  // *******************************************************************
  // FRICTION - EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if there is no friction
  if(friction_==false)
    return true; 
  
  // get out of here if slip set is empty
  if (slipnodes_==null)
  {
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  if (slipnodes_->NumGlobalElements()==0)
  {
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countslipT=0;
  vector<int> myslipTgids((Dim()-1)*slipnodes_->NumMyElements());

  // dimension check
  dimcheck =(slipdofs_->NumGlobalElements())/(slipnodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add dofs to slipTmap
    for (int j=1;j<numdof;++j)
    {
      myslipTgids[countslipT] = cnode->Dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  Comm().SumAll(&countslipT,&gcountslipT,1);

  // create Tslipmap objects
  slipt_   = rcp(new Epetra_Map(gcountslipT,countslipT,&myslipTgids[0],0,Comm()));

  return true;
}

#endif  // #ifdef CCADISCRET
