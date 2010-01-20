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
selfcontact_(selfcontact)
{
  // empty constructor
  
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

    //reset nodal normal and tangents and jumps
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((node->GetDerivTeta()).size());++j)
      (node->GetDerivTeta())[j].clear();
    (node->GetDerivTeta()).resize(0);

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();

    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((node->GetDerivZ()).size()); ++j)
      (node->GetDerivZ())[j].clear();
    (node->GetDerivZ()).resize(0);

    // reset nodal weighted gap and derivative
    node->Getg() = 1.0e12;
    (node->GetDerivG()).clear();

    // reset feasible projection status
    node->HasProj() = false;
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

  return true;
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
    double gap = cnode->Getg();

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    //double modgap = cnode->Getg();
    //if (abs(modgap) < 1.0e-10) modgap=0.0;

    double kappa = cnode->Kappa();

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->lmuzawa()[k]*cnode->n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitely to constant
    // and corresponding derivatives to zero

    for( int j=0;j<dim;++j)
      cnode->lm()[j] = i*j;

    cnode->GetDerivZ().clear();

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

      double* normal = cnode->n();

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      map<int,double>& derivg = cnode->GetDerivG();
      map<int,double>::iterator gcurr;

      // contribution of derivative of normal
      vector<map<int,double> >& derivn = cnode->GetDerivN();
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
      for( int j=0;j<dim;++j) cnode->lm()[j] = 0;

      // clear derivz
      cnode->GetDerivZ().clear();

    } // Macauley-Bracket
  } // loop over slave nodes

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
    vector<map<int,double> >& derivz = cnode->GetDerivZ();

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
      double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(1,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->txi()[j];
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
      double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(2,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->txi()[j];
        Tnode(1,j) = cnode->teta()[j];
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
    map<int,double>& dgmap = cnode->GetDerivG();
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
      vector<map<int,double> >& dtmap = cnode->GetDerivTxi();
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
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
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
      vector<map<int,double> >& dtximap = cnode->GetDerivTxi();
      vector<map<int,double> >& dtetamap = cnode->GetDerivTeta();
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
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
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
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
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
    map<int,map<int,double> >& dderiv = cnode->GetDerivD();

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

        double* lm = csnode->lm();

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
    map<int,double>& dderiv = cnode->GetDerivD();
    map<int,double>::iterator colcurr;

    // current Lagrange multipliers
    double* lm = cnode->lm();

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

    // Mortar matrix M derivatives
    map<int,map<int,double> >& mderiv = cnode->GetDerivM();
    map<int,double>::iterator colcurr;
    map<int,map<int,double> >::iterator mcurr;

    // current Lagrange multipliers
    double* lm = cnode->lm();

    // inter-proc. communication
    // (we want to keep all procs around, although at the moment only
    // the cnode owning proc does relevant work!)
    int mastersize = 0;
    if (Comm().MyPID()==cnode->Owner())
    {
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
    for (int j=0;j<mastersize;++j)
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
      int mnumdof = cmnode->NumDof();

      // Mortar matrix M derivatives
      map<int,double>& thismderiv = cnode->GetDerivM()[mgid];
      int mapsize = 0;
      if (Comm().MyPID()==cnode->Owner())
        mapsize = (int)(thismderiv.size());
      lComm()->Broadcast(&mapsize,1,procmap_[cnode->Owner()]);

      // loop over all master node dofs
      for (int k=0;k<mnumdof;++k)
      {
        int row = cmnode->Dofs()[k];

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
            val = lm[k] * (colcurr->second);
            ++colcurr;
          }
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
    if (cnode->Getg()!=0.0)
    {
      double gap = cnode->Getg();

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
        cnode->Getg()=gap;
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
 |  initialize active set (nodes / dofs)                      popp 03/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::InitializeActiveSet()
{
  // define local variables
  int countnodes = 0;
  int countdofs = 0;
  vector<int> mynodegids(snoderowmap_->NumMyElements());
  vector<int> mydofgids(sdofrowmap_->NumMyElements());

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
  }

  // resize the temporary vectors
  mynodegids.resize(countnodes);
  mydofgids.resize(countdofs);

  // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
  int gcountnodes, gcountdofs;
  Comm().SumAll(&countnodes,&gcountnodes,1);
  Comm().SumAll(&countdofs,&gcountdofs,1);

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(gcountdofs,countdofs,&mydofgids[0],0,Comm()));

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
    return true;
  }

  else if (activenodes_->NumGlobalElements()==0)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
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

  return true;
}

#endif  // #ifdef CCADISCRET
