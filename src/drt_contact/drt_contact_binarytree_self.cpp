/*!----------------------------------------------------------------------
\file drt_contact_binarytree_self.cpp
\brief A class for performing self contact search in 2D / 3D based
       on binary search trees and dual graphs

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

#include "drt_contact_binarytree_self.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

                     
/*----------------------------------------------------------------------*
 |  ctor BinaryTreeNode for self contact (public)             popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::BinaryTreeSelfNode::BinaryTreeSelfNode(
                       BinaryTreeSelfNodeType type,
                       DRT::Discretization& discret,
                       RCP<BinaryTreeSelfNode> parent,
                       vector<int> elelist,
                       const Epetra_SerialDenseMatrix& dopnormals,
                       const Epetra_SerialDenseMatrix& samplevectors,
                       const int& kdop, const int& dim, const int& nvectors,
                       const int layer,
                       vector<vector<RCP<BinaryTreeSelfNode> > > & treenodes) :
type_(type),
idiscret_(discret),
parent_(parent),
elelist_(elelist),
dopnormals_(dopnormals),
samplevectors_(samplevectors),
kdop_(kdop),
dim_(dim),
nvectors_(nvectors),
layer_(layer),
treenodes_(treenodes)
{
  // reshape slabs matrix
  if (dim_==2)      slabs_.Reshape(kdop_/2,2);
  else if (dim_==3) slabs_.Reshape(kdop_/2,2);
  else              dserror("ERROR: Problem dimension must be 2D or 3D!");
  
  return;
}

/*----------------------------------------------------------------------*
 | complete the tree storage in a top down way (public)       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::CompleteTree(int layer,double& enlarge)
{ 
  // calculate bounding volume
  CalculateSlabsDop(true);
  EnlargeGeometry(enlarge);
  
  // check root node (layer 0) for sanity
  if (layer == 0)
  {
    SetLayer(0);
    if (type_ == SELF_INNER) treenodes_[layer].push_back(rcp(this));
    else dserror("ERROR: root must be inner node in treenodes scheme");
  }
  
  // build treenode storage recursively
  if (type_==SELF_INNER)
  {
    leftchild_->SetLayer(layer_+1);
    rightchild_->SetLayer(layer_+1);
    
    // if map of treenodes does not have enogh rows-->resize!
    if ((int)(treenodes_.size()) <= (layer_+1))
      treenodes_.resize((layer_+2));  
    
    // put new pointers to children into map
    treenodes_[(layer_+1)].push_back(leftchild_); 
    treenodes_[(layer_+1)].push_back(rightchild_);
         
    rightchild_->CompleteTree(layer_+1,enlarge);
    leftchild_->CompleteTree(layer_+1,enlarge);
  }
  
  // do nothing if arrived at leaf level
 
 return;
}

/*----------------------------------------------------------------------*
 | Calculate booleans for qualified sample vectors (public)   popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::CalculateQualifiedVectors()
{
  if (type_!=SELF_LEAF) dserror("ERROR: Calculate qual. vec. called for non-leaf node!");

  // resize qualified vectors
  qualifiedvectors_.resize(nvectors_);
  
  // we first need the element center:
  // for line2, line3, quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  CElement* celement= static_cast<CElement*>(idiscret_.gElement(elelist_[0]));
  double loccenter[2];

  DRT::Element::DiscretizationType dt = celement->Shape();
  if (dt==CElement::tri3 || dt==CElement::tri6)
  {
    loccenter[0] = 1.0/3;
    loccenter[1] = 1.0/3;
  }
  else if (dt==CElement::line2 || dt==CElement::line3 || dt==CElement::quad4 ||
           dt==CElement::quad8 || dt==CElement::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else dserror("ERROR: CalculateQualifiedVectors called for unknown element type");
    
  // now get the element center normal
  double normal[3] = {0.0, 0.0, 0.0};
  celement->ComputeUnitNormalAtXi(loccenter,normal);
  
  // check normal against sample vectors
  for (int i=0;i<(int)qualifiedvectors_.size();++i)
  {
    double scalar = (double)samplevectors_(i,0)*normal[0]
                  + (double)samplevectors_(i,1)*normal[1]
                  + (double)samplevectors_(i,2)*normal[2];
    if (scalar > 0) qualifiedvectors_[i] = true;
    else            qualifiedvectors_[i] = false;
  }
  
  return;  
}

/*----------------------------------------------------------------------*
 | Update qualified sample vectors for inner node (public)    popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::UpdateQualifiedVectorsBottomUp()
{
  if(type_==SELF_LEAF) dserror("ERROR: Update qual. vec. called for leaf node!");

  // calculate the qualified vectors (= valid samplevectors) of an inner
  // treenode by comparing the qualified vectors of the children
  qualifiedvectors_.resize(nvectors_);
  
  for (int i=0;i<(int)qualifiedvectors_.size();++i)
    qualifiedvectors_.at(i) = ( (rightchild_->QualifiedVectors()).at(i)
                             && (leftchild_->QualifiedVectors()).at(i));
  
  return;
}

/*----------------------------------------------------------------------*
 | Update endnodes of one tree node (only 2D) (public)        popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::UpdateEndnodes()
{
  if(type_==SELF_LEAF) dserror("ERROR: Update endnodes. called for leaf node!");
  
  // reset endnodes
  endnodes_.clear();
 
  // find out which nodes the children have in commen, save others as endnodes
  if (leftchild_->endnodes_[0] == rightchild_->endnodes_[0]
   && leftchild_->endnodes_[1] != rightchild_->endnodes_[1])
  {
    endnodes_.push_back(rightchild_->endnodes_[1]);
    endnodes_.push_back(leftchild_->endnodes_[1]);
  }
  
  else if (leftchild_->endnodes_[1] == rightchild_->endnodes_[1] 
        && leftchild_->endnodes_[0] != rightchild_->endnodes_[0])
  {
    endnodes_.push_back(rightchild_->endnodes_[0]);
    endnodes_.push_back(leftchild_->endnodes_[0]);
  }
  
  else if (leftchild_->endnodes_[0] == rightchild_->endnodes_[1]
        && leftchild_->endnodes_[1] != rightchild_->endnodes_[0])
  {
    endnodes_.push_back(rightchild_->endnodes_[0]);
    endnodes_.push_back(leftchild_->endnodes_[1]);
  }
  
  else if (leftchild_->endnodes_[1] == rightchild_->endnodes_[0] 
        && leftchild_->endnodes_[0] != rightchild_->endnodes_[1] )
  {
    endnodes_.push_back(rightchild_->endnodes_[1]);
    endnodes_.push_back(leftchild_->endnodes_[0]);
  }
  
  else // the treenode is a closed surface (ring) -> no endnodes
  {
    endnodes_.push_back(-1);
    endnodes_.push_back(-1);
  }
    
  return;
}

/*----------------------------------------------------------------------*
 | Calculate slabs of DOP out of node postions (public)       popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::CalculateSlabsDop(bool isinit)
{
  // initialize slabs
  for (int j=0; j<kdop_/2; j++)
  {
    slabs_(j,0) =  1.0e12;
    slabs_(j,1) = -1.0e12;
  } 

  // calculate slabs for every element
  for (int i=0; i<(int)elelist_.size();++i)
  {
    int gid = Elelist()[i];
    DRT::Element* element= idiscret_.gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n",gid);
    DRT::Node** nodes = element->Nodes();
    if (!nodes) dserror("ERROR: Null pointer!");
    
    // calculate slabs for every node on every element
    for (int k=0;k<element->NumNode();++k)
    {
      CNode* cnode=static_cast<CNode*>(nodes[k]);
      if (!cnode) dserror("ERROR: Null pointer!");
      
      // decide which position is relevant (initial or current)
      double pos[3] = {0.0, 0.0, 0.0};
      for (int j=0;j<dim_;++j)
      {
        if (isinit) pos[j] = cnode->X()[j];
        else        pos[j] = cnode->xspatial()[j];
      }

      // calculate slabs
      for(int j=0;j<kdop_/2;++j)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        double num = dopnormals_(j,0)*pos[0]
                   + dopnormals_(j,1)*pos[1]
                   + dopnormals_(j,2)*pos[2];
        double denom = sqrt((dopnormals_(j,0)*dopnormals_(j,0))
                           +(dopnormals_(j,1)*dopnormals_(j,1))
                           +(dopnormals_(j,2)*dopnormals_(j,2)));
        double dcurrent = num/denom;
        
        if (dcurrent > slabs_(j,1)) slabs_(j,1) = dcurrent;
        if (dcurrent < slabs_(j,0)) slabs_(j,0) = dcurrent;
      }
      
      // if update for contactsearch --> add auxiliary positions
      if ((!isinit) && (type_==SELF_INNER || type_==SELF_LEAF))
      {
        double auxpos[3] = {0.0, 0.0, 0.0};
        double scalar = 0.0;
        for (int j=0;j<dim_;++j)
          scalar = scalar + (cnode->X()[j]+cnode->uold()[j]-cnode->xspatial()[j])*cnode->n()[j];
        for (int j=0;j<dim_;++j)
          auxpos[j] = cnode->xspatial()[j] + scalar*cnode->n()[j];

        for (int j=0;j<kdop_/2;++j)
        {
          //= ax+by+cz=d/sqrt(aa+bb+cc)
          double num = dopnormals_(j,0)*auxpos[0]
                     + dopnormals_(j,1)*auxpos[1]
                     + dopnormals_(j,2)*auxpos[2];
          double denom = sqrt((dopnormals_(j,0)*dopnormals_(j,0))
                             +(dopnormals_(j,1)*dopnormals_(j,1))
                             +(dopnormals_(j,2)*dopnormals_(j,2)));
          double dcurrent = num/denom;
                  
          if (dcurrent > slabs_(j,1)) slabs_(j,1) = dcurrent;
          if (dcurrent < slabs_(j,0)) slabs_(j,0) = dcurrent;
        }
      }       
    }
  }

  //Prints slabs to std::cout
  //PrintSlabs();
  
  return;
}

/*----------------------------------------------------------------------*
 | Update slabs bottom-up (public)                            popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::UpdateSlabsBottomUp(double & enlarge)
{
  // if current treenode is inner node
  if (type_==SELF_INNER)
  {
    for (int k=0;k<kdop_/2;++k)
    {
      // for minimum
      if (leftchild_->Slabs()(k,0) <= rightchild_->Slabs()(k,0))
        slabs_(k,0) = leftchild_->Slabs()(k,0);
      else 
        slabs_(k,0) = rightchild_->Slabs()(k,0);
        
      // for maximum
      if (leftchild_->Slabs()(k,1) >= rightchild_->Slabs()(k,1))
        slabs_(k,1) = leftchild_->Slabs()(k,1);
      else
        slabs_(k,1) = rightchild_->Slabs()(k,1);    
    }
  }
  
  // if current treenode is leaf node
  if (type_==SELF_LEAF)
  {
    // initialize slabs
    for (int j=0;j<kdop_/2;++j)
    {
      slabs_(j,0) =  1.0e12;
      slabs_(j,1) = -1.0e12;
    }
    
    int gid = Elelist()[0];
    DRT::Element* element= idiscret_.gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n",gid);
    DRT::Node** nodes = element->Nodes();
    if (!nodes) dserror("ERROR: Null pointer!");
    
    // update slabs for every node
    for (int k=0;k<element->NumNode();++k)
    {
      CNode* cnode=static_cast<CNode*>(nodes[k]);
      if (!cnode) dserror("ERROR: Null pointer!");
      
      // decide which position is relevant (initial or current)
      double pos[3] = {0.0, 0.0, 0.0};
      for (int j=0;j<dim_;++j) pos[j] = cnode->xspatial()[j];

      // calculate slabs
      for (int j=0;j<kdop_/2;++j)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        double num = dopnormals_(j,0)*pos[0]
                   + dopnormals_(j,1)*pos[1]
                   + dopnormals_(j,2)*pos[2];
        double denom = sqrt((dopnormals_(j,0)*dopnormals_(j,0))
                           +(dopnormals_(j,1)*dopnormals_(j,1))
                           +(dopnormals_(j,2)*dopnormals_(j,2)));
        double dcurrent = num/denom;
        
        if (dcurrent > slabs_(j,1)) slabs_(j,1) = dcurrent;
        if (dcurrent < slabs_(j,0)) slabs_(j,0) = dcurrent;
      }
      
      // enlarge slabs with aux. position   
      double auxpos[3] = {0.0, 0.0, 0.0};
      double scalar = 0.0;
      for (int j=0;j<dim_;++j)
        scalar = scalar + (cnode->X()[j]+cnode->uold()[j]-cnode->xspatial()[j])*cnode->n()[j];
      for (int j=0;j<dim_;++j)
        auxpos[j] = cnode->xspatial()[j] + scalar*cnode->n()[j];

      for (int j=0;j<kdop_/2;++j)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        double num = dopnormals_(j,0)*auxpos[0]
                   + dopnormals_(j,1)*auxpos[1]
                   + dopnormals_(j,2)*auxpos[2];
        double denom = sqrt((dopnormals_(j,0)*dopnormals_(j,0))
                           +(dopnormals_(j,1)*dopnormals_(j,1))
                           +(dopnormals_(j,2)*dopnormals_(j,2)));
        double dcurrent = num/denom;
                
        if (dcurrent > slabs_(j,1)) slabs_(j,1) = dcurrent;
        if (dcurrent < slabs_(j,0)) slabs_(j,0) = dcurrent;
      }
    }
   
    for (int i=0;i<kdop_/2;++i)
    {
      slabs_(i,0) = slabs_(i,0) - enlarge;  
      slabs_(i,1) = slabs_(i,1) + enlarge;
    }
   
    //Prints slabs to std::cout
    //PrintSlabs(); 
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | Print type of treenode to std::cout (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::PrintType()
{
  if (type_ == SELF_INNER)
    cout << endl << "SELF_INNER ";
  else if (type_ == SELF_LEAF)
    cout << endl << "SELF_LEAF ";
  else if (type_ == SELF_NO_ELEMENTS)
    cout << endl << "TreeNode contains no elements = SELF_NO_ELEMENTS ";
  else
    cout << endl << "SELF_UNDEFINED ";
  
  return;
}

/*----------------------------------------------------------------------*
 | Print slabs to std::cout (public)                          popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::PrintSlabs()
{
   cout << endl << Comm().MyPID() << "************************************************************";
   PrintType();
   cout << "slabs:";
   for (int i=0;i<slabs_.M();++i)
     cout << "\nslab: " << i << " min: " << slabs_.operator ()(i,0) << " max: " << slabs_.operator ()(i,1);
   cout << "\n**********************************************************\n";
   
   return;
}

/*----------------------------------------------------------------------*
 | Set slabs of current treenode with new slabs (public)      popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::SetSlabs(Epetra_SerialDenseMatrix& newslabs)
{
  for (int i=0;i<kdop_/2;++i)
  {
    slabs_(i,0) = newslabs(i,0);
    slabs_(i,1) = newslabs(i,1);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | Set children of current treenode (public)                  popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::SetChildren(RCP<BinaryTreeSelfNode> leftchild,
                                              RCP<BinaryTreeSelfNode> rightchild)
{
  leftchild_= leftchild;
  rightchild_ = rightchild;
  
  return;
}

/*----------------------------------------------------------------------*
 | Enlarge geometry of treenode (public)                      popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelfNode::EnlargeGeometry(double& enlarge)
{
  // scale slabs with scalar enlarge
  for (int i=0;i<kdop_/2;++i)
  {
    slabs_(i,0) = slabs_(i,0) - enlarge;  
    slabs_(i,1) = slabs_(i,1) + enlarge;
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | Constructor DualEdge (public)                              popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::DualEdge::DualEdge(RCP<BinaryTreeSelfNode> node1,
                            RCP<BinaryTreeSelfNode> node2,
                            const int& dim) :                                                     
node1_(node1),
node2_(node2),
dim_(dim)
{
  // directly move on to cost function
  CalculateCosts();
  
  return;
} 

/*----------------------------------------------------------------------*
 | Calculate cost function value (public)                     popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::DualEdge::CalculateCosts()
{ 
  // update slabs of both dual edge nodes
  double enlarge = 0.0;
  node1_->UpdateSlabsBottomUp(enlarge);
  node2_->UpdateSlabsBottomUp(enlarge);
  
  // build parent slab for dual edge
  Epetra_SerialDenseMatrix parentslabs;
  int n = node1_->Kdop();
  parentslabs.Reshape(n/2,2);
  
  for (int k=0;k<n/2;++k)
  {
    //for minimum
    if (node1_->Slabs()(k,0) <= node2_->Slabs()(k,0))
      parentslabs(k,0) = node1_->Slabs()(k,0);
    else 
      parentslabs(k,0) = node2_->Slabs()(k,0);
      
    // for maximum
    if (node1_->Slabs()(k,1) >= node2_->Slabs()(k,1))
      parentslabs(k,1) = node1_->Slabs()(k,1);
    else
      parentslabs(k,1) = node2_->Slabs()(k,1);    
  } 

  // compute maximal k-DOP length for dual edge
  double lmaxdop = 0.0;
  int slab = 0;
  for (int i=0;i<n/2;++i)
  {
     double lcurrent = abs(parentslabs(i,1)-parentslabs(i,0));
     if (lmaxdop < lcurrent)
     {
        lmaxdop = lcurrent;
        slab=i;
     }
  }
  
  // two-dimensional case
  if (dim_==2)
  {
    // compute total length of dual edge
    double lele = 0.0;
    
    for (int l=0;l<(int)(node1_->Elelist().size());++l)
    {
      int gid = (node1_->Elelist()).at(l);
      CElement* celement= static_cast<CElement*> (node1_->Discret().gElement(gid));
      lele = lele + celement->MaxEdgeSize();
    }
    
    for (int l=0;l<(int)(node2_->Elelist().size());++l)
    {
      int gid = node2_->Elelist()[l];
      CElement* celement= static_cast<CElement*> (node2_->Discret().gElement(gid));
      lele = lele + celement->MaxEdgeSize();
    }
    
    // cost function = nele * ( L / L_max )
    int nele = node1_->Elelist().size() + node2_->Elelist().size();
    costs_= nele * (lele/lmaxdop);
  }
  
  // three-dimensional case
  else 
  {
    // compute total area of dual edge
    double area = 0.0;
    
    for (int l=0;l<(int)(node1_->Elelist().size());++l)
    {
      int gid = (node1_->Elelist()).at(l);
      CElement* celement= static_cast<CElement*> (node1_->Discret().gElement(gid));
      area = area + celement->ComputeArea();
    }
    
    for (int l=0;l<(int)(node2_->Elelist().size());++l)
    {
      int gid = node2_->Elelist()[l];
      CElement* celement= static_cast<CElement*> (node2_->Discret().gElement(gid));
      area = area + celement->ComputeArea();
    }

    // compute maximal k-DOP area for dual edge
    double lmaxdop2 = 0.0;
    Epetra_SerialDenseMatrix dopnormals = node1_->Dopnormals();
    for (int j=0;j<n/2;++j)
    {
      double scalar = dopnormals(j,0)*dopnormals(slab,0)
                    + dopnormals(j,1)*dopnormals(slab,1)
                    + dopnormals(j,2)*dopnormals(slab,2);

      if (scalar==0)
      {
         double lcurrent2 = abs(parentslabs(j,1)-parentslabs(j,0));
         if (lmaxdop2 < lcurrent2) lmaxdop2 = lcurrent2;
      }
    }
    double doparea = lmaxdop * lmaxdop2;

    // cost function = nele * ( A / A_max )
    int nele = node1_->Elelist().size() + node2_->Elelist().size();
    costs_ = nele * nele * (area/doparea);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor BinaryTreeSelf (public)                              popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::BinaryTreeSelf::BinaryTreeSelf(DRT::Discretization& discret,
                                        RCP<Epetra_Map> elements,
                                        int dim, double eps) :
idiscret_(discret),
elements_(elements),
dim_(dim),
eps_(eps)
{
  // initialize sizes
  treenodes_.resize(1);
  leafsmap_.clear();

  // calculate min. element length and set enlargement accordingly
  SetEnlarge(true);

  //**********************************************************************
  // check for problem dimension
  //**********************************************************************
  // two-dimensional case
  if (dim_==2)
  {
    // set number of DOP sides
    kdop_=8;
    
    // set number of sample vectors
    nvectors_=16;
    
    // setup normals for DOP
    dopnormals_.Reshape(4,3); 
    dopnormals_(0,0)= 1; dopnormals_(0,1)= 0; dopnormals_(0,2)= 0;
    dopnormals_(1,0)= 0; dopnormals_(1,1)= 1; dopnormals_(1,2)= 0;
    dopnormals_(2,0)= 1; dopnormals_(2,1)= 1; dopnormals_(2,2)= 0;
    dopnormals_(3,0)=-1; dopnormals_(3,1)= 1; dopnormals_(3,2)= 0;
    
    // setup sample vectors
    samplevectors_.Reshape(16,3);
    samplevectors_(0,0)= 1.0; samplevectors_(0,1)= 0.0; samplevectors_(0,2)= 0.0;
    samplevectors_(8,0)=-1.0; samplevectors_(8,1)= 0.0; samplevectors_(8,2)= 0.0;
    
    for (int i=1;i<4;++i)
    {
      samplevectors_(i,0)= cos(PI*(double)i/8);
      samplevectors_(i,1)= sin(PI*(double)i/8);
      samplevectors_(i,2)= 0;
      
      samplevectors_(i+8,0)=   cos(PI*(double)i/8);
      samplevectors_(i+8,1)=-1*sin(PI*(double)i/8);
      samplevectors_(i+8,2)= 0;
    }
    
    samplevectors_(4,0)=  0.0; samplevectors_(4,1)=  1.0; samplevectors_(4,2)=  0.0;
    samplevectors_(12,0)= 0.0; samplevectors_(12,1)=-1.0; samplevectors_(12,2)= 0.0;
    
    for (int i=5;i<8;++i)
    {
       samplevectors_(i,0)= cos(PI*(double)i/8);
       samplevectors_(i,1)= sin(PI*(double)i/8);
       samplevectors_(i,2)= 0;
          
       samplevectors_(i+8,0)=   cos(PI*(double)i/8);
       samplevectors_(i+8,1)=-1*sin(PI*(double)i/8);
       samplevectors_(i+8,2)= 0;
    }
  }
  
  // three-dimensional case
  else if (dim_==3)
  {
    // set number of DOP sides
    kdop_=18;
    
    // set number of sample vectors
    nvectors_=50;
    
    // setup normals for DOP
    dopnormals_.Reshape(9,3); 
    dopnormals_(0,0)= 1; dopnormals_(0,1)= 0; dopnormals_(0,2)= 0;
    dopnormals_(1,0)= 0; dopnormals_(1,1)= 1; dopnormals_(1,2)= 0;
    dopnormals_(2,0)= 0; dopnormals_(2,1)= 0; dopnormals_(2,2)= 1;
    dopnormals_(3,0)= 1; dopnormals_(3,1)= 1; dopnormals_(3,2)= 0;
    dopnormals_(4,0)= 1; dopnormals_(4,1)= 0; dopnormals_(4,2)= 1;
    dopnormals_(5,0)= 0; dopnormals_(5,1)= 1; dopnormals_(5,2)= 1;
    dopnormals_(6,0)= 1; dopnormals_(6,1)= 0; dopnormals_(6,2)=-1;
    dopnormals_(7,0)= 1; dopnormals_(7,1)=-1; dopnormals_(7,2)= 0;
    dopnormals_(8,0)= 0; dopnormals_(8,1)= 1; dopnormals_(8,2)=-1;
     
    // setup sample vectors
    samplevectors_.Reshape(50,3);
    samplevectors_(0,0)= 0; samplevectors_(0,1)= 0.0; samplevectors_(0,2)= 1.0;
    samplevectors_(1,0)= 0; samplevectors_(1,1)= 0.0; samplevectors_(1,2)=-1.0;
    
    for (int i=0;i<8;++i)
    {
      for (int j=1;j<4;++j)
      {
        samplevectors_(1+6*i+j,0)= sin(PI*(double)j/8)*cos(PI*(double)i/8);
        samplevectors_(1+6*i+j,1)= sin(PI*(double)j/8)*sin(PI*(double)i/8);
        samplevectors_(1+6*i+j,2)= cos(PI*(double)j/8);
      }
      for (int j=5;j<8;++j)
      {
        samplevectors_(6*i+j,0)= sin(PI*(double)j/8)*cos(PI*(double)i/8);
        samplevectors_(6*i+j,1)= sin(PI*(double)j/8)*sin(PI*(double)i/8);
        samplevectors_(6*i+j,2)= cos(PI*(double)j/8);
      }
    }
  }
  
  else dserror("ERROR: Problem dimension must be 2D or 3D!");

  //**********************************************************************
  // initialize binary tree leaf nodes
  //**********************************************************************
  // create element lists
  vector<int> elelist;
  vector<int> localelelist;
  
  // build global element list
  for (int i=0;i<elements_->NumMyElements();++i)
  {
    int gid = elements_->GID(i);
    elelist.push_back(gid);
  }

  if (elelist.size()<=1) dserror("ERROR: Less than 2 elements for binary tree initialization!");
  

  // build local element list and create leaf nodes
  for (int i=0;i<(int)(elelist.size());++i)
  {
    localelelist.clear();
    localelelist.push_back(elelist[i]);
    RCP<BinaryTreeSelfNode> leaf = rcp(new BinaryTreeSelfNode(SELF_LEAF,idiscret_,null,localelelist,
       DopNormals(),SampleVectors(),Kdop(),Dim(),Nvectors(),-1,treenodes_));
    leafsmap_[elelist[i]]=leaf;
  }
  
  // double-check if there is at the least one leaf node in tree now
  if (leafsmap_.size()==0) dserror("ERROR: BinaryTreeSelf: No contact elements defined!");
  
  //**********************************************************************
  // find adjacent tree nodes and build dual graph
  //**********************************************************************
  // initialize dual graph
  map<RCP<DualEdge>, vector<RCP<DualEdge> > > dualgraph; 
  
  // loop over all self contact elements
  for (int i=0;i<(int)elelist.size();++i)
  {
    // global id of current element
    int gid = elelist[i];

    // initialize
    // ... vector of adjacent treenodes (elements) of current element
    // ... vector of global IDs of adjacent elements of current element
    // ... vector of adjacent dual edges containing current element
    // ... vector of global IDs of possible adjacent elements 
    vector<RCP<BinaryTreeSelfNode> >  adjtreenodes;
    vector<int>                       adjids;
    vector<RCP<DualEdge> >            adjdualedges;
    vector<int>                       possadjids; 

    // get curent elements and its nodes
    DRT::Element* element= idiscret_.gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n",gid);
    DRT::Node** nodes = element->Nodes();
    if (!nodes) dserror("ERROR: Null pointer!");
     
    // first treenode of one dual edge which includes current element
    // is the element itself saved as a treenode
    RCP<BinaryTreeSelfNode> node1=leafsmap_[gid];  
    
    // for 2D only: get the finite element nodes of the  element
    // and save them as endnodes of the treenode    
    if(dim_==2)
    {
      vector<int> nodeIds;
      nodeIds.push_back(element->NodeIds()[0]);
      nodeIds.push_back(element->NodeIds()[1]);
      node1->SetEndnodes(nodeIds);
    }
    
    // loop over all nodes of current element
    for(int j=0; j<element->NumNode();j++)
    {
      DRT::Node* node =nodes[j];
      if (!node) dserror("ERROR: Null pointer!");
        
      // adjacent elements of current node
      int numE =node->NumElement();
      DRT::Element** adjElements = node->Elements();
      if (!adjElements) dserror("ERROR: Null pointer!");
       
      // loop over all adjacent elements of current node
      for(int k=0;k<numE;++k)
      {
        int eleID=adjElements[k]->Id();
           
        // if eleID isn't the currently considered element, it could be a neighbour
        if (eleID!=gid)
        {
          // if there are not yet any possibly adjacent elements
          if (possadjids.size()==0)
          {
            // in 2D one common node implies adjacency
            if(dim_==2)
            {
              // get second node from leafsmap
              RCP<BinaryTreeSelfNode> node2 = leafsmap_[eleID];
              if (node2 == null) dserror("adjacent leaf tree node not found in leafsmap!!");
              
              // get the finite element nodes of the element equal to treenode 2
              // and save them as endnodes of the treenode
              vector<int> nodeIds;
              nodeIds.clear();
              nodeIds.push_back(adjElements[k]->NodeIds()[0]);
              nodeIds.push_back(adjElements[k]->NodeIds()[1]);
              node2->SetEndnodes( nodeIds );
              
              // add treenode in list of adjacent treenodes of current treenode
              adjtreenodes.push_back(node2);
              adjids.push_back(eleID);
              
              // create edge and add it to the list
              RCP<DualEdge> edge = rcp(new DualEdge(node1,node2,Dim()));
              adjdualedges.push_back(edge);
            }
            // in 3D adjacency is more complicated
            else
            {
              // get second node from leafsmap
              RCP<BinaryTreeSelfNode> node2 = leafsmap_[eleID];
              adjtreenodes.push_back(node2);
              possadjids.push_back(eleID);
            }
          }
          
          // if there are already possibly adjacent elements in possasdjids
          else
          {
            bool saved=false;
             
            // in 2D one common node implies adjacency
            if(dim_==2)
            {
              //get second node from leafsmap
              RCP<BinaryTreeSelfNode> node2 = leafsmap_[eleID];
              if (node2 == null) dserror("adjacent tree node not found in leafsmap!!");
              
              // get the finite element nodes of the element equal to treenode 2
              // and save them as endnodes of the treenode
              vector<int> nodeIds;
              nodeIds.push_back(adjElements[k]->NodeIds()[0]);
              nodeIds.push_back(adjElements[k]->NodeIds()[1]);
              node2->SetEndnodes(nodeIds);
              
              // add treenode in list of adjacent treenodes of current treenode
              adjtreenodes.push_back(node2);
              adjids.push_back(eleID);
              
              // create edge and add it to the list
              RCP<DualEdge> edge = rcp(new DualEdge(node1,node2,Dim()));
              adjdualedges.push_back(edge);
            }            
            // in 3D adjacency is more complicated
            // (adjacent elements have at least 2 common nodes)
            else
            {
              // get second node from leafsmap
              RCP<BinaryTreeSelfNode> node2 = leafsmap_[eleID];
              
              for (int l=0;l<(int)possadjids.size();++l)
              {
                // check if possibly adjacent element is already in the list
                // if true, there are 2 common nodes, which means it is a neighbour
                if (eleID==possadjids[l])
                {
                  saved=true;
                  if (node2 == null) dserror("adjacent tree node not found in leafsmap!!");
                  
                  // add treenode in list of adjacent treenodes of current treenode
                  // adjtreenodes.push_back(node2);
                  adjids.push_back(eleID);
                  
                  // create edge and add it to the list
                  RCP<DualEdge> edge = rcp(new DualEdge(node1,node2,Dim()));
                  adjdualedges.push_back(edge);
                  break;
                }
              }
              
              // possibly adjacent element is not yet in the list --> add
              if (!saved)
              {
                possadjids.push_back(eleID);
                adjtreenodes.push_back(node2);
              }
            } // else 3D
          } // else possadjids empty
        } // if eleID!=gid
      } // all adjacent elements
    } // all nodes
     
    // add the vector of adjacent tree nodes to the adjcencymatrix
    // we only need the matrix in 3D, because in 2D the adjacencytest works by
    // comparing endnodes only
    if (dim_==3) adjacencymatrix_[gid]=adjtreenodes;
    
    //get adjacent dual edges
    for(int k=0;k<(int)adjdualedges.size();++k)
    {
      for(int j=0;j<(int)adjdualedges.size();++j)
      {
        if (j!=k)
        {
          adjdualedges[k]->CommonNode(adjdualedges[j]);
          dualgraph[adjdualedges[k]].push_back(adjdualedges[j]);
        }
      }
    }
  } // all elements
   
  // plot adjacencymatrix
  /*
  map<int,vector<RCP<BinaryTreeSelfNode> > > ::iterator iter2 = adjacencymatrix_.begin();
  map<int,vector<RCP<BinaryTreeSelfNode> > > ::iterator iter2_end = adjacencymatrix_.end();
    
  cout << "\n" << leafsmap_.size() << " elements in leafsmap\n"; 
  cout << adjacencymatrix_.size() << " elements in adjazencymatrix\n";
  
  
  while (iter2 != iter2_end)
  {
    cout << "element " << (*iter2).first << ": ";
    //int adj_cnt= 0 ;
    
    vector<RCP<BinaryTreeSelfNode > >  adj_ = (*iter2).second;
    cout << adj_.size() << " elements: ";
    for (int i=0;i<(int)adj_.size();++i)
      cout << adj_[i]->Elelist()[0] << " ";
    cout << "\n";
    ++iter2;
  }
  */

  // plot dual graph
  /*
  cout << "\n" <<leafsmap_.size() << " elements in leafmap\n";
  cout << dualgraph.size() << " edges in dual graph\n";

  map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator iter3 = (*dualGraph)_.begin();
  map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator iter3_end = (*dualGraph)_.end(); 

  cout << (*dualGraph)_.max_size() << " maximal\n";
  int cnt = 0;
  
  while (iter3 != iter3_end)
  {
    cout << "\n Kante " << cnt << ": " << ((*iter3).first)->GetNode1()->Elelist()[0] << " "
         << ((*iter3).first)->GetNode2()->Elelist()[0] << "\n";
    cout << "Kosten: " << ((*iter3).first)->Costs() << "\n";

    vector<RCP<DualEdge> > edges = (*iter3).second;
    cout << edges.size() << " NachbarKanten: ";
    for (int i=0;i<(int)edges.size();++i)
    {
      cout << edges[i]->GetNode1()->Elelist()[0] << " ";
      cout << edges[i]->GetNode2()->Elelist()[0] << " ";
      cout << "Kosten: " << edges[i]->Costs() << " \n";
    }

    cout << "\n";
    ++iter3;
    ++cnt;
  }
  */

  // now initialze BinaryTreeSelf in a bottom-up way based on dualgraph
  InitializeTreeBottomUp(&dualgraph);
   
  return;
}

/*----------------------------------------------------------------------*
 | Find minimal length of contact elements (public)           popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::SetEnlarge(bool isinit)
{
  // minimal length of finite elements
  double lmin = 1.0e12;
  
  // calculate minimal length
  for (int i=0;i<elements_->NumMyElements();++i)
  {
    int gid = elements_->GID(i);
    DRT::Element* element = idiscret_.gElement(gid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n",gid);
    CONTACT::CElement* celement = (CElement*) element;
    double mincurrent=celement->MinEdgeSize(isinit);
    if (mincurrent < lmin) lmin = mincurrent;
  }
  
  if (lmin<=0.0) dserror("ERROR: Minimal element length < 0!");
  
  // set the class variable
  enlarge_ = eps_ * lmin;
  
  return;
}

/*----------------------------------------------------------------------*
 | Initialize tree bottum-up based on dual graph (public)     popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::InitializeTreeBottomUp(map <RCP<DualEdge>,vector <RCP<DualEdge> > > *dualGraph)
{
  // temporary vector collecting root nodes  
  vector<RCP<BinaryTreeSelfNode> > roots(0);
  
  // the idea is to empty the dual graph step by step
  while (!(*dualGraph).empty())
  {
    // get the edge with lowest costs (= the first edge in the dual graph as
    // the map is automatically sorted by the costs)  to contract it
    map <RCP<DualEdge>,vector <RCP<DualEdge> > > ::iterator iter = (*dualGraph).begin();
    map <RCP<DualEdge>,vector <RCP<DualEdge> > > ::iterator iter_end = (*dualGraph).end();
    RCP<DualEdge> contractedEdge = (*iter).first;
    RCP<BinaryTreeSelfNode> node1 = contractedEdge->GetNode1();
    RCP<BinaryTreeSelfNode> node2 = contractedEdge->GetNode2();
      
    // combine list of elements of both tree nodes to get new list
    vector<int> list = node1->Elelist();
    vector<int> list2 = node2->Elelist();
    for (int i=0;i<(int)(node2->Elelist().size());++i)
      list.push_back(list2[i]);
  
    // define new tree node
    RCP<BinaryTreeSelfNode> newNode = rcp(new BinaryTreeSelfNode(SELF_INNER,idiscret_,null,list,
                            DopNormals(),SampleVectors(),Kdop(),Dim(),Nvectors(),-1,treenodes_)); 
    newNode->SetChildren(node1,node2);
    node1->SetParent(newNode);
    node2->SetParent(newNode);

    // in 2D we simply save the endnodes as adjacency criterion
    if (dim_==2) newNode->UpdateEndnodes();
           
    // update dualGraph
    // this means we have to create new edges, which include the new treenode and
    // delete the edges, which are adjacent to the contracted edge and update
    // their neighbours
    // additionally we have to check if the tree is nearly complete
    
    // get the adjacent edges of the contracted edge
    vector<RCP<DualEdge> > adjEdges = (*iter).second;
    
    // check if the new treenode includes whole selfcontact-surface
    // in this case the treenode has saved itself as adjacent edge
    if (adjEdges[0] == contractedEdge)
    {
      // save the tree node as root and leave the loop
      roots.push_back(newNode);
      (*dualGraph).erase(contractedEdge);
      break;
    }
    
    // vector of all new edges
    vector<RCP<DualEdge> > newAdjEdges;
    
    // when contracting an edge, we need to update all adjacent edges
    for (int j=0;j<(int)adjEdges.size();++j)
    {
      // define new edge
      RCP<DualEdge> newEdge;
      
      if (adjEdges[j]->GetNode1() != node1 && adjEdges[j]->GetNode1() !=  node2)
        newEdge = rcp(new DualEdge(newNode,adjEdges[j]->GetNode1(),Dim()));
      
      else if (adjEdges[j]->GetNode2() != node1 && adjEdges[j]->GetNode2() !=  node2)
        newEdge = rcp(new DualEdge(newNode,adjEdges[j]->GetNode2(),Dim()));
      
      else dserror("ERROR: Tried to contract identical tree nodes!!");
   
      // get the neighbours of the new edges
      vector<RCP<DualEdge> > adjEdgesOfNeighbour = (*dualGraph)[adjEdges[j]];
      
      // if there are two adjacent surfaces left,
      // the new edge contains the whole self-contact surface.
      // in this case we save the edge as neighbour of itself
      if (adjEdges.size()==1 && adjEdgesOfNeighbour.size()==1)
        (*dualGraph)[newEdge].push_back(newEdge);
      
      // otherwise we need to update the neighbours of the adjacent edges:
      // first add all neighbours of the current adjacent edge -except the contracted edge-
      // as neighbour of newedge and delete the old adjacent edge
      else
      {
        // two-dimensional case
        if (dim_==2)
        {
          // in 2D every edge has 2 neighbours at the most
          newAdjEdges.push_back(newEdge);
          for (int k=0;k<(int)adjEdgesOfNeighbour.size();++k)
          {
            if (adjEdgesOfNeighbour[k] != contractedEdge)
            {
              // save the neighbours of the old edge as neighbours of the new one
              (*dualGraph)[newEdge].push_back(adjEdgesOfNeighbour[k]);
              
              // we have to update all neighbours of our old edge
              map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator edge_iter = (*dualGraph).find(adjEdgesOfNeighbour[k]);
              
              // find the old edge in the list of neighbours and replace it by the new edge
              for (int a=0;a<(int)edge_iter->second.size();++a)
              {
                if (edge_iter->second.at(a) == adjEdges[j])
                  edge_iter->second.at(a) = newEdge;
              }
            } 
          }
        }
        // three-dimensional case
        else
        {
          // in 3D every edge can have a arbirary number of neighbours which 
          // could also be adjacent to each other (3 edges form a "ring").
          // We have to check if the new edge has already been created. 
          // (in 2D this occurs only when the tree is almost finished, in 3D
          // this could happen any time)
 
          bool issaved = false;
          for (int n=0;n<(int)newAdjEdges.size();++n)
          {
            if (newAdjEdges[n]==newEdge)
            {
              issaved=true;
              break;
            }
          }
          
          // save the neighbours of the old edge as neighbours of the new one
          if (!issaved) newAdjEdges.push_back(newEdge);
          
          // we have to update all neighbours of the old/new edge
          // we just update those edges wich aren't adjacent edges of the contracted
          // edge, as we have to update those seperately
          
          // get the common (tree)node of the contracted edge and the old edge 
          RCP<BinaryTreeSelfNode> commonnode = contractedEdge->CommonNode(adjEdges[j]);
          
          // loop over all neighbours of old edge
          for(int k=0;k<(int)adjEdgesOfNeighbour.size();++k)
          {
            // check if it the current edge is not the contracted edge
            // or a neigbour of the contracted edge (we do not need to update
            // these edges now)
            if (adjEdgesOfNeighbour[k] != contractedEdge &&
                adjEdges[j]->CommonNode(adjEdgesOfNeighbour[k]) != commonnode )
            {
              // if the current edge (=neighbour of the new and old edge) is not a
              // neighbour of the contracted edge, they do not have a common node
              RCP<BinaryTreeSelfNode> commonnode2 = contractedEdge->CommonNode(adjEdgesOfNeighbour[k]);
                
              // now we want to save the current edge as neighbour of the new edge
              if (commonnode2==null)
              {
                // first find the new edge in the dual graph
                map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator edge_iter1 = (*dualGraph).find(newEdge);
                map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator end =(*dualGraph).end();
                
                // if the edge has been found, check if the current edge has
                // already been saved as neighbour
                if (edge_iter1 != end)
                {
                  bool edgesaved = false;
                  for(int z=0;z<(int)edge_iter1->second.size();++z)
                  {
                    if (edge_iter1->second.at(z) == adjEdgesOfNeighbour[k])
                      edgesaved=true;
                  }
                  // if not yet saved, save it
                  if (!edgesaved) edge_iter1->second.push_back(adjEdgesOfNeighbour[k]);
                }
                
                // if the new edge itself hasn't been saved yet,
                // the neighbour hasn't been saved either
                else  (*dualGraph)[newEdge].push_back(adjEdgesOfNeighbour[k]);
              
                // find the old edge in the list of neighbours and replace it
                // by the new edge
                
                // find the current edge in the dual graph
                map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator edge_iter2 = (*dualGraph).find(adjEdgesOfNeighbour[k]);
         
                bool egdeerased = false;
                bool newedgesaved = false;
                
                // loop over all neighbours of the current edge
                vector<RCP<DualEdge> > ::iterator adjIter = edge_iter2->second.begin();
                while (adjIter != edge_iter2->second.end())
                {
                  if (*adjIter == adjEdges[j])
                  {
                    // erase the old edge
                    adjIter = edge_iter2->second.erase(adjIter);
                    egdeerased = true;
                  }
                  else
                  {
                    if (*adjIter == newEdge) newedgesaved=true;
                    ++adjIter;
                  }
                }
                
                // as we could update the same edge several times (only in 3D), 
                // we only save the new edge as neighbour if not already done so
                if (egdeerased && !newedgesaved) edge_iter2->second.push_back(newEdge);
              }
            } 
          } // loop over all adjacent edges of neighbours
          
          // if three dual edges are forming a "ring" there will be just one edge left, so
          // we save this egde as neighbour of itself (in 2D we do not need to consider
          // this case seperately)
          if ((*dualGraph).size()==3 && newAdjEdges.size()==1)
          {
            map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator edge_iter3 = (*dualGraph).find(newAdjEdges[0]);
            map<RCP<DualEdge>,vector<RCP<DualEdge> > > ::iterator end =(*dualGraph).end();
            
            if (edge_iter3 == end) (*dualGraph)[newEdge].push_back(newEdge);
          }
        } // 3D
      } // else-block (not 2 adjacent treenodes left)
    } // loop over all adjacent edges
    
    // delete all adjacent edges of contracted edge
    for (int j=0;j<(int)adjEdges.size();++j)
      (*dualGraph).erase(adjEdges[j]);
    
    // now all new adjacent edges have been created.
    // save all adjacent edges as neigbour, respectively
    for (int l=0;l<(int)newAdjEdges.size();++l)
    {
      for (int m=0;m<(int)newAdjEdges.size();++m)
        if (l!=m && newAdjEdges[l]!=newAdjEdges[m])
          (*dualGraph)[ newAdjEdges[l] ].push_back(newAdjEdges[m]); 
    }

    // delete the contraceted edge    
    (*dualGraph).erase(contractedEdge);
   
  } // while(!(*dualGraph).empty())
  
  // complete the tree starting from its root (top-down)
  if ((int)roots.size() == 0) dserror("ERROR: No root treenode found!");
  if ((int)roots.size() > 1)  dserror("ERROR: Disconnected self contact surface!"); 
  root_=roots[0];
  root_->CompleteTree(0,enlarge_);
    
  // in 3D we have to calculate adjacent treenodes
  if (dim_== 3)
  {
    CalculateAdjacentLeaves();
    CalculateAdjacentTnodes();
  }
 
  /*
  // print root
  vector<int> rootElelist = root_->Elelist();
  cout << "\nWurzel: ";
  for (int d=0;d<(int)rootElelist.size();++d)
    cout << rootElelist[d] << " ";
  cout << "\n";
  
  // print tree
  cout << "SelfContactTree: \n";
  for (int i=0;i<(int)(treenodes_.size());++i)
  {
    cout << "\n Tree at layer: " << i << " Elements: ";
    for (int k=0;k<(int)(treenodes_[i].size());++k)
    {
      RCP<BinaryTreeSelfNode> currentnode = treenodes_[i][k];
      cout << " (";
      for (int l=0;l<(int)(currentnode->Elelist().size());++l)
      {
        cout << currentnode->Elelist().at(l) << " ";
        if (currentnode->Type()== SELF_LEAF) cout << "(Leaf) ";
      }
      cout << ") ";
    }
  }
  */
  
  return;
}

/*----------------------------------------------------------------------*
 | set adjacent tree nodes of leaf-nodes in the lowest layer (3D-SC) (public)    popp 07/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::CalculateAdjacentLeaves()
{
  //get the adjacent treenodes of each treenode in the lowest layer
  //and save the adjacent leafs which are in the same layer
  int maxlayer = treenodes_.size()-1;
  map<int, RCP<BinaryTreeSelfNode> > ::iterator leafiter = leafsmap_.begin(),
                            leafiter_end = leafsmap_.end();
  while(leafiter != leafiter_end)
  {
    if(leafiter->second->Layer() == maxlayer)
    {
      vector<RCP<BinaryTreeSelfNode> > adjtnodessamelayer_;
      vector<RCP<BinaryTreeSelfNode> > adjtnodes_=adjacencymatrix_[leafiter->first];
      for(int i=0; i<(int)adjtnodes_.size();i++)
      {
        if(adjtnodes_[i]->Layer() == maxlayer)
          adjtnodessamelayer_.push_back(adjtnodes_[i]);
      }
      leafiter->second->SetAdjacentTnodes(adjtnodessamelayer_);
    }
    leafiter++;
  }
    
  return;
}

/*----------------------------------------------------------------------*
 | set adjacent tree nodes of the tree (3D-SC) (public)    popp 07/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::CalculateAdjacentTnodes()
{
  //calculate adjacent treenodes in the same layer of the each treenode above the leaf-layer
  //they are calculated in a bottom up way, so the adjacent treenodes of the
  //lowest layer must have been calculated 
  
  int maxlayer = treenodes_.size();
  for(int i = maxlayer-2; i >= 0;i--)
  {
    for(int j=0; j< (int)treenodes_.at(i).size();j++)
    {
      //vector of adjacent treenodes
      vector<RCP<BinaryTreeSelfNode> > adjtnodes_;
      
      if(treenodes_[i][j]->Type() != SELF_LEAF)
      {
        //get the adjacent treenodes of the children
        vector<RCP<BinaryTreeSelfNode> > adjofleftchild_ = treenodes_[i][j]
                                                       ->Leftchild()->AdjacentTreenodes(),
                                    adjofrightchild_ = treenodes_[i][j]
                                                       ->Rightchild()->AdjacentTreenodes();
        
        //check the adjacent treenodes of the left child
        for(int k=0;k< (int)adjofleftchild_.size();k++)
        {
          //check if the parent of the adjacent node of the child has already been saved
          if(adjofleftchild_[k] != treenodes_[i][j]->Rightchild())
          {
            bool issaved=false;
            for(int l=0;l<(int)adjtnodes_.size();l++)
            {
              if(adjtnodes_[l] == null)
                dserror("null pointer");
              if(adjofleftchild_[k]->Parent() == adjtnodes_[l])
              {
                issaved=true;
                break;
              }
            }
            
            if(!issaved)
              adjtnodes_.push_back(adjofleftchild_[k]->Parent());
          }
        }
        
        //check the adjacent treenodes of the right child
        for(int k=0;k< (int)adjofrightchild_.size();k++)
        {
          //check if the parent of the adjacent node of the child has already been saved
          if(adjofrightchild_[k] != treenodes_[i][j]->Leftchild())
          {
            bool issaved=false;
            for(int m=0;m<(int)adjtnodes_.size();m++)
            {
              if(adjofrightchild_[k]->Parent() == adjtnodes_[m])
              {
                issaved=true;
                break;
              }
            }
            if(!issaved)
              adjtnodes_.push_back(adjofrightchild_[k]->Parent());
          }
        }
        
        treenodes_[i][j]->SetAdjacentTnodes(adjtnodes_);  
        
      }
      
      else //tree node is a leaf, which is not in the lowest layer
      {
        //get the adjacent leaf nodes from the adjacencymatrix
        int gid = treenodes_[i][j]->Elelist()[0];
        if(adjacencymatrix_.find(gid)==adjacencymatrix_.end())
          dserror("element not in adjazencymatrix!!");
        vector<RCP<BinaryTreeSelfNode> > adjleafs_=adjacencymatrix_[gid];
        
        for(int n=0; n<(int)adjleafs_.size();n++)
        {
          //for each adjacent leaf find the parent, which is on the same layer
          //as the current treenode
            RCP<BinaryTreeSelfNode>  adjtnode = adjleafs_[n];
            int diff= adjleafs_[n]->Layer()-i;
            if(diff>=0)
            {
              while(diff > 0)
              {
                adjtnode = adjtnode->Parent();
                diff--;
              }
              
              //check if the treenode has already been saved as adjacent node
              bool issaved=false;
              for(int p=0;p<(int)adjtnodes_.size();p++)
              {
                
                if(adjtnode==null)
                  dserror("null vector!!");
                
                if(adjtnode == adjtnodes_[p])
                {
                  issaved=true;
                  break;
                }
              }
              if(!issaved)
                adjtnodes_.push_back(adjtnode);
            }
        }
        //save the found treenodes as adjacent treenodes
        treenodes_[i][j]->SetAdjacentTnodes(adjtnodes_);  
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | Search for self contact (public)                                popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::SearchSelfContact(RCP<BinaryTreeSelfNode> treenode)
{
  
  if (treenode->QualifiedVectors().size()==0)
    dserror("no test vectors defined!");
  
  //if there is a qualified sample vector, there is no self contact
  for(int i=0; i < (int)treenode->QualifiedVectors().size();i++)
    if(treenode->QualifiedVectors()[i] == true)
    {
      return;
    }
  
  if (treenode->Type() != SELF_LEAF)
  {
    treenode->Leftchild()->CalculateSlabsDop(false);
    treenode->Leftchild()->EnlargeGeometry(enlarge_);
    treenode->Rightchild()->CalculateSlabsDop(false);
    treenode->Rightchild()->EnlargeGeometry(enlarge_);
    SearchSelfContact(treenode->Leftchild());
    SearchSelfContact(treenode->Rightchild());
    EvaluateContactAndAdjacency(treenode->Leftchild(),
                                treenode->Rightchild(),true);
  }

  return;
}

/*----------------------------------------------------------------------*
 | find contact and test adjacency (public)                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::EvaluateContactAndAdjacency(RCP<BinaryTreeSelfNode> treenode1,
                                                      RCP<BinaryTreeSelfNode> treenode2,
                                                       bool isadjacent)
{
  
  // check if treenodes intercept
  // (they only intercept if ALL slabs intersect!)
  int nintercepts = 0;
  
  for (int i=0;i<kdop_/2;++i)
  {
    if (treenode1->Slabs()(i,0) <= treenode2->Slabs()(i,0))
    {
      if (treenode1->Slabs()(i,1) >= treenode2->Slabs()(i,0))
        nintercepts++;
      else if (treenode1->Slabs()(i,1) >= treenode2->Slabs()(i,1))
        nintercepts++;
    }
    else if (treenode1->Slabs()(i,0) >= treenode2->Slabs()(i,0))
    {
      if (treenode2->Slabs()(i,1) >= treenode1->Slabs()(i,1))
        nintercepts++;
      else if (treenode2->Slabs()(i,1) >= treenode1->Slabs()(i,0))
        nintercepts++;
    } 
  }
  
  if (nintercepts==kdop_/2)
  {
    //teenodes intersect
      if (isadjacent)
      {
        if(dim_==2)
          isadjacent = TestAdjacent2D(treenode1,treenode2);
        else
        {
          isadjacent = TestAdjacent3D(treenode1,treenode2);
        }
      
        if(isadjacent)
        {
            vector<bool> qualifiedvectors1 = treenode1->QualifiedVectors();
            vector<bool> qualifiedvectors2 = treenode2->QualifiedVectors();
            
            if((int)qualifiedvectors1.size() == 0 or (int)qualifiedvectors2.size()==0)
                 dserror("no test vectors defined!");
            
            if((int)qualifiedvectors1.size() != (int)qualifiedvectors2.size())
              dserror("not the same number of test vectors!");
           
            for(int i=0; i<(int)qualifiedvectors1.size();i++)
            {
              if(qualifiedvectors1[i] and qualifiedvectors2.at(i))
              {
                return;
              }
            }
        }
      }
      
      if((int)treenode1->Elelist().size() > (int)treenode2->Elelist().size())
      {
        if (treenode1->Type() != SELF_LEAF)
        {
          treenode1->Leftchild()->CalculateSlabsDop(false);
          treenode1->Leftchild()->EnlargeGeometry(enlarge_);
          treenode1->Rightchild()->CalculateSlabsDop(false);
          treenode1->Rightchild()->EnlargeGeometry(enlarge_);
          EvaluateContactAndAdjacency(treenode1->Leftchild(),treenode2, isadjacent);
          EvaluateContactAndAdjacency(treenode1->Rightchild(),treenode2, isadjacent);
        }
      }
      
      else 
      {
         if (treenode2->Type() != SELF_LEAF)
         {
           treenode2->Leftchild()->CalculateSlabsDop(false);
           treenode2->Leftchild()->EnlargeGeometry(enlarge_);
           treenode2->Rightchild()->CalculateSlabsDop(false);
           treenode2->Rightchild()->EnlargeGeometry(enlarge_);
           EvaluateContactAndAdjacency(treenode2->Leftchild(),treenode1, isadjacent);
           EvaluateContactAndAdjacency(treenode2->Rightchild(),treenode1, isadjacent);
         }
         else //both tree nodes are leaves
         {
           if(dim_==2)
             isadjacent = TestAdjacent2D(treenode1,treenode2);
           if(dim_==3)
             isadjacent = TestAdjacent3D(treenode1,treenode2);
           if(!isadjacent)
           {
           contactpairs_[treenode1->Elelist()[0]].push_back(treenode2->Elelist()[0]);
           contactpairs_[treenode2->Elelist()[0]].push_back(treenode1->Elelist()[0]);
           }
         }
      }
        
  }
  else    //dops do not intercept;
    return;
}

/*----------------------------------------------------------------------*
 | find contact and test adjacency (public)                    popp 06/09|
 *----------------------------------------------------------------------*/
bool CONTACT::BinaryTreeSelf::TestAdjacent2D(RCP<BinaryTreeSelfNode> treenode1,
                                                      RCP<BinaryTreeSelfNode> treenode2)
{
  if(dim_ != 2)
    dserror("TestAdjacent2D: problem must be 2D!!\n");
  
  vector<int> endnodes1 = treenode1->Endnodes();
  vector<int> endnodes2 = treenode2->Endnodes();
  
  if(endnodes1.size() != 2 or endnodes2.size() != 2)
    dserror("treenode has not 2 endnodes!!\n");
  
  for(int i=0; i<(int)endnodes1.size();i++)
  {
    if(endnodes1[i]== -1)
    {
      //treenode is a closed surface -> has no endnodes;
      return false;
    }
    
    for(int j=0; j<(int)endnodes2.size();j++)
    {
      if(endnodes2[j] == -1) 
      {
        //treenode is a closed surface -> has no endnodes;
        return false;
      }
      
      else if(endnodes1[i] == endnodes2[j])
      {
        //treenodes are adjacent
        return true;
      }

    }
    
  }
  //treenodes are not adjacent
  return false;
}

/*----------------------------------------------------------------------*
 | find contact and test adjacency (public)                    popp 06/09|
 *----------------------------------------------------------------------*/
bool CONTACT::BinaryTreeSelf::TestAdjacent3D(RCP<BinaryTreeSelfNode> treenode1,
                                                      RCP<BinaryTreeSelfNode> treenode2)
{
  if(dim_ != 3)
    dserror("TestAdjacent3D: problem must be 3D!!\n");
  
  //if the treenodes are in the same layer check the vector of adjacent treenodes
  if(treenode1->Layer() == treenode2->Layer())
  {
   
    vector<RCP<BinaryTreeSelfNode> > adjtnodes = treenode1->AdjacentTreenodes();
    for(int i=0; i<(int)adjtnodes.size();i++)
    {
      if(adjtnodes[i] == treenode2)
      {
        //   treenodes are adjacent
        return true;
      }
    }
  }
  
  else 
  {
    // check if bounding volumes overlap
    // (they only intercept if ALL slabs intercept!)
    int nintercepts = 0;
      
    for (int i=0;i<kdop_/2;++i)
    {
      if (treenode1->Slabs()(i,0) <= treenode2->Slabs()(i,0))
      {
        if (treenode1->Slabs()(i,1) >= treenode2->Slabs()(i,0))
          nintercepts++;
        else if (treenode1->Slabs()(i,1) >= treenode2->Slabs()(i,1))
          nintercepts++;
      }
      else if (treenode1->Slabs()(i,0) >= treenode2->Slabs()(i,0))
      {
        if (treenode2->Slabs()(i,1) >= treenode1->Slabs()(i,1))
          nintercepts++;
        else if (treenode2->Slabs()(i,1) >= treenode1->Slabs()(i,0))
          nintercepts++;
      } 
    }
    
    //if the bounding voumes overlap
    if (nintercepts==kdop_/2)
    {
      if (treenode1->Type()==SELF_LEAF and treenode2->Type()==SELF_LEAF)
      {
     // two leaves
        vector<RCP<BinaryTreeSelfNode> > adjleafs = adjacencymatrix_[treenode1->Elelist()[0]];
        for(int i=0; i<(int)adjleafs.size();i++)
        {
          if(treenode2 == adjleafs[i])
          {
           // leaves are adjacent
            return true;
          }
        }
      }
    
      //one leaf and one inner treenode
      else if (treenode1->Type()==SELF_LEAF and treenode2->Type()!=SELF_LEAF)
        return ( TestAdjacent3D(treenode1,treenode2->Leftchild())
                  or TestAdjacent3D(treenode1,treenode2->Rightchild()));
      
      else if(treenode1->Type()!=SELF_LEAF and treenode2->Type()==SELF_LEAF)
        return ( TestAdjacent3D(treenode1->Leftchild(),treenode2)
                or TestAdjacent3D(treenode1->Rightchild(),treenode2));
      
      else if((treenode1->Layer()) > treenode2->Layer())
          return ( TestAdjacent3D(treenode1,treenode2->Leftchild())
                    or TestAdjacent3D(treenode1,treenode2->Rightchild()));
      
      else if((treenode1->Layer()) < treenode2->Layer())
          return ( TestAdjacent3D(treenode1->Leftchild(),treenode2)
                    or TestAdjacent3D(treenode1->Rightchild(),treenode2));
    }
  }
  //treenodes do not overlap;
  return false;
}
                               
/*----------------------------------------------------------------------*
 | do Master/Self facet sorting (self contact) (public)        popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::MasterSlaveSorting(int eleID,bool isslave)
{
  
  if(contactpairs_.find(eleID) != contactpairs_.end() && !contactpairs_.empty())
  {
    //set the element as isslave
    DRT::Element* element= idiscret_.gElement(eleID);
    CONTACT::CElement* celement = static_cast<CONTACT::CElement*>(element);
    celement->SetSlave()=isslave;
    
    //if the element is a slave, set its node to slave (otherwise the nodes
    //are master-nodes already and nodes between a master and a slave element
    //should be slave nodes)
    if(celement->IsSlave())
      for(int i=0; i<(int)element->NumNode();i++)
      {
        DRT::Node* node = element->Nodes()[i];
        CONTACT::CNode* cnode = static_cast<CONTACT::CNode*>(node);
        cnode->SetSlave()=isslave;
      }
    
    //get the ID of elements in contact with current one
    vector<int> contacteleID = contactpairs_[eleID];
    //erase the current element of the list
    contactpairs_.erase(eleID);
    
    
    for(int i =0; i<(int)contacteleID.size();i++)
    {
      
     if(celement->IsSlave())
        celement->AddSearchElements(contacteleID[i]); 
     if(contactpairs_.find(contacteleID[i])!=contactpairs_.end())
     {
        MasterSlaveSorting(contacteleID[i], !isslave);
      }
    }
  }

  return;
}
                                                      
/*----------------------------------------------------------------------*
 | Search for contact between TreeNodes and automatically update 
 | treeEvaluate Binary search tree for combined search and update (public) popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::SearchContactCombined()
{
  contactpairs_.clear();

  if (root_==null) dserror("ERROR: No root node for search!");
  
  root_->CalculateSlabsDop(false);
  root_->EnlargeGeometry(enlarge_);
  
  UpdateNormals();

  SearchSelfContact(root_);
  
  //Slave/master facet sorting
 
  //set all contact elements and nodes on the surface to master
  map<int, RCP<BinaryTreeSelfNode> > ::iterator leafiter = leafsmap_.begin(),
                               leafiter_end = leafsmap_.end();
  while(leafiter != leafiter_end)
  {
    int gid = leafiter->first;
    DRT::Element* element= idiscret_.gElement(gid);
    CONTACT::CElement* celement = static_cast<CONTACT::CElement*>(element);
    
    if(celement->IsSlave() == true)
    {
      celement->SetSlave()=false;
      for(int i=0; i<(int)element->NumNode();i++)
      {
        DRT::Node* node = element->Nodes()[i];
        CONTACT::CNode* cnode = static_cast<CONTACT::CNode*>(node);
        cnode->SetSlave()=false;
      }
    }

    leafiter++;
  }
  
  while(!contactpairs_.empty())
  { 
    DRT::Element* element= idiscret_.gElement(contactpairs_.begin()->first);
    CONTACT::CElement* celement = static_cast<CONTACT::CElement*>(element);
    MasterSlaveSorting(contactpairs_.begin()->first,celement->IsSlave());
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update normals and qualified sample vectors of the whole tree(public) popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::UpdateNormals()
{
  //first update normals and sample vectors of all leaf-nodes
  
  map<int, RCP<BinaryTreeSelfNode> > ::iterator iter = leafsmap_.begin(),
                               iter_end = leafsmap_.end();
  while(iter != iter_end)
  {
    iter->second->CalculateQualifiedVectors();
    iter++;
  }
  
  //update the rest of the tree in a bottom up way
  for(int i=(int)treenodes_.size()-1; i>=0; i--)
  {
    for(int j=0; j<(int)treenodes_[i].size();j++)
      if(treenodes_[i][j]->Type() != SELF_LEAF)
        treenodes_[i][j]->UpdateQualifiedVectorsBottomUp();
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | Update normals and qualified sample vectors of one treenode(public) popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::UpdateQualifiedVectors(RCP<BinaryTreeSelfNode> treenode)
{
  vector<bool> qualifiedvectors;

  qualifiedvectors.resize(nvectors_);
  
  map<int, RCP<BinaryTreeSelfNode> >   leafsmap =   leafsmap_; 
  for(int j=0; j<(int)treenode->Elelist().size();j++)
  {
    leafsmap[treenode->Elelist()[j]]->CalculateQualifiedVectors();
  }
  
  for(int i =0; i<(int)qualifiedvectors.size(); i++)
  {
    bool isvalid=true;
    for(int j=0; j<(int)treenode->Elelist().size();j++)
    {
      RCP<BinaryTreeSelfNode> leaf = leafsmap[treenode->Elelist()[j]];
      if( !(leaf->QualifiedVectors()[i]) )
      {
        isvalid = false;
        break;
      }
    }
    qualifiedvectors.at(i) = isvalid;
  }
  treenode->SetQualifiedVectors(qualifiedvectors);
  return;
}

#endif //#ifdef CCADISCRET
