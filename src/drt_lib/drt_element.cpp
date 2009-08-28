/*!----------------------------------------------------------------------
\file drt_element.cpp
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_element.H"
#include "drt_node.H"
#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_utils.H"
#include "drt_linedefinition.H"

#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::StringToDistype(std::string name)
{
  static std::map<std::string,DRT::Element::DiscretizationType> gid2distype;
  if (gid2distype.size()==0)
  {
    gid2distype["HEX8"]     = DRT::Element::hex8;
    gid2distype["HEX20"]    = DRT::Element::hex20;
    gid2distype["HEX27"]    = DRT::Element::hex27;
    gid2distype["TET4"]     = DRT::Element::tet4;
    gid2distype["TET10"]    = DRT::Element::tet10;
    gid2distype["WEDGE6"]   = DRT::Element::wedge6;
    gid2distype["WEDGE15"]  = DRT::Element::wedge15;
    gid2distype["PYRAMID5"] = DRT::Element::pyramid5;
    gid2distype["QUAD4"]    = DRT::Element::quad4;
    gid2distype["QUAD8"]    = DRT::Element::quad8;
    gid2distype["QUAD9"]    = DRT::Element::quad9;
    gid2distype["TRI3"]     = DRT::Element::tri3;
    gid2distype["TRI6"]     = DRT::Element::tri6;
    gid2distype["NURBS2"]   = DRT::Element::nurbs2;
    gid2distype["NURBS3"]   = DRT::Element::nurbs3;
    gid2distype["NURBS4"]   = DRT::Element::nurbs4;
    gid2distype["NURBS8"]   = DRT::Element::nurbs8;
    gid2distype["NURBS9"]   = DRT::Element::nurbs9;
    gid2distype["NURBS27"]  = DRT::Element::nurbs27;
    gid2distype["LINE2"]    = DRT::Element::line2;
    gid2distype["LINE3"]    = DRT::Element::line3;
    gid2distype["POINT1"]   = DRT::Element::point1;
    gid2distype["DIS_NONE"] = DRT::Element::dis_none;
    gid2distype["MAX_DISTYPE"] = DRT::Element::max_distype;
  }

  std::map<std::string,DRT::Element::DiscretizationType>::iterator i;
  i = gid2distype.find(name);
  if (i!=gid2distype.end())
    return i->second;
  dserror("unsupported distype '%s'",name.c_str());
  return DRT::Element::dis_none;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(int id, ElementType etype, int owner) :
ParObject(),
id_(id),
lid_(-1),
owner_(owner),
etype_(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(const DRT::Element& old) :
ParObject(old),
id_(old.id_),
lid_(old.lid_),
owner_(old.owner_),
etype_(old.etype_),
nodeid_(old.nodeid_),
node_(old.node_)
{
  // we do NOT want a deep copy of the condition_ as the condition
  // is only a reference in the elements anyway
  map<string,RefCountPtr<Condition> >::const_iterator fool;
  for (fool=old.condition_.begin(); fool!=old.condition_.end(); ++fool)
    SetCondition(fool->first,fool->second);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::~Element()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Element& element)
{
  element.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::Print(ostream& os) const
{
  os << setw(12) << Id() << " Owner " << setw(5) << Owner() << " ";
  const int nnode = NumNode();
  const int* nodeids = NodeIds();
  if (nnode > 0)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << setw(10) << nodeids[i] << " ";
  }

#if 0
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << endl << numcond << " Conditions:\n";
    map<string,RefCountPtr<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << curr->first << " ";
      os << *(curr->second) << endl;
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  read element input dummy (public)                        mwgee 11/06|
 |  this is a base class dummy for elements that do not need            |
 |  a reading method
 *----------------------------------------------------------------------*/
bool DRT::Element::ReadElement()
{
  cout << "DRT::Element::ReadElement:\n"
       << "Base class dummy routine Element::ReadElement() called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::Element::ReadElement(const std::string& eletype,
                               const std::string& distype,
                               DRT::INPUT::LineDefinition* linedef)
{
  dserror("subclass implementations missing");
  return false;
}


/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::SetNodeIds(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Element::SetNodeIds(const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  linedef->ExtractIntVector(distype,nodeid_);
  for (unsigned i=0; i<nodeid_.size(); ++i)
    nodeid_[i] -= 1;
  node_.resize(0);
}


/*----------------------------------------------------------------------*
 |  create material class (public)                                 05/07|
 *----------------------------------------------------------------------*/
void DRT::Element::SetMaterial(int matnum)
{
  mat_ = MAT::Material::Factory(matnum);
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Pack(vector<char>& data) const
{
  data.resize(0);
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add id
  AddtoPack(data,id_);
  // add owner
  AddtoPack(data,owner_);
  // add type of element
  AddtoPack(data,etype_);
  // add vector nodeid_
  AddtoPack(data,nodeid_);
  // add material
  vector<char> tmp;
  if (mat_!=Teuchos::null)
  {
    mat_->Pack(tmp);
  }
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // id_
  ExtractfromPack(position,data,id_);
  // owner_
  ExtractfromPack(position,data,owner_);
  // etype_
  ExtractfromPack(position,data,etype_);
  // nodeid_
  ExtractfromPack(position,data,nodeid_);
  // mat_
  vector<char> tmp;
  ExtractfromPack(position,data,tmp);
  if (tmp.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(tmp);
    MAT::Material* mat = dynamic_cast<MAT::Material*>(o);
    if (mat==NULL)
      dserror("failed to unpack material");
    mat_ = rcp(mat);
  }
  else
  {
    mat_ = Teuchos::null;
  }

  // node_ is NOT communicated
  node_.resize(0);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(map<int,RefCountPtr<DRT::Node> >& nodes)
{
  int        nnode   = NumNode();
  const int* nodeids = NodeIds();
  node_.resize(nnode);
  for (int i=0; i<nnode; ++i)
  {
    map<int,RefCountPtr<DRT::Node> >::const_iterator curr = nodes.find(nodeids[i]);
    // this node is not on this proc
    if (curr==nodes.end()) dserror("Element %d cannot find node %d",Id(),nodeids[i]);
    else
      node_[i] = curr->second.get();
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(DRT::Node** nodes)
{
  node_.resize(NumNode());
  for (int i=0; i<NumNode(); ++i) node_[i] = nodes[i];
  return true;
}


/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::GetCondition(const string& name,vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  multimap<string,RefCountPtr<Condition> >::const_iterator startit =
                                         condition_.lower_bound(name);
  multimap<string,RefCountPtr<Condition> >::const_iterator endit =
                                         condition_.upper_bound(name);
  int count=0;
  multimap<string,RefCountPtr<Condition> >::const_iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Element::GetCondition(const string& name) const
{
  multimap<string,RefCountPtr<Condition> >::const_iterator curr =
                                         condition_.find(name);
  if (curr==condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const Discretization& dis,
                                  vector<int>& lm, vector<int>& lmdirich,
                                  vector<int>& lmowner) const
{
  const int numnode = NumNode();
  const DRT::Node*const* nodes = Nodes();

  lm.clear();
  lmdirich.clear();
  lmowner.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i=0; i<numnode; ++i)
    {
      DRT::Condition* dirich = nodes[i]->GetCondition("Dirichlet");
      const vector<int>* flag = NULL;
      if (dirich)
      {
        if (dirich->Type()!=DRT::Condition::PointDirichlet &&
            dirich->Type()!=DRT::Condition::LineDirichlet &&
            dirich->Type()!=DRT::Condition::SurfaceDirichlet &&
            dirich->Type()!=DRT::Condition::VolumeDirichlet)
          dserror("condition with name dirichlet is not of type Dirichlet");
        flag = dirich->Get<vector<int> >("onoff");
      }
      const int owner = nodes[i]->Owner();
      vector<int> dof = dis.Dof(nodes[i]);
      const int size = NumDofPerNode(*(nodes[i]));
      for (int j=0; j<size; ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
        lmowner.push_back(owner);
        lm.push_back(dof[j]);
      }
    }
  }

  // fill the vector with element dofs
  const vector<int>* flag = NULL;
  DRT::Condition* dirich = GetCondition("Dirichlet");
  if (dirich)
  {
    if (dirich->Type()!=DRT::Condition::PointDirichlet &&
        dirich->Type()!=DRT::Condition::LineDirichlet &&
        dirich->Type()!=DRT::Condition::SurfaceDirichlet &&
        dirich->Type()!=DRT::Condition::VolumeDirichlet)
      dserror("condition with name dirichlet is not of type Dirichlet");
    flag = dirich->Get<vector<int> >("onoff");
  }
  const int owner = Owner();
  vector<int> dof = dis.Dof(this);
  for (unsigned j=0; j<dof.size(); ++j)
  {
    if (flag && (*flag)[j])
      lmdirich.push_back(1);
    else
      lmdirich.push_back(0);
    lmowner.push_back(owner);
    lm.push_back(dof[j]);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const Discretization& dis, vector<int>& lm, vector<int>& lmowner) const
{
  const int numnode = NumNode();
  const DRT::Node*const* nodes = Nodes();

  lm.clear();
  lmowner.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i=0; i<numnode; ++i)
    {
      const Node* node = nodes[i];
      dis.Dof(this,node,lm);
      lmowner.resize(lm.size(),node->Owner());
    }
  }

  // fill the vector with element dofs
  dis.Dof(this,lm);
  lmowner.resize(lm.size(),Owner());
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Element::Evaluate(ParameterList& params,
                           DRT::Discretization&      discretization,
                           vector<int>&              lm,
                           Epetra_SerialDenseMatrix& elemat1,
                           Epetra_SerialDenseMatrix& elemat2,
                           Epetra_SerialDenseVector& elevec1,
                           Epetra_SerialDenseVector& elevec2,
                           Epetra_SerialDenseVector& elevec3)
{
  cout << "DRT::Element::Evaluate:\n"
       << "Base class dummy routine DRT::Element::Evaluate(...) called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return -1;
}

#if 0 // this no longer is a dummy (but pure virtual) to check on the
      // parameter list. It can be a dummy again once everything with
      // EvaluateNeumann is fixed
/*----------------------------------------------------------------------*
 |  evaluate Neumann BC dummy (public)                       mwgee 01/07|
 *----------------------------------------------------------------------*/
int DRT::Element::EvaluateNeumann(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  DRT::Condition&           condition,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseVector& elevec1)
{
  cout << "DRT::Element::EvaluateNeumann:\n"
       << "Base class dummy routine DRT::Element::EvaluateNeumann(...) called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return -1;
}
#endif




#endif  // #ifdef CCADISCRET
