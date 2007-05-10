/*!----------------------------------------------------------------------
\file shell8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "shell8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(int id, int owner) :
DRT::Element(id,element_shell8,owner),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0),
data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(const DRT::Elements::Shell8& old) :
DRT::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_),
data_(old.data_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Shell8::Clone() const
{
  DRT::Elements::Shell8* newelement = new DRT::Elements::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Shell8::Shape() const
{
  switch (NumNode())
  {
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // forcetype_
  AddtoPack(data,forcetype_);
  // thickness_
  AddtoPack(data,thickness_);
  // ngp_
  AddtoPack(data,ngp_,3*sizeof(int));
  // ngptri_
  AddtoPack(data,ngptri_);
  // nhyb_
  AddtoPack(data,nhyb_);
  // eas_
  AddtoPack(data,eas_,5*sizeof(int));
  // ans_
  AddtoPack(data,ans_);
  // sdc_
  AddtoPack(data,sdc_);
  // material_
  AddtoPack(data,material_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // forcetype_
  ExtractfromPack(position,data,forcetype_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeof(int));
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::Print(ostream& os) const
{
  os << "Shell8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Shell8Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Shell8::ElementRegister() const
{
  return rcp(new DRT::Elements::Shell8Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Shell8::Lines()
{
  const int nline = NumLine();
  const int numnode = NumNode();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  int nodeids[100];
  DRT::Node* nodes[100];
  if (nline==4)
  {
    if (numnode==4)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      lines_[0] =
        rcp(new DRT::Elements::Shell8Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] =
        rcp(new DRT::Elements::Shell8Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      lines_[2] =
        rcp(new DRT::Elements::Shell8Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      lines_[3] =
        rcp(new DRT::Elements::Shell8Line(3,Owner(),2,nodeids,nodes,this,3));
      lineptrs_[3] = lines_[3].get();
    }
    else if (numnode==9)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[4];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[4];
      lines_[0] =
        rcp(new DRT::Elements::Shell8Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[5];
      lines_[1] =
        rcp(new DRT::Elements::Shell8Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodeids[2] = NodeIds()[6];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      nodes[2] = Nodes()[6];
      lines_[2] =
        rcp(new DRT::Elements::Shell8Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[7];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[7];
      lines_[3] =
        rcp(new DRT::Elements::Shell8Line(3,Owner(),3,nodeids,nodes,this,3));
      lineptrs_[3] = lines_[3].get();
    }
    else dserror("Number of nodes not supported");
  }
  else if (nline==3)
  {
    if (numnode==3)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      lines_[0] =
        rcp(new DRT::Elements::Shell8Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] =
        rcp(new DRT::Elements::Shell8Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[2] =
        rcp(new DRT::Elements::Shell8Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    else if (numnode==6)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[3];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[3];
      lines_[0] =
        rcp(new DRT::Elements::Shell8Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[4];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[4];
      lines_[1] =
        rcp(new DRT::Elements::Shell8Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[5];
      lines_[2] =
        rcp(new DRT::Elements::Shell8Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    else dserror("Number of nodes not supported");
  }
  else dserror("Number of lines not supported");
  return (DRT::Element**)(&(lineptrs_[0]));
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Shell8::Surfaces()
{
  surfaces_.resize(1);
  surfaces_[0] = this;
  return &surfaces_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::Shell8Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::Shell8Register(
                               const DRT::Elements::Shell8Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register* DRT::Elements::Shell8Register::Clone() const
{
  return new DRT::Elements::Shell8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::~Shell8Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Register::Print(ostream& os) const
{
  os << "Shell8Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
