/*!----------------------------------------------------------------------
\file wall1.cpp
\brief

<pre>
Maintainer: Maarkus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "wall1.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::Wall1(int id, int owner) :
DRT::Element(id,element_wall1,owner),
material_(0)
{
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::Wall1(const DRT::Elements::Wall1& old) :
DRT::Element(old),
material_(old.material_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Wall1::Clone() const
{
  DRT::Elements::Wall1* newelement = new DRT::Elements::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Wall1::Shape() const
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
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // material_
  AddtoPack(data,material_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Unpack(const vector<char>& data)
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
  // material_
  ExtractfromPack(position,data,material_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::~Wall1()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Print(ostream& os) const
{
  os << "Wall1 ";
  Element::Print(os);
  os << endl;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Wall1Register (public)              mgit 03/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Wall1::ElementRegister() const
{
  return rcp(new DRT::Elements::Wall1Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Wall1::Lines()
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
        rcp(new DRT::Elements::Wall1Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] = 
        rcp(new DRT::Elements::Wall1Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      lines_[2] = 
        rcp(new DRT::Elements::Wall1Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      lines_[3] = 
        rcp(new DRT::Elements::Wall1Line(3,Owner(),2,nodeids,nodes,this,3));
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
        rcp(new DRT::Elements::Wall1Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[5];
      lines_[1] = 
        rcp(new DRT::Elements::Wall1Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodeids[2] = NodeIds()[6];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      nodes[2] = Nodes()[6];
      lines_[2] = 
        rcp(new DRT::Elements::Wall1Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[7];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[7];
      lines_[3] = 
        rcp(new DRT::Elements::Wall1Line(3,Owner(),3,nodeids,nodes,this,3));
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
        rcp(new DRT::Elements::Wall1Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] = 
        rcp(new DRT::Elements::Wall1Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[2] = 
        rcp(new DRT::Elements::Wall1Line(2,Owner(),2,nodeids,nodes,this,2));
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
        rcp(new DRT::Elements::Wall1Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[4];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[4];
      lines_[1] = 
        rcp(new DRT::Elements::Wall1Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[5];
      lines_[2] = 
        rcp(new DRT::Elements::Wall1Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    else dserror("Number of nodes not supported");
  }
  else dserror("Number of lines not supported");
  return (DRT::Element**)(&(lineptrs_[0]));
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Wall1::Surfaces()
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
 |  ctor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::Wall1Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::Wall1Register(
                               const DRT::Elements::Wall1Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register* DRT::Elements::Wall1Register::Clone() const
{
  return new DRT::Elements::Wall1Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Pack(vector<char>& data) const
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
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::~Wall1Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Print(ostream& os) const
{
  os << "Wall1Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::Elements::Wall1Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
