/*!----------------------------------------------------------------------**##
\file so_tet4.cpp
\brief 

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet4.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::So_tet4(int id, int owner) :
DRT::Element(id,element_so_tet4,owner),
material_(0),
data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;  //whatis ngp_ ???????
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::So_tet4(const DRT::ELEMENTS::So_tet4& old) :
DRT::Element(old),
material_(old.material_),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_tet4::Clone() const
{
  DRT::ELEMENTS::So_tet4* newelement = new DRT::ELEMENTS::So_tet4(*this);  
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_tet4::Shape() const
{
  return tet4;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // ngp_
  //AddtoPack(data,ngp_,3*sizeof(int));
  // material_
  AddtoPack(data,material_);
  // stresstype_
  AddtoPack(data,stresstype_);
  // kintype_
  AddtoPack(data,kintype_);
  // rewind flags
  AddtoPack(data,donerewinding_);

  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Unpack(const vector<char>& data)
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
  // ngp_
  //ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // material_
  ExtractfromPack(position,data,material_);
  // stresstype_
  ExtractfromPack(position,data,stresstype_);
  // kintype_
  ExtractfromPack(position,data,kintype_);
  // rewinding flags
  ExtractfromPack(position,data,donerewinding_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::~So_tet4()
{
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Print(ostream& os) const
{
  os << "So_tet4 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*------------------------------------------------------------------------*
 |  allocate and return So_tet4Register (public)               volf 10/07|
 *------------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_tet4::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Sotet4Register(Type()));
}

  /*====================================================================*/
  /* 4-node tetrahedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
   * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
   *  4-node hexahedron: node 0,1,...,3
   *          
   * -----------------------
   *- this is the numbering used in GiD & EXODUS!!
   *      3-
   *      |\ ---
   *      |  \    ---
   *      |    \      ---
   *      |      \        -2
   *      |        \       /\
   *      |          \   /   \
   *      |            X      \
   *      |          /   \     \
   *      |        /       \    \
   *      |      /           \   \
   *      |    /               \  \
   *      |  /                   \ \
   *      |/                       \\
   *      0--------------------------1
   */ 
  /*====================================================================*/

/*----------------------------------------------------------------------***
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::So_tet4::Volumes()
{
  volume_.resize(1);
  return 0;
}


 /*----------------------------------------------------------------------**#
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::So_tet4::Surfaces()
{
  
  const int nsurf = NumSurface();
  surfaces_.resize(nsurf);
  surfaceptrs_.resize(nsurf);
  int nodeids[100];
  DRT::Node* nodes[100];

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[1];
  nodeids[2] = NodeIds()[3];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[1];
  nodes[2] = Nodes()[3];
  surfaces_[0] =
    rcp(new DRT::ELEMENTS::Sotet4Surface(0,Owner(),3,nodeids,nodes,this,0));
  surfaceptrs_[0] = surfaces_[0].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[2];
  nodeids[2] = NodeIds()[3];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[2];
  nodes[2] = Nodes()[3];
  surfaces_[1] =
    rcp(new DRT::ELEMENTS::Sotet4Surface(1,Owner(),3,nodeids,nodes,this,1));
  surfaceptrs_[1] = surfaces_[1].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[3];
  nodeids[2] = NodeIds()[2];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[3];
  nodes[2] = Nodes()[2];
  surfaces_[2] =
    rcp(new DRT::ELEMENTS::Sotet4Surface(2,Owner(),3,nodeids,nodes,this,2));
  surfaceptrs_[2] = surfaces_[2].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[2];
  nodeids[2] = NodeIds()[1];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[2];
  nodes[2] = Nodes()[1];
  surfaces_[3] =
    rcp(new DRT::ELEMENTS::Sotet4Surface(3,Owner(),3,nodeids,nodes,this,3));
  surfaceptrs_[3] = surfaces_[3].get();

  return (DRT::Element**)(&(surfaceptrs_[0]));
  
  return 0;
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::So_tet4::Lines()
{
  const int nline = NumLine();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  int nodeids[100];
  DRT::Node* nodes[100];

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[1];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[1];
  lines_[0] =
    rcp(new DRT::ELEMENTS::Sotet4Line(0,Owner(),2,nodeids,nodes,this,0));
  lineptrs_[0] = lines_[0].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[2];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[2];
  lines_[1] =
    rcp(new DRT::ELEMENTS::Sotet4Line(1,Owner(),2,nodeids,nodes,this,1));
  lineptrs_[1] = lines_[1].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[2];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[2];
  lines_[2] =
    rcp(new DRT::ELEMENTS::Sotet4Line(2,Owner(),2,nodeids,nodes,this,2));
  lineptrs_[2] = lines_[2].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[3];
  lines_[3] =
    rcp(new DRT::ELEMENTS::Sotet4Line(3,Owner(),2,nodeids,nodes,this,3));
  lineptrs_[3] = lines_[3].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[3];
  lines_[4] =
    rcp(new DRT::ELEMENTS::Sotet4Line(4,Owner(),2,nodeids,nodes,this,4));
  lineptrs_[4] = lines_[4].get();

  nodeids[0] = NodeIds()[2];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[2];
  nodes[1] = Nodes()[3];
  lines_[5] =
    rcp(new DRT::ELEMENTS::Sotet4Line(5,Owner(),2,nodeids,nodes,this,5));
  lineptrs_[5] = lines_[5].get();

  return (DRT::Element**)(&(lineptrs_[0]));

  return 0;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sotet4Register::Sotet4Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sotet4Register::Sotet4Register(
                               const DRT::ELEMENTS::Sotet4Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sotet4Register* DRT::ELEMENTS::Sotet4Register::Clone() const
{
//  return new DRT::ELEMENTS::Soh8Register(*this);
  return new DRT::ELEMENTS::Sotet4Register(*this);
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sotet4Register::Pack(vector<char>& data) const
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


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
//void DRT::ELEMENTS::Soh8Register::Unpack(const vector<char>& data)
void DRT::ELEMENTS::Sotet4Register::Unpack(const vector<char>& data)
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


/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sotet4Register::~Sotet4Register()
{
  return;
}

/*----------------------------------------------------------------------***
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sotet4Register::Print(ostream& os) const
{
  os << "Sotet4Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
