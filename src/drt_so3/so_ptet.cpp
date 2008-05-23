/*!----------------------------------------------------------------------**##
\file so_ptet.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_ptet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Ptet::Ptet(int id, int owner) :
DRT::Element(id,element_ptet,owner),
material_(0),
V_(-1.0),
nxyz_(0,0),
FisNew_(false),
F_(0,0)
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Ptet::Ptet(const DRT::ELEMENTS::Ptet& old) :
DRT::Element(old),
material_(old.material_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_),
V_(old.V_),
nxyz_(old.nxyz_),
FisNew_(old.FisNew_),
F_(old.F_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Ptet::~Ptet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::Pack(vector<char>& data) const
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
  // stresstype_
  AddtoPack(data,stresstype_);
  // V_
  AddtoPack(data,V_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::Unpack(const vector<char>& data)
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
  // stresstype_
  ExtractfromPack(position,data,stresstype_);
  // V_
  ExtractfromPack(position,data,V_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::so_tet4_expol(Epetra_SerialDenseMatrix& stresses,
                                        Epetra_SerialDenseMatrix& nodalstresses)
{
  static Epetra_SerialDenseMatrix expol(NUMNOD_PTET, NUMGPT_PTET);
  static bool isfilled;

  if (isfilled==true)
  {
    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);
  }
  else
  {
    expol(0,0)=1.0;
    expol(1,0)=1.0;
    expol(2,0)=1.0;
    expol(3,0)=1.0;

    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);

    isfilled=true;
  }
}




/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::Print(ostream& os) const
{
  os << "Ptet ";
  Element::Print(os);
  return;
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

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  gee 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Ptet::Volumes()
{
  dserror("volume not impl. yet");
  return NULL;
}


 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             gee 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Ptet::Surfaces()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)surfaces_.size()    == NumSurface() &&
      (int)surfaceptrs_.size() == NumSurface() &&
      dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(surfaceptrs_[0]) )
    return (DRT::Element**)(&(surfaceptrs_[0]));

  const int nsurf = NumSurface();
  surfaces_.resize(nsurf);
  surfaceptrs_.resize(nsurf);
  int nodeids[3];
  DRT::Node* nodes[3];

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[1];
  nodeids[2] = NodeIds()[3];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[1];
  nodes[2] = Nodes()[3];
  surfaces_[0] =
    rcp(new DRT::ELEMENTS::StructuralSurface(0,Owner(),3,nodeids,nodes,this,0));
  surfaceptrs_[0] = surfaces_[0].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[2];
  nodeids[2] = NodeIds()[3];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[2];
  nodes[2] = Nodes()[3];
  surfaces_[1] =
    rcp(new DRT::ELEMENTS::StructuralSurface(1,Owner(),3,nodeids,nodes,this,1));
  surfaceptrs_[1] = surfaces_[1].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[3];
  nodeids[2] = NodeIds()[2];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[3];
  nodes[2] = Nodes()[2];
  surfaces_[2] =
    rcp(new DRT::ELEMENTS::StructuralSurface(2,Owner(),3,nodeids,nodes,this,2));
  surfaceptrs_[2] = surfaces_[2].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[2];
  nodeids[2] = NodeIds()[1];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[2];
  nodes[2] = Nodes()[1];
  surfaces_[3] =
    rcp(new DRT::ELEMENTS::StructuralSurface(3,Owner(),3,nodeids,nodes,this,3));
  surfaceptrs_[3] = surfaces_[3].get();

  return (DRT::Element**)(&(surfaceptrs_[0]));

  return 0;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Ptet::Lines()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)lines_.size()    == NumLine() &&
      (int)lineptrs_.size() == NumLine() &&
      dynamic_cast<DRT::ELEMENTS::StructuralLine*>(lineptrs_[0]) )
    return (DRT::Element**)(&(lineptrs_[0]));

  const int nline = NumLine();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  int nodeids[2];
  DRT::Node* nodes[2];

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[1];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[1];
  lines_[0] =
    rcp(new DRT::ELEMENTS::StructuralLine(0,Owner(),2,nodeids,nodes,this,0));
  lineptrs_[0] = lines_[0].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[2];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[2];
  lines_[1] =
    rcp(new DRT::ELEMENTS::StructuralLine(1,Owner(),2,nodeids,nodes,this,1));
  lineptrs_[1] = lines_[1].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[2];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[2];
  lines_[2] =
    rcp(new DRT::ELEMENTS::StructuralLine(2,Owner(),2,nodeids,nodes,this,2));
  lineptrs_[2] = lines_[2].get();

  nodeids[0] = NodeIds()[0];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[0];
  nodes[1] = Nodes()[3];
  lines_[3] =
    rcp(new DRT::ELEMENTS::StructuralLine(3,Owner(),2,nodeids,nodes,this,3));
  lineptrs_[3] = lines_[3].get();

  nodeids[0] = NodeIds()[1];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[1];
  nodes[1] = Nodes()[3];
  lines_[4] =
    rcp(new DRT::ELEMENTS::StructuralLine(4,Owner(),2,nodeids,nodes,this,4));
  lineptrs_[4] = lines_[4].get();

  nodeids[0] = NodeIds()[2];
  nodeids[1] = NodeIds()[3];
  nodes[0] = Nodes()[2];
  nodes[1] = Nodes()[3];
  lines_[5] =
    rcp(new DRT::ELEMENTS::StructuralLine(5,Owner(),2,nodeids,nodes,this,5));
  lineptrs_[5] = lines_[5].get();

  return (DRT::Element**)(&(lineptrs_[0]));

  return 0;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PtetRegister::PtetRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PtetRegister::PtetRegister(const DRT::ELEMENTS::PtetRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PtetRegister::~PtetRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PtetRegister* DRT::ELEMENTS::PtetRegister::Clone() const
{
  return new DRT::ELEMENTS::PtetRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::Pack(vector<char>& data) const
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
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::Unpack(const vector<char>& data)
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
 |  print (public)                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PtetRegister::Print(ostream& os) const
{
  os << "PtetRegister ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
