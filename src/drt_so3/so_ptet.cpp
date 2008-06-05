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
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Ptet::Ptet(const DRT::ELEMENTS::Ptet& old) :
DRT::Element(old),
material_(old.material_),
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::Ptet::Volumes()
{
  dserror("volume not impl. yet");
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             gee 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Ptet::Surfaces()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Ptet::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
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
