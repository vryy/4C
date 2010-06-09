/*!----------------------------------------------------------------------**##
\file so_nstet.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nstet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"



DRT::ELEMENTS::NStetRegisterType DRT::ELEMENTS::NStetRegisterType::instance_;


DRT::ParObject* DRT::ELEMENTS::NStetRegisterType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NStetRegister* object =
    new DRT::ELEMENTS::NStetRegister(DRT::Element::element_nstet);
  object->Unpack(data);
  return object;
}


DRT::ELEMENTS::NStetType DRT::ELEMENTS::NStetType::instance_;


DRT::ParObject* DRT::ELEMENTS::NStetType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NStet* object = new DRT::ELEMENTS::NStet(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStetType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="NSTET4" )
  {
    RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NStet(id,owner));
    return ele;
  }
  return Teuchos::null;
}


void DRT::ELEMENTS::NStetType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::NStetType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}


/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(int id, int owner) :
DRT::Element(id,owner),
material_(0),
V_(-1.0),
nxyz_(),
FisNew_(false),
F_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(const DRT::ELEMENTS::NStet& old) :
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
 |  dtor (public)                                              gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::~NStet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Pack(vector<char>& data) const
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
  // stabtype_
  AddtoPack(data,stabtype_);
  // V_
  AddtoPack(data,V_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
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
  // stabtype_
  ExtractfromPack(position,data,stabtype_);
  // V_
  ExtractfromPack(position,data,V_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::so_nstet_expol(LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>& stresses,
                                        LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET>& nodalstresses)
{
  LINALG::Matrix<NUMNOD_NSTET, NUMGPT_NSTET> expol;
  expol(0,0)=1.0;
  expol(1,0)=1.0;
  expol(2,0)=1.0;
  expol(3,0)=1.0;
  nodalstresses.Multiply(expol,stresses);
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Print(ostream& os) const
{
  os << "NStet ";
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
 |  get vector of volumes (length 1) (public)                  gee 12/09|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Volumes()
{
  dserror("volume not impl. yet");
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             gee 12/09|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Surfaces()
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
 |  get vector of lines (public)                               gee 12/09|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Lines()
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
 |  ctor (public)                                              gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStetRegister::NStetRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStetRegister::NStetRegister(const DRT::ELEMENTS::NStetRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStetRegister::~NStetRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStetRegister* DRT::ELEMENTS::NStetRegister::Clone() const
{
  return new DRT::ELEMENTS::NStetRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetRegister::Pack(vector<char>& data) const
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
 |                                                             gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetRegister::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  print (public)                                             gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetRegister::Print(ostream& os) const
{
  os << "NStetRegister ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
