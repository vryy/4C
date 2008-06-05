/*!----------------------------------------------------------------------**##
\file so_ctet10.cpp
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

#include "so_ctet10.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_ctet10::So_ctet10(int id, int owner) :
DRT::Element(id,element_so_ctet10,owner),
material_(0),
data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0; //whatis ngp_ ???????
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_ctet10::So_ctet10(const DRT::ELEMENTS::So_ctet10& old) :
DRT::Element(old),
material_(old.material_),
data_(old.data_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_ctet10::Clone() const
{
  DRT::ELEMENTS::So_ctet10* newelement = new DRT::ELEMENTS::So_ctet10(*this);  
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_ctet10::Shape() const
{
  return tet10;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::Pack(vector<char>& data) const
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

  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::So_ctet10::~So_ctet10()
{
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::Print(ostream& os) const
{
  os << "So_ctet10 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------***
 |  extrapolate stresses to the nodes   (public)               vlf 04/08|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::So_ctet10::so_ctet10_expol(Epetra_SerialDenseMatrix& stresses,
                                        Epetra_SerialDenseMatrix& nodalstresses)
{
   //since stresses are evaluated at nodes, they dont need to be extrapolated
   nodalstresses = stresses;
}

/*-------------------------------------------------------------------------*
 |  allocate and return So_ctet10Register (public)               volf 06/07|
 *-------------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_ctet10::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Soctet10Register(Type()));
}

  /*====================================================================*/
  /* 10-node tetrahedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
   * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
   *  10-node hexahedron: node 0,1,...,9
   *          
   * -----------------------
   *- this is the numbering used in GiD & EXODUS!!
   *      3-
   *      |\ ---
   *      |  \    --9
   *      |    \      ---
   *      |      \        -2
   *      |        \       /\
   *      |          \   /   \
   *      7            8      \
   *      |          /   \     \
   *      |        6       \    5
   *      |      /           \   \
   *      |    /               \  \
   *      |  /                   \ \
   *      |/                       \\
   *      0------------4-------------1
   */ 
  /*====================================================================*/

/*----------------------------------------------------------------------***
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_ctet10::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------**#
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_ctet10::Surfaces()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Soctet10Surface,So_ctet10>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_ctet10::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Soctet10Line,So_ctet10>(DRT::UTILS::buildLines,this);
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soctet10Register::Soctet10Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soctet10Register::Soctet10Register(
                               const DRT::ELEMENTS::Soctet10Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soctet10Register* DRT::ELEMENTS::Soctet10Register::Clone() const
{
//  return new DRT::ELEMENTS::Soh8Register(*this);
  return new DRT::ELEMENTS::Soctet10Register(*this);
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soctet10Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Soctet10Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Soctet10Register::~Soctet10Register()
{
  return;
}

/*----------------------------------------------------------------------***
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soctet10Register::Print(ostream& os) const
{
  os << "Soctet10Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
