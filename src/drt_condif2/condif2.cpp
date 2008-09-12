/*!----------------------------------------------------------------------
\file condif2.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "condif2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/matlist.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 05/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2::Condif2(int id, int owner) :
DRT::Element(id,element_condif2,owner),
gaussrule_(intrule2D_undefined),
data_(),
numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          vg 05/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2::Condif2(const DRT::ELEMENTS::Condif2& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Condif2 and return pointer to it (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Condif2::Clone() const
{
  DRT::ELEMENTS::Condif2* newelement = new DRT::ELEMENTS::Condif2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::SetMaterial(int matnum)
{
  // the standard part:
  //mat_ = MAT::Material::Factory(matnum);  // not allowed since mat_ is private
  DRT::Element::SetMaterial(matnum);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  RefCountPtr<MAT::Material> mat = Material();
  MATERIAL* actmat = NULL;
  if(mat->MaterialType()== m_condif)
  {
    numdofpernode_=1; // we only have a single scalar
  }
  else if (mat->MaterialType()== m_matlist) // we have a system of scalars
  {
    actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
    numdofpernode_=actmat->m.matlist->nummat;
  }
  else
    dserror("condif material expected but got type %d", mat->MaterialType());

  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Condif2::Shape() const
{
  switch (NumNode())
  {
  case  3: return tri3;
  case  4: return quad4;
  case  6: return tri6;
  case  8: return quad8;
  case  9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // Gaussrule
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  // numdofpernode
  AddtoPack(data,numdofpernode_);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // Gaussrule
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = GaussRule2D(gausrule_integer); //explicit conversion from integer to enum
  // numdofpernode
  ExtractfromPack(position,data,numdofpernode_);

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2::~Condif2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                 vg 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2::Print(ostream& os) const
{
  os << "Condif2 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Condif2Register (public)                vg 05/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Condif2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Condif2Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Condif2::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Condif2Line,Condif2>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)                  vg 05/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Condif2::Surfaces()
{
  vector<RCP<Element> > surfaces(1);
  surfaces[0]= rcp(this, false);
  return surfaces;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2Register::Condif2Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          vg 05/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2Register::Condif2Register(
                               const DRT::ELEMENTS::Condif2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2Register* DRT::ELEMENTS::Condif2Register::Clone() const
{
  return new DRT::ELEMENTS::Condif2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2Register::Pack(vector<char>& data) const
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
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2Register::~Condif2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                              vg 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif2Register::Print(ostream& os) const
{
  os << "Condif2Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
