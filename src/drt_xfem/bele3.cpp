/*!----------------------------------------------------------------------
\file bele3.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "bele3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(int id, int owner) :
DRT::Element(id,element_bele3,owner),
data_(),
is_moving_(true)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(const DRT::ELEMENTS::Bele3& old) :
DRT::Element(old),
data_(old.data_),
is_moving_(old.is_moving_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Bele3 and return pointer to it (public) |
 |                                                          gammi 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele3::Clone() const
{
  DRT::ELEMENTS::Bele3* newelement = new DRT::ELEMENTS::Bele3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele3::Shape() const
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
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Unpack(const vector<char>& data)
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

  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::~Bele3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Print(ostream& os) const
{
  os << "Bele3 ";
  Element::Print(os);
  os << endl;
  os << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Bele3Register (public)              gammi 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Bele3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Bele3Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Bele3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Bele3Line,Bele3>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Bele3::Surfaces()
{
  vector<RCP<DRT::Element> > surfaces(1);
  surfaces[0]=rcp(this,false);
  return surfaces;
}


DRT::UTILS::GaussRule2D DRT::ELEMENTS::Bele3::getOptimalGaussrule(const DRT::Element::DiscretizationType& distype)
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
    switch (distype)
    {
    case DRT::Element::quad4:
        rule = DRT::UTILS::intrule_quad_4point;
        break;
    case DRT::Element::quad8: case DRT::Element::quad9:
        rule = DRT::UTILS::intrule_quad_9point;
        break;
    case DRT::Element::tri3:
        rule = DRT::UTILS::intrule_tri_3point;
        break;
    case DRT::Element::tri6:
        rule = DRT::UTILS::intrule_tri_6point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register::Bele3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register::Bele3Register(
                               const DRT::ELEMENTS::Bele3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register* DRT::ELEMENTS::Bele3Register::Clone() const
{
  return new DRT::ELEMENTS::Bele3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Bele3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Bele3Register::~Bele3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Register::Print(ostream& os) const
{
  os << "Bele3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
