/*!----------------------------------------------------------------------
\file vele3.cpp
\brief

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "vele3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3::Vele3(int id, int owner) :
DRT::Element(id,element_vele3,owner)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3::Vele3(const DRT::ELEMENTS::Vele3& old) :
DRT::Element(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Vele3::Clone() const
{
  DRT::ELEMENTS::Vele3* newelement = new DRT::ELEMENTS::Vele3(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Vele3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);

  Element::Pack(basedata);

  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Unpack(const vector<char>& data)
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

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3::~Vele3()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Print(ostream& os) const
{
  os << "Vele3 " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Vele3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Vele3Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Vele3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Vele3Line,Vele3>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Vele3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Vele3Surface,Vele3>(DRT::UTILS::buildSurfaces,this);
}



/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Vele3::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule (public)                          u.may 05/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule3D DRT::ELEMENTS::Vele3::getOptimalGaussrule(const DRT::Element::DiscretizationType& distype) const
{
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
      case DRT::Element::hex8:
          rule = DRT::UTILS::intrule_hex_8point;
          break;
      case DRT::Element::hex20: case DRT::Element::hex27:
          rule = DRT::UTILS::intrule_hex_27point;
          break;
      case DRT::Element::tet4:
          rule = DRT::UTILS::intrule_tet_4point;
          break;
      case DRT::Element::tet10:
          rule = DRT::UTILS::intrule_tet_10point;
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
DRT::ELEMENTS::Vele3Register::Vele3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Register::Vele3Register(
                               const DRT::ELEMENTS::Vele3Register& old) :
ElementRegister(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Register* DRT::ELEMENTS::Vele3Register::Clone() const
{
  return new DRT::ELEMENTS::Vele3Register(*this);
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Vele3Register::Unpack(const vector<char>& data)
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
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Register::~Vele3Register()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Register::Print(ostream& os) const
{
  os << "Vele3Register ";
  ElementRegister::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
