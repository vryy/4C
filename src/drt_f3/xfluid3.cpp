/*!----------------------------------------------------------------------
\file xfluid3.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*/
// map to convert strings tao actions (stabilisation)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::XFluid3::StabilisationAction> DRT::ELEMENTS::XFluid3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3::XFluid3(int id, int owner) :
DRT::Element(id,element_xfluid3,owner),
is_ale_(false),
data_(),
eleDofManager_()
{
    surfaces_.resize(0);
    surfaceptrs_.resize(0);
    lines_.resize(0);
    lineptrs_.resize(0);
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3::XFluid3(const DRT::ELEMENTS::XFluid3& old) :
DRT::Element(old),
is_ale_(old.is_ale_),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_),
eleDofManager_(old.eleDofManager_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of XFluid3 and return pointer to it (public)|
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::XFluid3::Clone() const
{
  DRT::ELEMENTS::XFluid3* newelement = new DRT::ELEMENTS::XFluid3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::XFluid3::Shape() const
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
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // is_ale_
  AddtoPack(data,is_ale_);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::Unpack(const vector<char>& data)
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
  // is_ale_
  ExtractfromPack(position,data,is_ale_);

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3::~XFluid3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3::Print(ostream& os) const
{
  os << "XFluid3 ";
  Element::Print(os);
  std::cout << endl;
  std::cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 02/08|
 *----------------------------------------------------------------------*/
RCP<DRT::ElementRegister> DRT::ELEMENTS::XFluid3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::XFluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::XFluid3::Lines()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)lines_.size()    == NumLine() &&
      (int)lineptrs_.size() == NumLine() &&
      dynamic_cast<DRT::ELEMENTS::XFluid3Line*>(lineptrs_[0]) )
    return (DRT::Element**)(&(lineptrs_[0]));
  
  // so we have to allocate new line elements
  DRT::UTILS::ElementBoundaryFactory<XFluid3Line,XFluid3>(DRT::UTILS::buildLines,lines_,lineptrs_,this);

  return (DRT::Element**)(&(lineptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::XFluid3::Surfaces()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)surfaces_.size()    == NumSurface() &&
      (int)surfaceptrs_.size() == NumSurface() &&
      dynamic_cast<DRT::ELEMENTS::XFluid3Surface*>(surfaceptrs_[0]) )
    return (DRT::Element**)(&(surfaceptrs_[0]));

  // so we have to allocate new surface elements
  DRT::UTILS::ElementBoundaryFactory<XFluid3Surface,XFluid3>(DRT::UTILS::buildSurfaces,surfaces_,surfaceptrs_,this);
  
  return (DRT::Element**)(&(surfaceptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::XFluid3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Fluid3 element itself
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Register::XFluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Register::XFluid3Register(
                               const DRT::ELEMENTS::XFluid3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Register* DRT::ELEMENTS::XFluid3Register::Clone() const
{
  return new DRT::ELEMENTS::XFluid3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::XFluid3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::XFluid3Register::~XFluid3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3Register::Print(ostream& os) const
{
  os << "XFluid3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
