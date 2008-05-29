/*!----------------------------------------------------------------------
\file so_disp.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(int id, int owner) :
DRT::Element(id,element_sodisp,owner),
data_()
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(const DRT::ELEMENTS::SoDisp& old) :
DRT::Element(old),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoDisp::Clone() const
{
  DRT::ELEMENTS::SoDisp* newelement = new DRT::ELEMENTS::SoDisp(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::SoDisp::Shape() const
{
    switch (NumNode())
    {
    case  4: return tet4;
    case  8: return hex8;
    case 10: return tet10;
    case 20: return hex20;
    case 27: return hex27;
    case  6: return wedge6;
    case 15: return wedge15;
    case  5: return pyramid5;
    default:
      dserror("unexpected number of nodes %d", NumNode());
    }
    return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // stresstype_
  AddtoPack(data,stresstype_);
  // kintype_
  AddtoPack(data,kintype_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Unpack(const vector<char>& data)
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


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::~SoDisp()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Print(ostream& os) const
{
  os << "SoDisp ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return SoDispRegister (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::SoDisp::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::SoDispRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Lines()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)lines_.size()    == NumLine() &&
      (int)lineptrs_.size() == NumLine() &&
      dynamic_cast<DRT::ELEMENTS::SoDispLine*>(lineptrs_[0]) )
    return (DRT::Element**)(&(lineptrs_[0]));
  
  // so we have to allocate new line elements
  DRT::UTILS::ElementBoundaryFactory<SoDispLine,SoDisp>(DRT::UTILS::buildLines,lines_,lineptrs_,this);

  return (DRT::Element**)(&(lineptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Surfaces()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)surfaces_.size()    == NumSurface() &&
      (int)surfaceptrs_.size() == NumSurface() &&
      dynamic_cast<DRT::ELEMENTS::SoDispSurface*>(surfaceptrs_[0]) )
    return (DRT::Element**)(&(surfaceptrs_[0]));

  // so we have to allocate new surface elements
  DRT::UTILS::ElementBoundaryFactory<SoDispSurface,SoDisp>(DRT::UTILS::buildSurfaces,surfaces_,surfaceptrs_,this);
  
  return (DRT::Element**)(&(surfaceptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Volumes()
{
  volume_.resize(1);
  volume_[0] = this;
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister::SoDispRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister::SoDispRegister(
                               const DRT::ELEMENTS::SoDispRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister* DRT::ELEMENTS::SoDispRegister::Clone() const
{
  return new DRT::ELEMENTS::SoDispRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispRegister::Pack(vector<char>& data) const
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
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispRegister::Unpack(const vector<char>& data)
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
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister::~SoDispRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispRegister::Print(ostream& os) const
{
  os << "SoDispRegister ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                      12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDispRegister::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
