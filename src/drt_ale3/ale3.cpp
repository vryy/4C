//-----------------------------------------------------------------------
/*!
\file ale3.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE
#ifdef CCADISCRET

#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


using namespace DRT::UTILS;

DRT::ELEMENTS::Ale3::Ale3(int id, int owner)
  : DRT::Element(id,element_ale3,owner),
    data_()
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
}


DRT::ELEMENTS::Ale3::Ale3(const DRT::ELEMENTS::Ale3& old)
  : DRT::Element(old),
    data_(old.data_),
    surfaces_(old.surfaces_),
    surfaceptrs_(old.surfaceptrs_)
{
  return;
}


DRT::Element* DRT::ELEMENTS::Ale3::Clone() const
{
  DRT::ELEMENTS::Ale3* newelement = new DRT::ELEMENTS::Ale3(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::ELEMENTS::Ale3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


void DRT::ELEMENTS::Ale3::Pack(vector<char>& data) const
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
  //AddtoPack(data,gaussrule_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);
}


void DRT::ELEMENTS::Ale3::Unpack(const vector<char>& data)
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
  // Gaussrule
  //ExtractfromPack(position,data,gaussrule_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


DRT::ELEMENTS::Ale3::~Ale3()
{
}


void DRT::ELEMENTS::Ale3::Print(ostream& os) const
{
  os << "Ale3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Ale3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Ale3Register(Type()));
}


//
// get vector of surfaces
//
DRT::Element** DRT::ELEMENTS::Ale3::Surfaces()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)surfaces_.size()    == NumSurface() &&
      (int)surfaceptrs_.size() == NumSurface() &&
      dynamic_cast<DRT::ELEMENTS::Ale3Surface*>(surfaceptrs_[0]) )
    return (DRT::Element**)(&(surfaceptrs_[0]));

  // so we have to allocate new surface elements
  DRT::UTILS::ElementBoundaryFactory<Ale3Surface,Ale3>(false,surfaces_,surfaceptrs_,this);
  
  return (DRT::Element**)(&(surfaceptrs_[0]));
}


DRT::Element** DRT::ELEMENTS::Ale3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Ale3 element itself
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


DRT::ELEMENTS::Ale3Register::Ale3Register(DRT::Element::ElementType etype)
  : ElementRegister(etype)
{
}


DRT::ELEMENTS::Ale3Register::Ale3Register(const DRT::ELEMENTS::Ale3Register& old)
  : ElementRegister(old)
{
}


DRT::ELEMENTS::Ale3Register* DRT::ELEMENTS::Ale3Register::Clone() const
{
  return new DRT::ELEMENTS::Ale3Register(*this);
}


void DRT::ELEMENTS::Ale3Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);
}


void DRT::ELEMENTS::Ale3Register::Unpack(const vector<char>& data)
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
}


DRT::ELEMENTS::Ale3Register::~Ale3Register()
{
}


void DRT::ELEMENTS::Ale3Register::Print(ostream& os) const
{
  os << "Ale3Register ";
  ElementRegister::Print(os);
}


#endif
#endif
