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
#include "ale3_nurbs.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


using namespace DRT::UTILS;

DRT::ELEMENTS::Ale3Type::Ale3Type DRT::ELEMENTS::Ale3Type::instance_;

DRT::ParObject* DRT::ELEMENTS::Ale3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Ale3* object = new DRT::ELEMENTS::Ale3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Ale3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  Teuchos::RCP<DRT::Element> ele;

  if ( eletype=="ALE3" )
  {
    if(eledistype!="NURBS27")
    {
      ele = rcp(new DRT::ELEMENTS::Ale3(id,owner));
    }
  }

  return ele;
}


DRT::ELEMENTS::Ale3::Ale3(int id, int owner)
  : DRT::Element(id,element_ale3,owner),
    data_()
{
}


DRT::ELEMENTS::Ale3::Ale3(const DRT::ELEMENTS::Ale3& old)
  : DRT::Element(old),
    data_(old.data_)
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
	vector<char>::size_type position = 0;
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

  if (position != data.size())
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
  return Teuchos::null;
}


//
// get vector of surfaces
//
vector<RCP<DRT::Element> > DRT::ELEMENTS::Ale3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Ale3Surface,Ale3>(DRT::UTILS::buildSurfaces,this);
}


vector<RCP<DRT::Element> > DRT::ELEMENTS::Ale3::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


#endif
#endif
