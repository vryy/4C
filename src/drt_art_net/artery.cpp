
/*! \file
\brief One-dimensional artery element

\level 3

\maintainer Johannes Kremheller

*----------------------------------------------------------------------*/

#include "artery.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;

DRT::ELEMENTS::ArteryType DRT::ELEMENTS::ArteryType::instance_;

DRT::ELEMENTS::ArteryType& DRT::ELEMENTS::ArteryType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::ArteryType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Artery* object = new DRT::ELEMENTS::Artery(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ArteryType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ART")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Artery(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ArteryType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Artery(id, owner));
  return ele;
}


void DRT::ELEMENTS::ArteryType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["ART"];

  defs["LINE2"]
      .AddIntVector("LINE2", 2)
      .AddNamedInt("MAT")
      .AddNamedInt("GP")
      .AddNamedString("TYPE")
      .AddNamedDouble("DIAM");

  defs["LIN2"]
      .AddIntVector("LIN2", 2)
      .AddNamedInt("MAT")
      .AddNamedInt("GP")
      .AddNamedString("TYPE")
      .AddNamedDouble("DIAM");
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::Artery(int id, int owner)
    : DRT::Element(id, owner), impltype_(INPAR::ARTDYN::impltype_undefined), data_()
{
  gaussrule_ = intrule1D_undefined;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::Artery(const DRT::ELEMENTS::Artery& old)
    : DRT::Element(old), impltype_(old.impltype_), gaussrule_(old.gaussrule_), data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Artery and return pointer to it (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Artery::Clone() const
{
  DRT::ELEMENTS::Artery* newelement = new DRT::ELEMENTS::Artery(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Artery::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);
  // Gaussrule
  AddtoPack(data, gaussrule_);  // implicit conversion from enum to integer
  AddtoPack(data, impltype_);

  // data_
  AddtoPack(data, data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);

  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // Gaussrule
  int gausrule_integer;
  ExtractfromPack(position, data, gausrule_integer);
  gaussrule_ = GaussRule1D(gausrule_integer);  // explicit conversion from integer to enum
  impltype_ = static_cast<INPAR::ARTDYN::ImplType>(ExtractInt(position, data));

  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       kremheller 10/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Artery::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                           ismail 01/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::~Artery() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Print(std::ostream& os) const
{
  os << "Artery ";
  Element::Print(os);

  return;
}
