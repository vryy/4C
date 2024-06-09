
/*! \file
\brief One-dimensional artery element

\level 3


*----------------------------------------------------------------------*/

#include "4C_art_net_artery.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::ArteryType Discret::ELEMENTS::ArteryType::instance_;

Discret::ELEMENTS::ArteryType& Discret::ELEMENTS::ArteryType::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::ArteryType::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Artery* object = new Discret::ELEMENTS::Artery(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ArteryType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ART")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Artery(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ArteryType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Artery(id, owner));
  return ele;
}


void Discret::ELEMENTS::ArteryType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["ART"];

  defs["LINE2"] = Input::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedInt("GP")
                      .AddNamedString("TYPE")
                      .AddNamedDouble("DIAM")
                      .Build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Artery::Artery(int id, int owner)
    : Core::Elements::Element(id, owner), impltype_(Inpar::ArtDyn::impltype_undefined)
{
  gaussrule_ = Core::FE::GaussRule1D::undefined;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Artery::Artery(const Discret::ELEMENTS::Artery& old)
    : Core::Elements::Element(old), impltype_(old.impltype_), gaussrule_(old.gaussrule_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Artery and return pointer to it (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Artery::Clone() const
{
  Discret::ELEMENTS::Artery* newelement = new Discret::ELEMENTS::Artery(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Artery::Shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Artery::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  Element::Pack(data);
  // Gaussrule
  add_to_pack(data, gaussrule_);
  add_to_pack(data, impltype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Artery::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);
  // Gaussrule
  extract_from_pack(position, data, gaussrule_);
  impltype_ = static_cast<Inpar::ArtDyn::ImplType>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       kremheller 10/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Artery::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/09|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Artery::Print(std::ostream& os) const
{
  os << "Artery ";
  Element::Print(os);

  return;
}

FOUR_C_NAMESPACE_CLOSE
