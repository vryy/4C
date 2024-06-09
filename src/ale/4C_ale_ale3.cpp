/*----------------------------------------------------------------------------*/
/*! \file

\brief 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale3.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Ale3Type Discret::ELEMENTS::Ale3Type::instance_;

Discret::ELEMENTS::Ale3Type& Discret::ELEMENTS::Ale3Type::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::Ale3Type::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Ale3* object = new Discret::ELEMENTS::Ale3(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele;

  if (eletype == "ALE3")
  {
    if (eledistype != "NURBS27")
    {
      ele = Teuchos::rcp(new Discret::ELEMENTS::Ale3(id, owner));
    }
  }

  return ele;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale3Type::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Ale3(id, owner));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Ale3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["ALE3"];

  defs["HEX8"] =
      Input::LineDefinition::Builder().AddIntVector("HEX8", 8).AddNamedInt("MAT").Build();

  defs["HEX20"] =
      Input::LineDefinition::Builder().AddIntVector("HEX20", 20).AddNamedInt("MAT").Build();

  defs["HEX27"] =
      Input::LineDefinition::Builder().AddIntVector("HEX27", 27).AddNamedInt("MAT").Build();

  defs["TET4"] =
      Input::LineDefinition::Builder().AddIntVector("TET4", 4).AddNamedInt("MAT").Build();

  defs["TET10"] =
      Input::LineDefinition::Builder().AddIntVector("TET10", 10).AddNamedInt("MAT").Build();

  defs["WEDGE6"] =
      Input::LineDefinition::Builder().AddIntVector("WEDGE6", 6).AddNamedInt("MAT").Build();

  defs["WEDGE15"] =
      Input::LineDefinition::Builder().AddIntVector("WEDGE15", 15).AddNamedInt("MAT").Build();

  defs["PYRAMID5"] =
      Input::LineDefinition::Builder().AddIntVector("PYRAMID5", 5).AddNamedInt("MAT").Build();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale3SurfaceType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Ale3Surface( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale3::Ale3(int id, int owner) : Core::Elements::Element(id, owner) {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale3::Ale3(const Discret::ELEMENTS::Ale3& old) : Core::Elements::Element(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Ale3::Clone() const
{
  Discret::ELEMENTS::Ale3* newelement = new Discret::ELEMENTS::Ale3(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Ale3::Shape() const
{
  switch (num_node())
  {
    case 4:
      return Core::FE::CellType::tet4;
    case 5:
      return Core::FE::CellType::pyramid5;
    case 6:
      return Core::FE::CellType::wedge6;
    case 8:
      return Core::FE::CellType::hex8;
    case 10:
      return Core::FE::CellType::tet10;
    case 20:
      return Core::FE::CellType::hex20;
    case 27:
      return Core::FE::CellType::hex27;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3::Print(std::ostream& os) const
{
  os << "Ale3 ";
  Element::Print(os);
  std::cout << std::endl;
  // cout << data_;
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Ale3::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<Ale3Surface, Ale3>(
      Core::Communication::buildSurfaces, *this);
}

FOUR_C_NAMESPACE_CLOSE
