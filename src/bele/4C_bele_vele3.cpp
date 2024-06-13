/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_vele3.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Vele3Type Discret::ELEMENTS::Vele3Type::instance_;

Discret::ELEMENTS::Vele3Type& Discret::ELEMENTS::Vele3Type::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::Vele3Type::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Vele3* object = new Discret::ELEMENTS::Vele3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Vele3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "VELE3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Vele3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Vele3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(new Discret::ELEMENTS::Vele3(id, owner));
  return ele;
}


void Discret::ELEMENTS::Vele3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Vele3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented for element type vele3!");
  return nullspace;
}

void Discret::ELEMENTS::Vele3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["VELE3"];

  defs["HEX8"] = Input::LineDefinition::Builder().add_int_vector("HEX8", 8).Build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Vele3SurfaceType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Vele3Surface( id, owner ) );
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Vele3LineType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Vele3Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3::Vele3(int id, int owner) : Core::Elements::Element(id, owner) { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3::Vele3(const Discret::ELEMENTS::Vele3& old) : Core::Elements::Element(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Vele3::Clone() const
{
  Discret::ELEMENTS::Vele3* newelement = new Discret::ELEMENTS::Vele3(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Vele3::Shape() const
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
    case 15:
      return Core::FE::CellType::wedge15;
    case 20:
      return Core::FE::CellType::hex20;
    case 27:
      return Core::FE::CellType::hex27;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3::Print(std::ostream& os) const
{
  os << "Vele3 " << Core::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Vele3::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Vele3Line, Vele3>(
      Core::Communication::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Vele3::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<Vele3Surface, Vele3>(
      Core::Communication::buildSurfaces, *this);
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule (public)                          u.may 05/09|
 *----------------------------------------------------------------------*/
Core::FE::GaussRule3D Discret::ELEMENTS::Vele3::get_optimal_gaussrule(
    const Core::FE::CellType& distype) const
{
  Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      rule = Core::FE::GaussRule3D::hex_8point;
      break;
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
      rule = Core::FE::GaussRule3D::hex_27point;
      break;
    case Core::FE::CellType::tet4:
      rule = Core::FE::GaussRule3D::tet_4point;
      break;
    case Core::FE::CellType::tet10:
      rule = Core::FE::GaussRule3D::tet_10point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Vele3::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  return true;
}

FOUR_C_NAMESPACE_CLOSE
