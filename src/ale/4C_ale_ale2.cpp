/*----------------------------------------------------------------------------*/
/*! \file

\brief ALE element for 2D case


\level 1
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale2.hpp"

#include "4C_ale_ale2_nurbs.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Ale2Type Discret::ELEMENTS::Ale2Type::instance_;

Discret::ELEMENTS::Ale2Type& Discret::ELEMENTS::Ale2Type::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::Ale2Type::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Ale2* object = new Discret::ELEMENTS::Ale2(-1, -1);
  object->unpack(data);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale2Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele;

  if (eletype == "ALE2")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      ele = Teuchos::rcp(new Discret::ELEMENTS::Ale2(id, owner));
    }
  }

  return ele;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale2Type::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Ale2(id, owner));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Ale2Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["ALE2"];

  defs["QUAD4"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD4", 4).add_named_int("MAT").build();

  defs["QUAD8"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD8", 8).add_named_int("MAT").build();

  defs["QUAD9"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD9", 9).add_named_int("MAT").build();

  defs["TRI3"] =
      Input::LineDefinition::Builder().add_int_vector("TRI3", 3).add_named_int("MAT").build();

  defs["TRI6"] =
      Input::LineDefinition::Builder().add_int_vector("TRI6", 6).add_named_int("MAT").build();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Ale2LineType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Ale2Line( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale2::Ale2(int id, int owner) : Core::Elements::Element(id, owner) {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale2::Ale2(const Discret::ELEMENTS::Ale2& old) : Core::Elements::Element(old) {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Ale2::Clone() const
{
  Discret::ELEMENTS::Ale2* newelement = new Discret::ELEMENTS::Ale2(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Ale2::Shape() const
{
  switch (num_node())
  {
    case 3:
      return Core::FE::CellType::tri3;
    case 4:
      return Core::FE::CellType::quad4;
    case 6:
      return Core::FE::CellType::tri6;
    case 8:
      return Core::FE::CellType::quad8;
    case 9:
      return Core::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2::Print(std::ostream& os) const
{
  os << "Ale2 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Ale2::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Ale2Line, Ale2>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Ale2::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::GaussRule2D Discret::ELEMENTS::Ale2::get_optimal_gaussrule(
    const Core::FE::CellType& distype)
{
  Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::nurbs4:
      rule = Core::FE::GaussRule2D::quad_4point;
      break;
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::nurbs9:
      rule = Core::FE::GaussRule2D::quad_9point;
      break;
    case Core::FE::CellType::tri3:
      rule = Core::FE::GaussRule2D::tri_3point;
      break;
    case Core::FE::CellType::tri6:
      rule = Core::FE::GaussRule2D::tri_6point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
      break;
  }
  return rule;
}

FOUR_C_NAMESPACE_CLOSE
