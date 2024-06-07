/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_bele3.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::Bele3Type Discret::ELEMENTS::Bele3Type::instance_;


Discret::ELEMENTS::Bele3Type& Discret::ELEMENTS::Bele3Type::Instance() { return instance_; }


Core::Communication::ParObject* Discret::ELEMENTS::Bele3Type::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Bele3* object = new Discret::ELEMENTS::Bele3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Bele3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  // Search for "BELE3". If found, search for "_"
  // the number after "_" is numdof: so BELE3_4 is a BELE3 element
  // with numdof=4
  std::size_t pos = eletype.rfind("BELE3");
  if (pos != std::string::npos)
  {
    if (eletype.substr(pos + 5, 1) == "_")
    {
      std::istringstream is(eletype.substr(pos + 6, 1));

      int numdof = -1;
      is >> numdof;
      Teuchos::RCP<Discret::ELEMENTS::Bele3> ele =
          Teuchos::rcp(new Discret::ELEMENTS::Bele3(id, owner));
      ele->set_num_dof_per_node(numdof);
      return ele;
    }
    else
    {
      FOUR_C_THROW("ERROR: Found BELE3 element without specified number of dofs!");
    }
  }

  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Bele3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(new Discret::ELEMENTS::Bele3(id, owner));
  return ele;
}


void Discret::ELEMENTS::Bele3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Bele3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::Bele3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs3 = definitions["BELE3_3"];

  defs3["TRI3"] =
      Input::LineDefinition::Builder().AddIntVector("TRI3", 3).AddOptionalNamedInt("MAT").Build();

  defs3["TRI6"] =
      Input::LineDefinition::Builder().AddIntVector("TRI6", 6).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD4"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD8"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD9"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddOptionalNamedInt("MAT").Build();

  std::map<std::string, Input::LineDefinition>& defs4 = definitions["BELE3_4"];

  defs4["TRI3"] =
      Input::LineDefinition::Builder().AddIntVector("TRI3", 3).AddOptionalNamedInt("MAT").Build();

  defs4["TRI6"] =
      Input::LineDefinition::Builder().AddIntVector("TRI6", 6).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD4"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD8"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD9"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddOptionalNamedInt("MAT").Build();
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Bele3LineType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Bele3Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Bele3::Bele3(int id, int owner)
    : Core::Elements::Element(id, owner), numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Bele3::Bele3(const Discret::ELEMENTS::Bele3& old)
    : Core::Elements::Element(old), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Bele3::Clone() const
{
  Discret::ELEMENTS::Bele3* newelement = new Discret::ELEMENTS::Bele3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Bele3::Shape() const
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // numdofpernode_
  AddtoPack(data, numdofpernode_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3::Print(std::ostream& os) const
{
  os << "Bele3_" << numdofpernode_ << " " << Core::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Bele3::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Bele3Line, Bele3>(
      Core::Communication::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Bele3::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}


Core::FE::GaussRule2D Discret::ELEMENTS::Bele3::get_optimal_gaussrule() const
{
  Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  switch (Shape())
  {
    case Core::FE::CellType::quad4:
      rule = Core::FE::GaussRule2D::quad_4point;
      break;
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Bele3::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // check if material is defined
  if (linedef->HaveNamed("MAT"))
  {
    int material = 0;
    // read number of material model
    linedef->ExtractInt("MAT", material);
    SetMaterial(0, Mat::Factory(material));
  }
  return true;
}

FOUR_C_NAMESPACE_CLOSE
