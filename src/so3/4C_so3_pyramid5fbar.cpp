/*----------------------------------------------------------------------*/
/*! \file

\brief pyramid shaped solid element

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_pyramid5fbar.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::SoPyramid5fbarType Discret::ELEMENTS::SoPyramid5fbarType::instance_;


Discret::ELEMENTS::SoPyramid5fbarType& Discret::ELEMENTS::SoPyramid5fbarType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoPyramid5fbarType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoPyramid5fbar(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoPyramid5fbarType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoPyramid5fbar(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoPyramid5fbarType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoPyramid5fbar(id, owner));
  return ele;
}


void Discret::ELEMENTS::SoPyramid5fbarType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoPyramid5fbarType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::SoPyramid5fbarType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["PYRAMID5"] = Input::LineDefinition::Builder()
                         .add_int_vector("PYRAMID5", 5)
                         .add_named_int("MAT")
                         .add_named_string("KINEM")
                         .add_optional_named_double_vector("RAD", 3)
                         .add_optional_named_double_vector("AXI", 3)
                         .add_optional_named_double_vector("CIR", 3)
                         .add_optional_named_double_vector("FIBER1", 3)
                         .add_optional_named_double_vector("FIBER2", 3)
                         .add_optional_named_double_vector("FIBER3", 3)
                         .add_optional_named_double("GROWTHTRIG")
                         .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           seitz 03/15 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoPyramid5fbar::SoPyramid5fbar(int id, int owner)
    : Discret::ELEMENTS::SoPyramid5(id, owner)
{
  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }

  if (Prestress::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(NUMNOD_SOP5, NUMGPT_SOP5 + 1));
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      seitz 03/15 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoPyramid5fbar::SoPyramid5fbar(const Discret::ELEMENTS::SoPyramid5fbar& old)
    : Discret::ELEMENTS::SoPyramid5(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 03/15 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoPyramid5fbar::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoPyramid5fbar(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 03/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoPyramid5fbar::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class So_pyramid5 Element
  Discret::ELEMENTS::SoPyramid5::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 03/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoPyramid5fbar::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_pyramid5 Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Discret::ELEMENTS::SoPyramid5::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              seitz 03/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoPyramid5fbar::Print(std::ostream& os) const
{
  os << "So_pyramid5fbar ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
