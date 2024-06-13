/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_so3_hex8p1j1.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::SoHex8P1J1Type Discret::ELEMENTS::SoHex8P1J1Type::instance_;

Discret::ELEMENTS::SoHex8P1J1Type& Discret::ELEMENTS::SoHex8P1J1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8P1J1Type::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoHex8P1J1(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8P1J1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoHex8P1J1(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8P1J1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoHex8P1J1(id, owner));
  return ele;
}


void Discret::ELEMENTS::SoHex8P1J1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  //   numdf = 3;
  //   dimns = 6;
  //   nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoHex8P1J1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}

void Discret::ELEMENTS::SoHex8P1J1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder()
                     .add_int_vector("HEX8", 8)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                               lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8P1J1::SoHex8P1J1(int id, int owner) : Discret::ELEMENTS::SoHex8(id, owner)
{
  k_pu_.PutScalar(0.0);
  k_tu_.PutScalar(0.0);

  r_t_.PutScalar(0.0);
  r_p_.PutScalar(0.0);

  k_tt_ = 0.0;
  k_pt_ = 0.0;

  p_.PutScalar(0.0);
  p_o_.PutScalar(0.0);

  t_.PutScalar(1.0);
  t_o_.PutScalar(1.0);

  m_.PutScalar(0.0);
  for (int i = 0; i < 3; ++i)
  {
    m_(i, 0) = 1.0;
  }

  identity6_.PutScalar(0.0);
  for (int i = 0; i < 6; ++i)
  {
    identity6_(i, i) = 1.0;
  }

  i_d_ = identity6_;
  i_d_.MultiplyNT(-1.0 / 3.0, m_, m_, 1.0);

  i_0_.PutScalar(0.0);

  for (int i = 0; i < 3; ++i)
  {
    i_0_(i, i) = 1.0;
  }
  for (int i = 3; i < 6; ++i)
  {
    i_0_(i, i) = 0.5;
  }

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8P1J1::SoHex8P1J1(const Discret::ELEMENTS::SoHex8P1J1& old)
    : Discret::ELEMENTS::SoHex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoHex8P1J1::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoHex8P1J1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8P1J1::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class So_hex8 Element
  Discret::ELEMENTS::SoHex8::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8P1J1::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Discret::ELEMENTS::SoHex8::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                 lw 12/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8P1J1::Print(std::ostream& os) const
{
  os << "So_Hex8P1J1 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
