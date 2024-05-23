/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_so3_hex8p1j1.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SoHex8P1J1Type DRT::ELEMENTS::SoHex8P1J1Type::instance_;

DRT::ELEMENTS::SoHex8P1J1Type& DRT::ELEMENTS::SoHex8P1J1Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8P1J1Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoHex8P1J1(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8P1J1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoHex8P1J1(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8P1J1Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoHex8P1J1(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoHex8P1J1Type::nodal_block_information(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  //   numdf = 3;
  //   dimns = 6;
  //   nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoHex8P1J1Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}

void DRT::ELEMENTS::SoHex8P1J1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .Build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                               lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1::SoHex8P1J1(int id, int owner) : DRT::ELEMENTS::SoHex8(id, owner)
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
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1::SoHex8P1J1(const DRT::ELEMENTS::SoHex8P1J1& old)
    : DRT::ELEMENTS::SoHex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoHex8P1J1::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoHex8P1J1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_hex8 Element
  DRT::ELEMENTS::SoHex8::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::SoHex8::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                 lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1::Print(std::ostream& os) const
{
  os << "So_Hex8P1J1 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
