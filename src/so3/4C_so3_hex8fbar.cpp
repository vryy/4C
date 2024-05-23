/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element with F-bar modification

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_so3_hex8fbar.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SoHex8fbarType DRT::ELEMENTS::SoHex8fbarType::instance_;

DRT::ELEMENTS::SoHex8fbarType& DRT::ELEMENTS::SoHex8fbarType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8fbarType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoHex8fbar(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8fbarType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoHex8fbar(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8fbarType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoHex8fbar(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoHex8fbarType::nodal_block_information(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoHex8fbarType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::SoHex8fbarType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("GROWTHTRIG")
                     .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8fbar::SoHex8fbar(int id, int owner) : DRT::ELEMENTS::SoHex8(id, owner)
{
  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8, NUMGPT_SOH8 + 1));

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
 |  copy-ctor (public)                                        popp 07/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8fbar::SoHex8fbar(const DRT::ELEMENTS::SoHex8fbar& old)
    : DRT::ELEMENTS::SoHex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoHex8fbar::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoHex8fbar(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8fbar::Pack(CORE::COMM::PackBuffer& data) const
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
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8fbar::Unpack(const std::vector<char>& data)
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
 |  print this element (public)                               popp 07/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8fbar::Print(std::ostream& os) const
{
  os << "So_hex8fbar ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
