/*----------------------------------------------------------------------*/
/*! \file

\brief solid shell8 element formulation

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_sh8.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::SoSh8Type DRT::ELEMENTS::SoSh8Type::instance_;

DRT::ELEMENTS::SoSh8Type& DRT::ELEMENTS::SoSh8Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoSh8Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoSh8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoSh8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoSh8(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoSh8Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoSh8(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoSh8Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoSh8Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::SoSh8Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("EAS")
                     .AddNamedString("ANS")
                     .AddNamedString("THICKDIR")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("STRENGTH")
                     .add_optional_named_double("GROWTHTRIG")
                     .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoSh8::SoSh8(int id, int owner)
    : DRT::ELEMENTS::SoHex8(id, owner),
      thickdir_(undefined),
      anstype_(ansnone),
      nodes_rearranged_(false),
      thickvec_(3, 0.0)
{
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
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoSh8::SoSh8(const DRT::ELEMENTS::SoSh8& old)
    : DRT::ELEMENTS::SoHex8(old),
      thickdir_(old.thickdir_),
      anstype_(old.anstype_),
      nodes_rearranged_(old.nodes_rearranged_),
      thickvec_(old.thickvec_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::SoSh8::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoSh8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_hex8 Element
  DRT::ELEMENTS::SoHex8::Pack(data);
  // thickdir
  AddtoPack(data, thickdir_);
  AddtoPack(data, thickvec_);
  AddtoPack(data, anstype_);
  AddtoPack(data, nodes_rearranged_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::SoHex8::Unpack(basedata);
  // thickdir
  thickdir_ = static_cast<ThicknessDirection>(ExtractInt(position, data));
  ExtractfromPack(position, data, thickvec_);
  anstype_ = static_cast<ANSType>(ExtractInt(position, data));
  nodes_rearranged_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::Print(std::ostream& os) const
{
  os << "So_sh8 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
