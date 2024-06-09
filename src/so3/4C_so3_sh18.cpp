/*----------------------------------------------------------------------*/
/*! \file

\brief ToDo Add meaningful comment.

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_sh18.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::SoSh18Type Discret::ELEMENTS::SoSh18Type::instance_;

Discret::ELEMENTS::SoSh18Type& Discret::ELEMENTS::SoSh18Type::Instance() { return instance_; }
namespace
{
  const std::string name = Discret::ELEMENTS::SoSh18Type::Instance().Name();
}

Core::Communication::ParObject* Discret::ELEMENTS::SoSh18Type::Create(const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoSh18(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoSh18Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoSh18(id, owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoSh18Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoSh18(id, owner));
  return ele;
}

void Discret::ELEMENTS::SoSh18Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX18"] = Input::LineDefinition::Builder()
                      .AddIntVector("HEX18", 18)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedString("TSL")
                      .AddNamedString("MEL")
                      .AddNamedString("CTL")
                      .AddNamedString("VOL")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           seitz 11/14 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoSh18::SoSh18(int id, int owner) : SoBase(id, owner), SoHex18(id, owner)
{
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
 |  copy-ctor (public)                                      seitz 11/14 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoSh18::SoSh18(const Discret::ELEMENTS::SoSh18& old)
    : SoBase(old),
      SoHex18(old),
      dsg_shear_(old.dsg_shear_),
      dsg_membrane_(old.dsg_membrane_),
      dsg_ctl_(old.dsg_ctl_),
      eas_(old.eas_)
{
  setup_dsg();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoSh18::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoSh18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  SoBase::Pack(data);

  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);

  // element technology bools
  add_to_pack(data, (int)dsg_shear_);
  add_to_pack(data, (int)dsg_membrane_);
  add_to_pack(data, (int)dsg_ctl_);
  add_to_pack(data, (int)eas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  SoBase::Unpack(basedata);

  // detJ_
  extract_from_pack(position, data, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(position, data, invJ_[i]);

  // element technology bools
  dsg_shear_ = ExtractInt(position, data);
  dsg_membrane_ = ExtractInt(position, data);
  dsg_ctl_ = ExtractInt(position, data);
  eas_ = ExtractInt(position, data);
  setup_dsg();

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::Print(std::ostream& os) const
{
  os << "So_sh18 ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
