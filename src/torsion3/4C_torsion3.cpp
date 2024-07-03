/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional torsion spring element

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_torsion3.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Torsion3Type Discret::ELEMENTS::Torsion3Type::instance_;

Discret::ELEMENTS::Torsion3Type& Discret::ELEMENTS::Torsion3Type::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::Torsion3Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Torsion3* object = new Discret::ELEMENTS::Torsion3(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Torsion3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TORSION3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Torsion3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Torsion3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Torsion3(id, owner));
  return ele;
}


void Discret::ELEMENTS::Torsion3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Torsion3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::Torsion3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["TORSION3"];

  defs["LINE3"] = Input::LineDefinition::Builder()
                      .add_int_vector("LINE3", 3)
                      .add_named_int("MAT")
                      .add_named_string("BENDINGPOTENTIAL")
                      .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Torsion3::Torsion3(int id, int owner) : Core::Elements::Element(id, owner)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Torsion3::Torsion3(const Discret::ELEMENTS::Torsion3& old)
    : Core::Elements::Element(old)
{
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Torsion3 and return pointer to it (public)|
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Torsion3::Clone() const
{
  Discret::ELEMENTS::Torsion3* newelement = new Discret::ELEMENTS::Torsion3(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Torsion3::print(std::ostream& os) const { return; }


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Torsion3::Shape() const { return Core::FE::CellType::line3; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Torsion3::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  Element::pack(data);
  add_to_pack(data, bendingpotential_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Torsion3::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);
  bendingpotential_ = static_cast<BendingPotential>(extract_int(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Torsion3::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Torsion3::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<Solid::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::ParamsInterface> Discret::ELEMENTS::Torsion3::ParamsInterfacePtr()
{
  return interface_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
