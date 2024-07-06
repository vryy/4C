/*----------------------------------------------------------------------*/
/*! \file
\brief A 2D constraint element with no physics attached
\level 2


*----------------------------------------------------------------------*/

#include "4C_constraint_element2.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::ConstraintElement2Type Discret::ELEMENTS::ConstraintElement2Type::instance_;


Discret::ELEMENTS::ConstraintElement2Type& Discret::ELEMENTS::ConstraintElement2Type::instance()
{
  return instance_;
}


Core::Communication::ParObject* Discret::ELEMENTS::ConstraintElement2Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::ConstraintElement2* object = new Discret::ELEMENTS::ConstraintElement2(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ConstraintElement2Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "CONSTRELE2")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::ConstraintElement2(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ConstraintElement2Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::ConstraintElement2(id, owner));
  return ele;
}


void Discret::ELEMENTS::ConstraintElement2Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::ConstraintElement2Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ConstraintElement2::ConstraintElement2(int id, int owner)
    : Core::Elements::Element(id, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ConstraintElement2::ConstraintElement2(
    const Discret::ELEMENTS::ConstraintElement2& old)
    : Core::Elements::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ConstraintElement2::clone() const
{
  Discret::ELEMENTS::ConstraintElement2* newelement =
      new Discret::ELEMENTS::ConstraintElement2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ConstraintElement2::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ConstraintElement2::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

FOUR_C_NAMESPACE_CLOSE
