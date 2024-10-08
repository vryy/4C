/*----------------------------------------------------------------------*/
/*! \file
\brief A 3D constraint element with no physics attached
\level 2


*----------------------------------------------------------------------*/

#include "4C_constraint_element3.hpp"

#include "4C_comm_pack_helpers.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::ConstraintElement3Type Discret::ELEMENTS::ConstraintElement3Type::instance_;


Discret::ELEMENTS::ConstraintElement3Type& Discret::ELEMENTS::ConstraintElement3Type::instance()
{
  return instance_;
}


Core::Communication::ParObject* Discret::ELEMENTS::ConstraintElement3Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::ELEMENTS::ConstraintElement3* object = new Discret::ELEMENTS::ConstraintElement3(-1, -1);
  object->unpack(buffer);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ConstraintElement3Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "CONSTRELE3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::RCP(new Discret::ELEMENTS::ConstraintElement3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ConstraintElement3Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::RCP(new Discret::ELEMENTS::ConstraintElement3(id, owner));
  return ele;
}


void Discret::ELEMENTS::ConstraintElement3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::ConstraintElement3Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ConstraintElement3::ConstraintElement3(int id, int owner)
    : Core::Elements::Element(id, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ConstraintElement3::ConstraintElement3(
    const Discret::ELEMENTS::ConstraintElement3& old)
    : Core::Elements::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ConstraintElement3::clone() const
{
  Discret::ELEMENTS::ConstraintElement3* newelement =
      new Discret::ELEMENTS::ConstraintElement3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ConstraintElement3::pack(Core::Communication::PackBuffer& data) const
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
void Discret::ELEMENTS::ConstraintElement3::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Element::unpack(base_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}

FOUR_C_NAMESPACE_CLOSE
