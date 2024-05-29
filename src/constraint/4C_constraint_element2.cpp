/*----------------------------------------------------------------------*/
/*! \file
\brief A 2D constraint element with no physics attached
\level 2


*----------------------------------------------------------------------*/

#include "4C_constraint_element2.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::ConstraintElement2Type DRT::ELEMENTS::ConstraintElement2Type::instance_;


DRT::ELEMENTS::ConstraintElement2Type& DRT::ELEMENTS::ConstraintElement2Type::Instance()
{
  return instance_;
}


CORE::COMM::ParObject* DRT::ELEMENTS::ConstraintElement2Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::ConstraintElement2* object = new DRT::ELEMENTS::ConstraintElement2(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::ConstraintElement2Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "CONSTRELE2")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::ConstraintElement2(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::ConstraintElement2Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::ConstraintElement2(id, owner));
  return ele;
}


void DRT::ELEMENTS::ConstraintElement2Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::ConstraintElement2Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement2::ConstraintElement2(int id, int owner)
    : CORE::Elements::Element(id, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement2::ConstraintElement2(const DRT::ELEMENTS::ConstraintElement2& old)
    : CORE::Elements::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::ConstraintElement2::Clone() const
{
  DRT::ELEMENTS::ConstraintElement2* newelement = new DRT::ELEMENTS::ConstraintElement2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement2::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement2::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

FOUR_C_NAMESPACE_CLOSE
