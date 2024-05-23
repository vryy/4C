/*----------------------------------------------------------------------*/
/*! \file
\brief A 3D constraint element with no physics attached
\level 2


*----------------------------------------------------------------------*/

#include "4C_constraint_element3.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::ConstraintElement3Type DRT::ELEMENTS::ConstraintElement3Type::instance_;


DRT::ELEMENTS::ConstraintElement3Type& DRT::ELEMENTS::ConstraintElement3Type::Instance()
{
  return instance_;
}


CORE::COMM::ParObject* DRT::ELEMENTS::ConstraintElement3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::ConstraintElement3* object = new DRT::ELEMENTS::ConstraintElement3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ConstraintElement3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "CONSTRELE3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::ConstraintElement3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ConstraintElement3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::ConstraintElement3(id, owner));
  return ele;
}


void DRT::ELEMENTS::ConstraintElement3Type::nodal_block_information(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::ConstraintElement3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(int id, int owner) : DRT::Element(id, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(const DRT::ELEMENTS::ConstraintElement3& old)
    : DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ConstraintElement3::Clone() const
{
  DRT::ELEMENTS::ConstraintElement3* newelement = new DRT::ELEMENTS::ConstraintElement3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3::Pack(CORE::COMM::PackBuffer& data) const
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
void DRT::ELEMENTS::ConstraintElement3::Unpack(const std::vector<char>& data)
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
