/*----------------------------------------------------------------------------*/
/*! \file

\brief Nurbs verison of 2D ALE element

\level 3

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale2_nurbs.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Nurbs::Ale2NurbsType Discret::ELEMENTS::Nurbs::Ale2NurbsType::instance_;

Discret::ELEMENTS::Nurbs::Ale2NurbsType& Discret::ELEMENTS::Nurbs::Ale2NurbsType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::Nurbs::Ale2NurbsType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Nurbs::Ale2Nurbs* object = new Discret::ELEMENTS::Nurbs::Ale2Nurbs(-1, -1);
  object->unpack(data);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Ale2NurbsType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ALE2")
  {
    if (eledistype == "NURBS4" || eledistype == "NURBS9")
    {
      return Teuchos::rcp(new Discret::ELEMENTS::Nurbs::Ale2Nurbs(id, owner));
    }
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Ale2NurbsType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Nurbs::Ale2Nurbs(id, owner));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Ale2Nurbs::Ale2Nurbs(int id, int owner)
    : Discret::ELEMENTS::Ale2::Ale2(id, owner)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Ale2Nurbs::Ale2Nurbs(const Discret::ELEMENTS::Nurbs::Ale2Nurbs& old)
    : Discret::ELEMENTS::Ale2::Ale2(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::Ale2Nurbs::print(std::ostream& os) const
{
  os << "Ale2Nurbs ";
  Element::print(os);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Nurbs::Ale2Nurbs::Shape() const
{
  switch (num_node())
  {
    case 4:
      return Core::FE::CellType::nurbs4;
    case 9:
      return Core::FE::CellType::nurbs9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
