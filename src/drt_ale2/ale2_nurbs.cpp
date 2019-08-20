/*----------------------------------------------------------------------------*/
/*! \file

\brief Nurbs verison of 2D ALE element

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale2_nurbs.H"

DRT::ELEMENTS::NURBS::Ale2_NurbsType DRT::ELEMENTS::NURBS::Ale2_NurbsType::instance_;

DRT::ELEMENTS::NURBS::Ale2_NurbsType& DRT::ELEMENTS::NURBS::Ale2_NurbsType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::NURBS::Ale2Nurbs* object = new DRT::ELEMENTS::NURBS::Ale2Nurbs(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ALE2")
  {
    if (eledistype == "NURBS4" || eledistype == "NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Ale2Nurbs(id, owner));
    }
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Ale2Nurbs(id, owner));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::Ale2Nurbs(int id, int owner) : DRT::ELEMENTS::Ale2::Ale2(id, owner)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::Ale2Nurbs(const DRT::ELEMENTS::NURBS::Ale2Nurbs& old)
    : DRT::ELEMENTS::Ale2::Ale2(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::~Ale2Nurbs() { return; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::Ale2Nurbs::Print(std::ostream& os) const
{
  os << "Ale2Nurbs ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Ale2Nurbs::Shape() const
{
  switch (NumNode())
  {
    case 4:
      return nurbs4;
    case 9:
      return nurbs9;
    default:
      dserror("unexpected number of nodes %d", NumNode());
      break;
  }

  return dis_none;
}
