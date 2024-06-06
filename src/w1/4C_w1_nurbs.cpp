/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D solid-wall elements using NURBS shape functions.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_nurbs.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Nurbs::Wall1NurbsType Discret::ELEMENTS::Nurbs::Wall1NurbsType::instance_;

Discret::ELEMENTS::Nurbs::Wall1NurbsType& Discret::ELEMENTS::Nurbs::Wall1NurbsType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::Nurbs::Wall1NurbsType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Nurbs::Wall1Nurbs* object = new Discret::ELEMENTS::Nurbs::Wall1Nurbs(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Wall1NurbsType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLNURBS")
  {
    if (eledistype == "NURBS4" || eledistype == "NURBS9")
    {
      return Teuchos::rcp(new Discret::ELEMENTS::Nurbs::Wall1Nurbs(id, owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Wall1NurbsType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Nurbs::Wall1Nurbs(id, owner));
}

void Discret::ELEMENTS::Nurbs::Wall1NurbsType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Nurbs::Wall1NurbsType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void Discret::ELEMENTS::Nurbs::Wall1NurbsType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLNURBS"];

  defs["NURBS4"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS4", 4)
                       .AddNamedInt("MAT")
                       .AddNamedString("KINEM")
                       .AddNamedString("EAS")
                       .AddNamedDouble("THICK")
                       .AddNamedString("STRESS_STRAIN")
                       .AddNamedIntVector("GP", 2)
                       .Build();

  defs["NURBS9"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS9", 9)
                       .AddNamedInt("MAT")
                       .AddNamedString("KINEM")
                       .AddNamedString("EAS")
                       .AddNamedDouble("THICK")
                       .AddNamedString("STRESS_STRAIN")
                       .AddNamedIntVector("GP", 2)
                       .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Wall1Nurbs::Wall1Nurbs(int id, int owner)
    : Discret::ELEMENTS::Wall1::Wall1(id, owner)
{
  SetNurbsElement() = true;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Wall1Nurbs::Wall1Nurbs(const Discret::ELEMENTS::Nurbs::Wall1Nurbs& old)
    : Discret::ELEMENTS::Wall1::Wall1(old)
{
  SetNurbsElement() = true;
  return;
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                          gammi 05/09|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Nurbs::Wall1Nurbs::Clone() const
{
  Discret::ELEMENTS::Nurbs::Wall1Nurbs* newelement =
      new Discret::ELEMENTS::Nurbs::Wall1Nurbs(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/09|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::Wall1Nurbs::Print(std::ostream& os) const
{
  os << "Wall1Nurbs ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Nurbs::Wall1Nurbs::Shape() const
{
  switch (num_node())
  {
    case 4:
      return Core::FE::CellType::nurbs4;
    case 9:
      return Core::FE::CellType::nurbs9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 05/09|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::Wall1Nurbs::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Wall1Line, Wall1>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          gammi 05/09|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::Wall1Nurbs::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
