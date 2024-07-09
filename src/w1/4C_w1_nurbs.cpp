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

Discret::ELEMENTS::Nurbs::Wall1NurbsType& Discret::ELEMENTS::Nurbs::Wall1NurbsType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::Nurbs::Wall1NurbsType::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Nurbs::Wall1Nurbs* object = new Discret::ELEMENTS::Nurbs::Wall1Nurbs(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Wall1NurbsType::create(
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


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::Wall1NurbsType::create(
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

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Nurbs::Wall1NurbsType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void Discret::ELEMENTS::Nurbs::Wall1NurbsType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLNURBS"];

  defs["NURBS4"] = Input::LineDefinition::Builder()
                       .add_int_vector("NURBS4", 4)
                       .add_named_int("MAT")
                       .add_named_string("KINEM")
                       .add_named_string("EAS")
                       .add_named_double("THICK")
                       .add_named_string("STRESS_STRAIN")
                       .add_named_int_vector("GP", 2)
                       .build();

  defs["NURBS9"] = Input::LineDefinition::Builder()
                       .add_int_vector("NURBS9", 9)
                       .add_named_int("MAT")
                       .add_named_string("KINEM")
                       .add_named_string("EAS")
                       .add_named_double("THICK")
                       .add_named_string("STRESS_STRAIN")
                       .add_named_int_vector("GP", 2)
                       .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Wall1Nurbs::Wall1Nurbs(int id, int owner)
    : Discret::ELEMENTS::Wall1::Wall1(id, owner)
{
  set_nurbs_element() = true;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::Wall1Nurbs::Wall1Nurbs(const Discret::ELEMENTS::Nurbs::Wall1Nurbs& old)
    : Discret::ELEMENTS::Wall1::Wall1(old)
{
  set_nurbs_element() = true;
  return;
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                          gammi 05/09|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Nurbs::Wall1Nurbs::clone() const
{
  Discret::ELEMENTS::Nurbs::Wall1Nurbs* newelement =
      new Discret::ELEMENTS::Nurbs::Wall1Nurbs(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/09|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::Wall1Nurbs::print(std::ostream& os) const
{
  os << "Wall1Nurbs ";
  Element::print(os);
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Nurbs::Wall1Nurbs::shape() const
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
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::Wall1Nurbs::lines()
{
  return Core::Communication::ElementBoundaryFactory<Wall1Line, Wall1>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          gammi 05/09|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::Wall1Nurbs::surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
