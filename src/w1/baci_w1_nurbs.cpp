/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D solid-wall elements using NURBS shape functions.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "baci_w1_nurbs.H"

#include "baci_comm_utils_factory.H"
#include "baci_io_linedefinition.H"
#include "baci_so3_nullspace.H"

DRT::ELEMENTS::NURBS::Wall1NurbsType DRT::ELEMENTS::NURBS::Wall1NurbsType::instance_;

DRT::ELEMENTS::NURBS::Wall1NurbsType& DRT::ELEMENTS::NURBS::Wall1NurbsType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::NURBS::Wall1NurbsType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::NURBS::Wall1Nurbs* object = new DRT::ELEMENTS::NURBS::Wall1Nurbs(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Wall1NurbsType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLNURBS")
  {
    if (eledistype == "NURBS4" || eledistype == "NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Wall1Nurbs(id, owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Wall1NurbsType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Wall1Nurbs(id, owner));
}

void DRT::ELEMENTS::NURBS::Wall1NurbsType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::NURBS::Wall1NurbsType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void DRT::ELEMENTS::NURBS::Wall1NurbsType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALLNURBS"];

  defs["NURBS4"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS4", 4)
                       .AddNamedInt("MAT")
                       .AddNamedString("KINEM")
                       .AddNamedString("EAS")
                       .AddNamedDouble("THICK")
                       .AddNamedString("STRESS_STRAIN")
                       .AddNamedIntVector("GP", 2)
                       .Build();

  defs["NURBS9"] = INPUT::LineDefinition::Builder()
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
DRT::ELEMENTS::NURBS::Wall1Nurbs::Wall1Nurbs(int id, int owner)
    : DRT::ELEMENTS::Wall1::Wall1(id, owner)
{
  SetNurbsElement() = true;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Wall1Nurbs::Wall1Nurbs(const DRT::ELEMENTS::NURBS::Wall1Nurbs& old)
    : DRT::ELEMENTS::Wall1::Wall1(old)
{
  SetNurbsElement() = true;
  return;
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                          gammi 05/09|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::NURBS::Wall1Nurbs::Clone() const
{
  DRT::ELEMENTS::NURBS::Wall1Nurbs* newelement = new DRT::ELEMENTS::NURBS::Wall1Nurbs(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::Wall1Nurbs::Print(std::ostream& os) const
{
  os << "Wall1Nurbs ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::NURBS::Wall1Nurbs::Shape() const
{
  switch (NumNode())
  {
    case 4:
      return CORE::FE::CellType::nurbs4;
    case 9:
      return CORE::FE::CellType::nurbs9;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 05/09|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NURBS::Wall1Nurbs::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Wall1Line, Wall1>(CORE::COMM::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          gammi 05/09|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NURBS::Wall1Nurbs::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}
