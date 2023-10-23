/*----------------------------------------------------------------------------*/
/*! \file
\brief wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "baci_w1.H"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_lib_discret.H"
#include "baci_lib_linedefinition.H"
#include "baci_lib_utils_factory.H"
#include "baci_so3_nullspace.H"
#include "baci_utils_exceptions.H"

DRT::ELEMENTS::Wall1Type DRT::ELEMENTS::Wall1Type::instance_;

DRT::ELEMENTS::Wall1Type& DRT::ELEMENTS::Wall1Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Wall1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1* object = new DRT::ELEMENTS::Wall1(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALL")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::Wall1(id, owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1Type::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Wall1(id, owner));
}


void DRT::ELEMENTS::Wall1Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Wall1Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, int const numdof, int const dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
  ;
}

void DRT::ELEMENTS::Wall1Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALL"];

  defs["QUAD4"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD4", 4)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedString("EAS")
                      .AddNamedDouble("THICK")
                      .AddNamedString("STRESS_STRAIN")
                      .AddNamedIntVector("GP", 2)
                      .Build();

  defs["QUAD8"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD8", 8)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedString("EAS")
                      .AddNamedDouble("THICK")
                      .AddNamedString("STRESS_STRAIN")
                      .AddNamedIntVector("GP", 2)
                      .Build();

  defs["QUAD9"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD9", 9)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedString("EAS")
                      .AddNamedDouble("THICK")
                      .AddNamedString("STRESS_STRAIN")
                      .AddNamedIntVector("GP", 2)
                      .Build();

  defs["TRI3"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI3", 3)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("EAS")
                     .AddNamedDouble("THICK")
                     .AddNamedString("STRESS_STRAIN")
                     .AddNamedIntVector("GP", 2)
                     .Build();

  defs["TRI6"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI6", 6)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("EAS")
                     .AddNamedDouble("THICK")
                     .AddNamedString("STRESS_STRAIN")
                     .AddNamedIntVector("GP", 2)
                     .Build();

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


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1LineType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Wall1Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/08/|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(int id, int owner)
    : So_base(id, owner),
      data_(),
      material_(0),
      thickness_(0.0),
      old_step_length_(0.0),
      gaussrule_(CORE::DRT::UTILS::GaussRule2D::undefined),
      wtype_(plane_none),
      stresstype_(w1_none),
      iseas_(false),
      eastype_(eas_vague),
      structale_(false),
      distype_(DRT::Element::DiscretizationType::dis_none)
{
  if (DRT::Problem::Instance()->GetProblemType() == ProblemType::struct_ale) structale_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(const DRT::ELEMENTS::Wall1& old)
    : So_base(old),
      data_(old.data_),
      material_(old.material_),
      thickness_(old.thickness_),
      old_step_length_(old.old_step_length_),
      gaussrule_(old.gaussrule_),
      wtype_(old.wtype_),
      stresstype_(old.stresstype_),
      iseas_(old.iseas_),
      eastype_(old.eas_vague),
      structale_(old.structale_),
      distype_(old.distype_)
// tsi_couptyp_(old.tsi_couptyp_)

{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1::Clone() const
{
  DRT::ELEMENTS::Wall1* newelement = new DRT::ELEMENTS::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Wall1::Shape() const { return distype_; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);
  // material_
  AddtoPack(data, material_);
  // thickness
  AddtoPack(data, thickness_);
  // plane strain or plane stress information
  AddtoPack(data, wtype_);
  // gaussrule_
  AddtoPack(data, gaussrule_);
  // stresstype
  AddtoPack(data, stresstype_);
  // eas
  AddtoPack(data, iseas_);
  // eas type
  AddtoPack(data, eastype_);
  // structale
  AddtoPack(data, structale_);
  // distype
  AddtoPack(data, distype_);
  // data
  AddtoPack(data, data_);
  // line search
  AddtoPack(data, old_step_length_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  So_base::Unpack(basedata);
  // material_
  ExtractfromPack(position, data, material_);
  // thickness_
  ExtractfromPack(position, data, thickness_);
  // plane strain or plane stress information_
  wtype_ = static_cast<DimensionalReduction>(ExtractInt(position, data));
  // gaussrule_
  ExtractfromPack(position, data, gaussrule_);
  // stresstype_
  stresstype_ = static_cast<StressType>(ExtractInt(position, data));
  // iseas_
  iseas_ = ExtractInt(position, data);
  // eastype_
  eastype_ = static_cast<EasType>(ExtractInt(position, data));
  // structale_
  structale_ = ExtractInt(position, data);
  // distype_
  distype_ = static_cast<DiscretizationType>(ExtractInt(position, data));
  // data
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);
  // line search
  ExtractfromPack(position, data, old_step_length_);
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::~Wall1() { return; }


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 07/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line, Wall1>(DRT::UTILS::buildLines, this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1::Surfaces()
{
  std::vector<Teuchos::RCP<Element>> surfaces(1);
  surfaces[0] = Teuchos::rcp(this, false);
  return surfaces;
}

/*-----------------------------------------------------------------------------*
| Map plane Green-Lagrange strains to 3d                       mayr.mt 05/2014 |
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::GreenLagrangePlane3d(
    const CORE::LINALG::SerialDenseVector& glplane, CORE::LINALG::Matrix<6, 1>& gl3d)
{
  gl3d(0) = glplane(0);               // E_{11}
  gl3d(1) = glplane(1);               // E_{22}
  gl3d(2) = 0.0;                      // E_{33}
  gl3d(3) = glplane(2) + glplane(3);  // 2*E_{12}=E_{12}+E_{21}
  gl3d(4) = 0.0;                      // 2*E_{23}
  gl3d(5) = 0.0;                      // 2*E_{31}

  return;
}
