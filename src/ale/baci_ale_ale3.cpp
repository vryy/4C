/*----------------------------------------------------------------------------*/
/*! \file

\brief 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_ale3.H"

#include "baci_lib_discret.H"
#include "baci_lib_linedefinition.H"
#include "baci_lib_utils_factory.H"
#include "baci_so3_nullspace.H"
#include "baci_utils_exceptions.H"

DRT::ELEMENTS::Ale3Type DRT::ELEMENTS::Ale3Type::instance_;

DRT::ELEMENTS::Ale3Type& DRT::ELEMENTS::Ale3Type::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Ale3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Ale3* object = new DRT::ELEMENTS::Ale3(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Ale3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele;

  if (eletype == "ALE3")
  {
    if (eledistype != "NURBS27")
    {
      ele = Teuchos::rcp(new DRT::ELEMENTS::Ale3(id, owner));
    }
  }

  return ele;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Ale3Type::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Ale3(id, owner));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Ale3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["ALE3"];

  defs["HEX8"] =
      INPUT::LineDefinition::Builder().AddIntVector("HEX8", 8).AddNamedInt("MAT").Build();

  defs["HEX20"] =
      INPUT::LineDefinition::Builder().AddIntVector("HEX20", 20).AddNamedInt("MAT").Build();

  defs["HEX27"] =
      INPUT::LineDefinition::Builder().AddIntVector("HEX27", 27).AddNamedInt("MAT").Build();

  defs["TET4"] =
      INPUT::LineDefinition::Builder().AddIntVector("TET4", 4).AddNamedInt("MAT").Build();

  defs["TET10"] =
      INPUT::LineDefinition::Builder().AddIntVector("TET10", 10).AddNamedInt("MAT").Build();

  defs["WEDGE6"] =
      INPUT::LineDefinition::Builder().AddIntVector("WEDGE6", 6).AddNamedInt("MAT").Build();

  defs["WEDGE15"] =
      INPUT::LineDefinition::Builder().AddIntVector("WEDGE15", 15).AddNamedInt("MAT").Build();

  defs["PYRAMID5"] =
      INPUT::LineDefinition::Builder().AddIntVector("PYRAMID5", 5).AddNamedInt("MAT").Build();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Ale3SurfaceType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Ale3Surface( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3::Ale3(int id, int owner) : DRT::Element(id, owner) {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3::Ale3(const DRT::ELEMENTS::Ale3& old) : DRT::Element(old) { return; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Ale3::Clone() const
{
  DRT::ELEMENTS::Ale3* newelement = new DRT::ELEMENTS::Ale3(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Ale3::Shape() const
{
  switch (NumNode())
  {
    case 4:
      return tet4;
    case 5:
      return pyramid5;
    case 6:
      return wedge6;
    case 8:
      return hex8;
    case 10:
      return tet10;
    case 20:
      return hex20;
    case 27:
      return hex27;
    default:
      dserror("unexpected number of nodes %d", NumNode());
      break;
  }
  return dis_none;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // Gaussrule
  // AddtoPack(data,gaussrule_);
  // data_
  // vector<char> tmp(0);
  // data_.Pack(tmp);
  // AddtoPack(data,tmp);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // Gaussrule
  // ExtractfromPack(position,data,gaussrule_);
  // data_
  // vector<char> tmp(0);
  // ExtractfromPack(position,data,tmp);
  // data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3::Print(std::ostream& os) const
{
  os << "Ale3 ";
  Element::Print(os);
  std::cout << std::endl;
  // cout << data_;
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Ale3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Ale3Surface, Ale3>(DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Ale3::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}
