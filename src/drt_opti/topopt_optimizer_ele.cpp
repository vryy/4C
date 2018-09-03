/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele.cpp

\brief element routines of the topology optimizer

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer_ele.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_mat/material.H"
#include "../drt_lib/drt_linedefinition.H"


DRT::ELEMENTS::TopOptType DRT::ELEMENTS::TopOptType::instance_;

DRT::ELEMENTS::TopOptType& DRT::ELEMENTS::TopOptType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::TopOptType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::TopOpt* object = new DRT::ELEMENTS::TopOpt(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TopOptType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TOPOPT")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::TopOpt(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TopOptType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::TopOpt(id, owner));
  return ele;
}


void DRT::ELEMENTS::TopOptType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["TOPOPT"];

  defs["HEX8"].AddIntVector("HEX8", 8).AddNamedInt("MAT");

  defs["HEX20"].AddIntVector("HEX20", 20).AddNamedInt("MAT");

  defs["HEX27"].AddIntVector("HEX27", 27).AddNamedInt("MAT");

  defs["TET4"].AddIntVector("TET4", 4).AddNamedInt("MAT");

  defs["TET10"].AddIntVector("TET10", 10).AddNamedInt("MAT");

  defs["QUAD4"].AddIntVector("QUAD4", 4).AddNamedInt("MAT");

  defs["QUAD8"].AddIntVector("QUAD8", 8).AddNamedInt("MAT");

  defs["QUAD9"].AddIntVector("QUAD9", 9).AddNamedInt("MAT");

  defs["TRI3"].AddIntVector("TRI3", 3).AddNamedInt("MAT");

  defs["TRI6"].AddIntVector("TRI6", 6).AddNamedInt("MAT");

  defs["LINE2"].AddIntVector("LINE2", 2).AddNamedInt("MAT");

  defs["LINE3"].AddIntVector("LINE3", 3).AddNamedInt("MAT");
}


DRT::ELEMENTS::TopOptBoundaryType DRT::ELEMENTS::TopOptBoundaryType::instance_;

DRT::ELEMENTS::TopOptBoundaryType& DRT::ELEMENTS::TopOptBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TopOptBoundaryType::Create(const int id, const int owner)
{
  dserror("is this used?");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOpt::TopOpt(int id, int owner) : DRT::Element(id, owner), distype_(dis_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOpt::TopOpt(const DRT::ELEMENTS::TopOpt& old)
    : DRT::Element(old), distype_(old.distype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Transport and return pointer to it (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::TopOpt::Clone() const
{
  DRT::ELEMENTS::TopOpt* newelement = new DRT::ELEMENTS::TopOpt(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOpt::SetMaterial(int matnum)
{
  // the standard part:
  // mat_ = MAT::Material::Factory(matnum);  // not allowed since mat_ is private
  DRT::Element::SetMaterial(matnum);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType() != INPAR::MAT::m_opti_dens)
    dserror("Topology optimization element got unsupported material type %d", mat->MaterialType());

  // TODO if this works, only the basic function is required

  return;
}


/*----------------------------------------------------------------------*
 |  Return the shape of a Transport element                      (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::TopOpt::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOpt::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  // distype
  AddtoPack(data, distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOpt::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // distype
  distype_ = static_cast<DiscretizationType>(ExtractInt(position, data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)           gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TopOpt::NumLine() const { return DRT::UTILS::getNumberOfElementLines(distype_); }


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)        gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TopOpt::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)          gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TopOpt::NumVolume() const
{
  return DRT::UTILS::getNumberOfElementVolumes(distype_);
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOpt::~TopOpt() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOpt::Print(std::ostream& os) const
{
  os << "Topology optimization element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << distype_ << std::endl;
  std::cout << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                  g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TopOpt::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumLine() > 1)  // 3D and 2D
    return DRT::UTILS::ElementBoundaryFactory<TopOptBoundary, TopOpt>(DRT::UTILS::buildLines, this);
  else
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> lines(1);
    lines[0] = Teuchos::rcp(this, false);
    return lines;
  }
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TopOpt::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  if (NumSurface() > 1)  // 3D
    return DRT::UTILS::ElementBoundaryFactory<TopOptBoundary, TopOpt>(
        DRT::UTILS::buildSurfaces, this);
  else if (NumSurface() == 1)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> surfaces(1);
    surfaces[0] = Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    // 1D
    dserror("Surfaces() for 1D-Transport element not implemented");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TopOpt::Volumes()
{
  if (NumVolume() == 1)
  {
    std::vector<Teuchos::RCP<Element>> volumes(1);
    volumes[0] = Teuchos::rcp(this, false);
    return volumes;
  }
  else
  {
    dserror("Volumes() for 1D-/2D-Transport element not implemented");
    return DRT::Element::Volumes();
  }
}


/*----------------------------------------------------------------------*
 | read element input (public)                                          |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::TopOpt::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SetDisType(DRT::StringToDistype(distype));

  return true;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundary::TopOptBoundary(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::TopOpt* parent, const int lbeleid)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lbeleid);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundary::TopOptBoundary(const DRT::ELEMENTS::TopOptBoundary& old)
    : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::TopOptBoundary::Clone() const
{
  DRT::ELEMENTS::TopOptBoundary* newelement = new DRT::ELEMENTS::TopOptBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::TopOptBoundary::Shape() const
{
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptBoundary::Pack(DRT::PackBuffer& data) const
{
  dserror("This TransportBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptBoundary::Unpack(const std::vector<char>& data)
{
  dserror("This TransportBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundary::~TopOptBoundary() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptBoundary::Print(std::ostream& os) const
{
  os << "Topology optimization boundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Shape() << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)        gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TopOptBoundary::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)    gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TopOptBoundary::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TopOptBoundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of TransportBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element>> lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TopOptBoundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of TransportBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element>> surfaces(0);
  return surfaces;
}
