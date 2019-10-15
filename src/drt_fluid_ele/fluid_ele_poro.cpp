/*-----------------------------------------------------------*/
/*! \file

\brief Fluid element for poroelasticity problems

\maintainer Johannes Kremheller

\level 2

*/
/*-----------------------------------------------------------*/


#include "fluid_ele_poro.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_inpar/inpar_fluid.H"

DRT::ELEMENTS::FluidPoroEleType DRT::ELEMENTS::FluidPoroEleType::instance_;

/*----------------------------------------------------------------------*
 *   return instance (public,static)                                    |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoroEleType& DRT::ELEMENTS::FluidPoroEleType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 *   create element (public)                                 vuong 06/13|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::FluidPoroEleType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::FluidPoro* object = new DRT::ELEMENTS::FluidPoro(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *   create element (public)                                 vuong 06/13|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroEleType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDPORO")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::FluidPoro(id, owner));
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *   create element (public)                                 vuong 06/13|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroEleType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::FluidPoro(id, owner));
}

/*----------------------------------------------------------------------*
 *   Setup Element Definition                               vuong 06/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoroEleType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_fluid;
  FluidType::SetupElementDefinition(definitions_fluid);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["FLUIDPORO"];

  // 3D
  defs["HEX8"] = defs_fluid["HEX8"];
  defs["HEX20"] = defs_fluid["HEX20"];
  defs["HEX27"] = defs_fluid["HEX27"];
  defs["TET4"] = defs_fluid["TET4"];
  defs["TET10"] = defs_fluid["TET10"];
  defs["WEDGE6"] = defs_fluid["WEDGE6"];
  defs["WEDGE15"] = defs_fluid["WEDGE15"];
  defs["PYRAMID5"] = defs_fluid["PYRAMID5"];
  defs["NURBS8"] = defs_fluid["NURBS8"];
  defs["NURBS27"] = defs_fluid["NURBS27"];

  // 2D
  defs["QUAD4"] = defs_fluid["QUAD4"];
  defs["QUAD8"] = defs_fluid["QUAD8"];
  defs["QUAD9"] = defs_fluid["QUAD9"];
  defs["TRI3"] = defs_fluid["TRI3"];
  defs["TRI6"] = defs_fluid["TRI6"];
  defs["NURBS4"] = defs_fluid["NURBS4"];
  defs["NURBS9"] = defs_fluid["NURBS9"];
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 06/13 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoro::FluidPoro(int id, int owner)
    : Fluid(id, owner), kintype_(INPAR::STR::kinem_vague)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                     vuong 06/13   |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoro::FluidPoro(const DRT::ELEMENTS::FluidPoro& old)
    : Fluid(old), kintype_(old.kintype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                     vuong 06/13    |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidPoro::Clone() const
{
  DRT::ELEMENTS::FluidPoro* newelement = new DRT::ELEMENTS::FluidPoro(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         vuong 06/13  |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // kinemtics type
  AddtoPack(data, kintype_);

  // add base class Element
  Fluid::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         vuong 06/13   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::ExtractfromPack(position, data, basedata);
  Fluid::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)           vuong 06/13     |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidPoro::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine() > 1)  // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<FluidPoroBoundary, FluidPoro>(
        DRT::UTILS::buildLines, this);
  }
  else if (NumLine() ==
           1)  // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> surfaces(1);
    surfaces[0] = Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                       vuong 06/13    |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidPoro::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)  // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<FluidPoroBoundary, FluidPoro>(
        DRT::UTILS::buildSurfaces, this);
  else if (NumSurface() ==
           1)  // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> surfaces(1);
    surfaces[0] = Teuchos::rcp(this, false);
    return surfaces;
  }
  else  // 1D elements
  {
    dserror("Surfaces() does not exist for 1D-element ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)           vuong 06/13      |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidPoro::Volumes()
{
  if (NumVolume() ==
      1)  // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element>> volumes(1);
    volumes[0] = Teuchos::rcp(this, false);
    return volumes;
  }
  else  //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Surfaces();
  }
}

/*----------------------------------------------------------------------*
 |  print this element (public)                         vuong 06/13    |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Print(std::ostream& os) const
{
  os << "FluidPoro " << (DistypeToString(distype_)).c_str();
  Element::Print(os);
  return;
}
