/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell.cpp

\brief fluid cell for meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

/*==========================================================================*\
 *                                                                          *
 * class MeshfreeFluidType                                                  *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 |  create instance of MeshfreeFluidType                 (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidType DRT::ELEMENTS::MeshfreeFluidType::instance_;


/*--------------------------------------------------------------------------*
 |  create parallel object of MeshfreeFluid cell         (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::MeshfreeFluidType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::MeshfreeFluid* object = new DRT::ELEMENTS::MeshfreeFluid(-1,-1);
  object->Unpack(data);
  return object;
}


/*--------------------------------------------------------------------------*
 |  create object of MeshfreeFluid type                  (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MeshfreeFluidType::Create(
  const std::string eletype,
  const std::string eledistype,
  const int id,
  const int owner)
{
  if (eletype=="MEFLUID")
    return Teuchos::rcp(new DRT::ELEMENTS::MeshfreeFluid(id,owner));

  return Teuchos::null;
}


/*--------------------------------------------------------------------------*
 |  create object of MeshfreeFluid type                  (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MeshfreeFluidType::Create(
  const int id,
  const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::MeshfreeFluid(id,owner));
}


/*--------------------------------------------------------------------------*
 |                                                       (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidType::NodalBlockInformation(
  DRT::Element * dwele,
  int & numdf,
  int & dimns,
  int & nv,
  int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf-1;
  np = 1;
}


/*--------------------------------------------------------------------------*
 |                                                       (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidType::ComputeNullSpace(
  DRT::Discretization & dis,
  std::vector<double> & ns,
  const double * x0,
  int numdf,
  int dimns )
{
  DRT::UTILS::ComputeFluidDNullSpace( dis, ns, x0, numdf, dimns );
}

/*--------------------------------------------------------------------------*
 |                                                       (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidType::SetupElementDefinition(
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions
  )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["MEFLUID"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  // 2D elements
  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;
}


/*==========================================================================*\
 *                                                                          *
 * class MeshfreeFluid                                                      *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluid::MeshfreeFluid(int id, int owner) :
  DRT::Element(id,owner),  // necessary due to virtual inheritance from DRT::Element
  DRT::MESHFREE::Cell(id,owner),
  distype_(DRT::Element::dis_none),
  is_ale_(false)
{
    return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluid::MeshfreeFluid(const DRT::ELEMENTS::MeshfreeFluid& old) :
  DRT::Element(old),  // necessary due to virtual inheritance from DRT::Element
  DRT::MESHFREE::Cell(old),
  distype_    (old.distype_),
  is_ale_     (old.is_ale_)
{
  return;
}

/*--------------------------------------------------------------------------*
 | returns pointer to deep copy of MeshfreeTransport     (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::MeshfreeFluid::Clone() const
{
  DRT::ELEMENTS::MeshfreeFluid* newelement = new DRT::ELEMENTS::MeshfreeFluid(*this);
  return newelement;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluid::~MeshfreeFluid()
{
  return;
}

/*--------------------------------------------------------------------------*
 |  Return the shape of a MeshfreeFluid cell             (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::MeshfreeFluid::Shape() const
{
  return distype_;
}

/*--------------------------------------------------------------------------*
 |  Return number of lines of this element               (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(distype_);
}


/*--------------------------------------------------------------------------*
 |  Return number of surfaces of this element            (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(distype_);
}


/*--------------------------------------------------------------------------*
 | Return number of volumes of this element              (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::NumVolume() const
{
  return DRT::UTILS::getNumberOfElementVolumes(distype_);
}


/*--------------------------------------------------------------------------*
 | Get vector of Teuchos::RCPs to lines of this cell     (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::MeshfreeFluid::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<MeshfreeFluidBoundary,MeshfreeFluid>(DRT::UTILS::buildLines,this);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > lines(1);
    lines[0]= Teuchos::rcp(this, false);
    return lines;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Lines();
  }
}


/*--------------------------------------------------------------------------*
 | Get vector of Teuchos::RCPs to surface of this cell   (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::MeshfreeFluid::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<MeshfreeFluidBoundary,MeshfreeFluid>(DRT::UTILS::buildSurfaces,this);
  }
  else if (NumSurface() == 1) // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else  // 1D elements
  {
    dserror("Surfaces() does not exist for 1D-element ");
    return DRT::Element::Surfaces();
  }
}


/*--------------------------------------------------------------------------*
 | Get vector of Teuchos::RCPs to volumes of this cell   (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::MeshfreeFluid::Volumes()
{
  if (NumVolume()==1) // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element> > volumes(1);
    volumes[0]= Teuchos::rcp(this, false);
    return volumes;
  }
  else //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Volumes();
  }
}

/*--------------------------------------------------------------------------*
 |  Pack data                                            (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  DRT::MESHFREE::Cell::Pack(data);
  // is_ale_
  AddtoPack(data,is_ale_);
  // Discretisation type
  AddtoPack(data,distype_);

  return;
}


/*--------------------------------------------------------------------------*
 |  Unpack data                                          (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);  ExtractfromPack(position,data,basedata);
  DRT::MESHFREE::Cell::Unpack(basedata);
  // is_ale_
  is_ale_ = ExtractInt(position,data);
  // distype
  distype_ = static_cast<DiscretizationType>( ExtractInt(position,data) );

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*--------------------------------------------------------------------------*
 |  Returns number of degrees of freedom at node         (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::NumDofPerNode(const DRT::Node& node) const
{
  // number of Dof's is fluid-specific.
  // Therefore, it is not pushed into drt_utilis_local_connectivity
  switch(distype_)
  {
  case DRT::Element::hex8:
  case DRT::Element::tet4:
    return 4;
    break;
  case DRT::Element::quad4:
  case DRT::Element::tri3:
    return 3;
    break;
  case DRT::Element::line2:
  case DRT::Element::line3:
    dserror("1D Fluid elements are not supported");
    break;
  default:
    dserror("discretization type %s not yet implemented", (DRT::DistypeToString(distype_)).c_str());
    break;
  }
  return 0;
}


/*--------------------------------------------------------------------------*
 |  Print this cell                                      (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluid::Print(std::ostream& os) const
{
  os << "MeshfreeFluid ";
  Print(os);
  return;
}

/*--------------------------------------------------------------------------*
 |  Read Input of this element                           (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
bool DRT::ELEMENTS::MeshfreeFluid::ReadElement(
  const std::string& eletype,
  const std::string& distype,
  DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(DRT::StringToDistype(distype));

  std::string na;
  linedef->ExtractString("NA",na);
  if (na=="ale" or na=="ALE" or na=="Ale")
  {
    is_ale_ = true;
  }
  else if (na=="euler" or na=="EULER" or na=="Euler")
    is_ale_ = false;
  else
    dserror("Reading of meshfree fluid element failed: Euler/Ale");

  return true;
}
