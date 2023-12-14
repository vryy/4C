/*----------------------------------------------------------------------*/
/*! \file
 \brief definition of porofluidmultiphase elements

   \level 3

 *----------------------------------------------------------------------*/


#include "baci_porofluidmultiphase_ele.H"

#include "baci_comm_utils_factory.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_fluid_ele_nullspace.H"
#include "baci_io_linedefinition.H"
#include "baci_mat_fluidporo_multiphase.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 ******************  PoroFluidMultiPhase ElementType *********************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | instantiate global instance                               vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseType DRT::ELEMENTS::PoroFluidMultiPhaseType::instance_;

/*----------------------------------------------------------------------*
 | instance access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseType& DRT::ELEMENTS::PoroFluidMultiPhaseType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create an element from data                              vuong 08/16 |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::PoroFluidMultiPhaseType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::PoroFluidMultiPhase* object = new DRT::ELEMENTS::PoroFluidMultiPhase(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |  create an element from a dat file specifier             vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::PoroFluidMultiPhaseType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "POROFLUIDMULTIPHASE")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::PoroFluidMultiPhase(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  create an empty element                                vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::PoroFluidMultiPhaseType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::PoroFluidMultiPhase(id, owner));
  return ele;
}

/*----------------------------------------------------------------------------*
 |  nodal block information to create a null space description    vuong 08/16 |
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}
/*----------------------------------------------------------------------*
 |  do the null space computation                            vuong 08/16 |
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::PoroFluidMultiPhaseType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

/*----------------------------------------------------------------------*
 |  setup the dat file input line definitions for this type of element   |
 |                                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["POROFLUIDMULTIPHASE"];

  defs["QUAD4"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddNamedInt("MAT").Build();

  defs["QUAD8"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddNamedInt("MAT").Build();

  defs["QUAD9"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddNamedInt("MAT").Build();

  defs["TRI3"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI3", 3).AddNamedInt("MAT").Build();

  defs["TRI6"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI6", 6).AddNamedInt("MAT").Build();

  defs["LINE2"] =
      INPUT::LineDefinition::Builder().AddIntVector("LINE2", 2).AddNamedInt("MAT").Build();

  defs["LINE3"] =
      INPUT::LineDefinition::Builder().AddIntVector("LINE3", 3).AddNamedInt("MAT").Build();

  defs["HEX8"] =
      INPUT::LineDefinition::Builder().AddIntVector("HEX8", 8).AddNamedInt("MAT").Build();

  defs["TET4"] =
      INPUT::LineDefinition::Builder().AddIntVector("TET4", 4).AddNamedInt("MAT").Build();

  defs["TET10"] =
      INPUT::LineDefinition::Builder().AddIntVector("TET10", 10).AddNamedInt("MAT").Build();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 ******************  PoroFluidMultiPhase BoundaryType ********************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | instantiate global instance                               vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryType
    DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryType::instance_;

/*----------------------------------------------------------------------*
 | instance access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryType&
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  create an empty element                                vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryType::Create(
    const int id, const int owner)
{
  // boundary elements are not created as stand-alone elements by the element type,
  // but they are build directly by the corresponding domain element instead.
  // See the ElementBoundaryFactory.
  // Hence, boundary type classes actually are not used as factories, but
  // only for type identification.
  // To make this clear, a null pointer is returned.

  // return Teuchos::rcp( new PoroFluidMultiPhaseBoundary( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::PoroFluidMultiPhase* actele =
        dynamic_cast<DRT::ELEMENTS::PoroFluidMultiPhase*>(dis.lColElement(i));
    if (!actele) dserror("cast to PoroFluidMultiPhase* failed");
    actele->Initialize();
  }
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 *********************  PoroFluidMultiPhase Element **********************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  create an empty element                                vuong 08/16 |
 *------------------------------------------------------ ----------------*/
DRT::ELEMENTS::PoroFluidMultiPhase::PoroFluidMultiPhase(int id, int owner)
    : DRT::Element(id, owner), distype_(CORE::FE::CellType::dis_none), numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhase::PoroFluidMultiPhase(
    const DRT::ELEMENTS::PoroFluidMultiPhase& old)
    : DRT::Element(old), distype_(old.distype_), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  Deep copy this instance of PoroFluidMultiPhase and return pointer to it  |
 |                                                 (public) vuong 08/16      |
 *---------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::PoroFluidMultiPhase::Clone() const
{
  DRT::ELEMENTS::PoroFluidMultiPhase* newelement = new DRT::ELEMENTS::PoroFluidMultiPhase(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a PoroFluidMultiPhase element          (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::PoroFluidMultiPhase::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Initialize element                                      (protected) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhase::Initialize()
{
  Teuchos::RCP<MAT::FluidPoroMultiPhase> actmat =
      Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(Material(), true);

  actmat->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhase::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);

  // add internal data
  AddtoPack(data, distype_);
  AddtoPack(data, numdofpernode_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // extract internal data
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  ExtractfromPack(position, data, numdofpernode_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)         vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhase::NumLine() const
{
  return CORE::DRT::UTILS::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhase::NumSurface() const
{
  return CORE::DRT::UTILS::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhase::NumVolume() const
{
  return CORE::DRT::UTILS::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhase::Print(std::ostream& os) const
{
  os << "PoroFluidMultiPhase element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                 vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::PoroFluidMultiPhase::Lines()
{
  return CORE::COMM::GetElementLines<PoroFluidMultiPhaseBoundary, PoroFluidMultiPhase>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::PoroFluidMultiPhase::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<PoroFluidMultiPhaseBoundary, PoroFluidMultiPhase>(*this);
}

/*----------------------------------------------------------------------*
 | read element input                                       vuong 08/16 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::PoroFluidMultiPhase::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set discretization type
  SetDisType(CORE::FE::StringToCellType(distype));

  return true;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhase::SetMaterial(int matnum)
{
  // the standard part:
  DRT::Element::SetMaterial(matnum);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType() == INPAR::MAT::m_fluidporo_multiphase or
      mat->MaterialType() == INPAR::MAT::m_fluidporo_multiphase_reactions)
  {
    const MAT::FluidPoroMultiPhase* actmat =
        dynamic_cast<const MAT::FluidPoroMultiPhase*>(mat.get());
    if (actmat == nullptr) dserror("cast failed");
    numdofpernode_ = actmat->NumMat();
  }
  else
    dserror("PoroFluidMultiPhase element got unsupported material type %d", mat->MaterialType());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 ***************  PoroFluidMultiPhase Boundary Element *******************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  ctor (public)                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::PoroFluidMultiPhaseBoundary(int id, int owner,
    int nnode, const int* nodeids, DRT::Node** nodes, DRT::ELEMENTS::PoroFluidMultiPhase* parent,
    const int lsurface)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::PoroFluidMultiPhaseBoundary(
    const DRT::ELEMENTS::PoroFluidMultiPhaseBoundary& old)
    : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it   (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Clone() const
{
  DRT::ELEMENTS::PoroFluidMultiPhaseBoundary* newelement =
      new DRT::ELEMENTS::PoroFluidMultiPhaseBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Shape() const
{
  return CORE::DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Pack(CORE::COMM::PackBuffer& data) const
{
  // boundary elements are rebuild by their parent element for each condition
  // after redistribution. This way we make sure, that the node ids always match.
  // -> no communication of boundary elements
  dserror("This PoroFluidMultiPhaseBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Unpack(const std::vector<char>& data)
{
  // boundary elements are rebuild by their parent element for each condition
  // after redistribution. This way we make sure, that the node ids always match.
  // -> no communication of boundary elements
  dserror("This PoroFluidMultiPhaseBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Print(std::ostream& os) const
{
  os << "PoroFluidMultiPhaseBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::NumLine() const
{
  return CORE::DRT::UTILS::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::NumSurface() const
{
  return CORE::DRT::UTILS::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Lines()
{
  dserror("Lines of PoroFluidMultiPhaseBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Surfaces()
{
  dserror("Surfaces of PoroFluidMultiPhaseBoundary not implemented");
}

BACI_NAMESPACE_CLOSE
