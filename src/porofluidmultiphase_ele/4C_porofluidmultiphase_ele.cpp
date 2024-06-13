/*----------------------------------------------------------------------*/
/*! \file
 \brief definition of porofluidmultiphase elements

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 ******************  PoroFluidMultiPhase ElementType *********************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | instantiate global instance                               vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseType Discret::ELEMENTS::PoroFluidMultiPhaseType::instance_;

/*----------------------------------------------------------------------*
 | instance access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseType& Discret::ELEMENTS::PoroFluidMultiPhaseType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create an element from data                              vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::PoroFluidMultiPhaseType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::PoroFluidMultiPhase* object =
      new Discret::ELEMENTS::PoroFluidMultiPhase(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |  create an element from a dat file specifier             vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::PoroFluidMultiPhaseType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "POROFLUIDMULTIPHASE")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::PoroFluidMultiPhase(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  create an empty element                                vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::PoroFluidMultiPhaseType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::PoroFluidMultiPhase(id, owner));
  return ele;
}

/*----------------------------------------------------------------------------*
 |  nodal block information to create a null space description    vuong 08/16 |
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhaseType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}
/*----------------------------------------------------------------------*
 |  do the null space computation                            vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::PoroFluidMultiPhaseType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

/*----------------------------------------------------------------------*
 |  setup the dat file input line definitions for this type of element   |
 |                                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhaseType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["POROFLUIDMULTIPHASE"];

  defs["QUAD4"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD4", 4).add_named_int("MAT").Build();

  defs["QUAD8"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD8", 8).add_named_int("MAT").Build();

  defs["QUAD9"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD9", 9).add_named_int("MAT").Build();

  defs["TRI3"] =
      Input::LineDefinition::Builder().add_int_vector("TRI3", 3).add_named_int("MAT").Build();

  defs["TRI6"] =
      Input::LineDefinition::Builder().add_int_vector("TRI6", 6).add_named_int("MAT").Build();

  defs["LINE2"] =
      Input::LineDefinition::Builder().add_int_vector("LINE2", 2).add_named_int("MAT").Build();

  defs["LINE3"] =
      Input::LineDefinition::Builder().add_int_vector("LINE3", 3).add_named_int("MAT").Build();

  defs["HEX8"] =
      Input::LineDefinition::Builder().add_int_vector("HEX8", 8).add_named_int("MAT").Build();

  defs["TET4"] =
      Input::LineDefinition::Builder().add_int_vector("TET4", 4).add_named_int("MAT").Build();

  defs["TET10"] =
      Input::LineDefinition::Builder().add_int_vector("TET10", 10).add_named_int("MAT").Build();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*
 ******************  PoroFluidMultiPhase BoundaryType ********************
 *---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | instantiate global instance                               vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryType
    Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryType::instance_;

/*----------------------------------------------------------------------*
 | instance access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryType&
Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  create an empty element                                vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryType::Create(
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
int Discret::ELEMENTS::PoroFluidMultiPhaseType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    Discret::ELEMENTS::PoroFluidMultiPhase* actele =
        dynamic_cast<Discret::ELEMENTS::PoroFluidMultiPhase*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to PoroFluidMultiPhase* failed");
    actele->initialize();
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
Discret::ELEMENTS::PoroFluidMultiPhase::PoroFluidMultiPhase(int id, int owner)
    : Core::Elements::Element(id, owner), distype_(Core::FE::CellType::dis_none), numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhase::PoroFluidMultiPhase(
    const Discret::ELEMENTS::PoroFluidMultiPhase& old)
    : Core::Elements::Element(old), distype_(old.distype_), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  Deep copy this instance of PoroFluidMultiPhase and return pointer to it  |
 |                                                 (public) vuong 08/16      |
 *---------------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::PoroFluidMultiPhase::Clone() const
{
  Discret::ELEMENTS::PoroFluidMultiPhase* newelement =
      new Discret::ELEMENTS::PoroFluidMultiPhase(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a PoroFluidMultiPhase element          (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::PoroFluidMultiPhase::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Initialize element                                      (protected) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhase::initialize()
{
  Teuchos::RCP<Mat::FluidPoroMultiPhase> actmat =
      Teuchos::rcp_dynamic_cast<Mat::FluidPoroMultiPhase>(Material(), true);

  actmat->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhase::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  Element::Pack(data);

  // add internal data
  add_to_pack(data, distype_);
  add_to_pack(data, numdofpernode_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);

  // extract internal data
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));
  extract_from_pack(position, data, numdofpernode_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)         vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::PoroFluidMultiPhase::NumLine() const
{
  return Core::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::PoroFluidMultiPhase::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::PoroFluidMultiPhase::NumVolume() const
{
  return Core::FE::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhase::Print(std::ostream& os) const
{
  os << "PoroFluidMultiPhase element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                 vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::PoroFluidMultiPhase::Lines()
{
  return Core::Communication::GetElementLines<PoroFluidMultiPhaseBoundary, PoroFluidMultiPhase>(
      *this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::PoroFluidMultiPhase::Surfaces()
{
  return Core::Communication::GetElementSurfaces<PoroFluidMultiPhaseBoundary, PoroFluidMultiPhase>(
      *this);
}

/*----------------------------------------------------------------------*
 | read element input                                       vuong 08/16 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::PoroFluidMultiPhase::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // set discretization type
  SetDisType(Core::FE::StringToCellType(distype));

  return true;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhase::SetMaterial(
    const int index, Teuchos::RCP<Core::Mat::Material> mat)
{
  // the standard part:
  Core::Elements::Element::SetMaterial(index, mat);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  if (mat->MaterialType() == Core::Materials::m_fluidporo_multiphase or
      mat->MaterialType() == Core::Materials::m_fluidporo_multiphase_reactions)
  {
    const Mat::FluidPoroMultiPhase* actmat =
        dynamic_cast<const Mat::FluidPoroMultiPhase*>(mat.get());
    if (actmat == nullptr) FOUR_C_THROW("cast failed");
    numdofpernode_ = actmat->NumMat();
  }
  else
    FOUR_C_THROW(
        "PoroFluidMultiPhase element got unsupported material type %d", mat->MaterialType());

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
Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::PoroFluidMultiPhaseBoundary(int id, int owner,
    int nnode, const int* nodeids, Core::Nodes::Node** nodes,
    Discret::ELEMENTS::PoroFluidMultiPhase* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::PoroFluidMultiPhaseBoundary(
    const Discret::ELEMENTS::PoroFluidMultiPhaseBoundary& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it   (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Clone() const
{
  Discret::ELEMENTS::PoroFluidMultiPhaseBoundary* newelement =
      new Discret::ELEMENTS::PoroFluidMultiPhaseBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Shape() const
{
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_element()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Pack(
    Core::Communication::PackBuffer& data) const
{
  // boundary elements are rebuild by their parent element for each condition
  // after redistribution. This way we make sure, that the node ids always match.
  // -> no communication of boundary elements
  FOUR_C_THROW("This PoroFluidMultiPhaseBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Unpack(const std::vector<char>& data)
{
  // boundary elements are rebuild by their parent element for each condition
  // after redistribution. This way we make sure, that the node ids always match.
  // -> no communication of boundary elements
  FOUR_C_THROW("This PoroFluidMultiPhaseBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Print(std::ostream& os) const
{
  os << "PoroFluidMultiPhaseBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::NumLine() const
{
  return Core::FE::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  vuong 08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Lines()
{
  FOUR_C_THROW("Lines of PoroFluidMultiPhaseBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            vuong 08/16 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::PoroFluidMultiPhaseBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of PoroFluidMultiPhaseBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
