/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::LubricationType DRT::ELEMENTS::LubricationType::instance_;

DRT::ELEMENTS::LubricationType& DRT::ELEMENTS::LubricationType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::LubricationType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Lubrication* object = new DRT::ELEMENTS::Lubrication(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::LubricationType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "LUBRICATION")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Lubrication(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::LubricationType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Lubrication(id, owner));
  return ele;
}


void DRT::ELEMENTS::LubricationType::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::LubricationType::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void DRT::ELEMENTS::LubricationType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["LUBRICATION"];

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::LubricationBoundaryType DRT::ELEMENTS::LubricationBoundaryType::instance_;

DRT::ELEMENTS::LubricationBoundaryType& DRT::ELEMENTS::LubricationBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::LubricationBoundaryType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new LubricationBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Lubrication::Lubrication(int id, int owner)
    : CORE::Elements::Element(id, owner), distype_(CORE::FE::CellType::dis_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Lubrication::Lubrication(const DRT::ELEMENTS::Lubrication& old)
    : CORE::Elements::Element(old), distype_(old.distype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Lubrication and return pointer to it        |
 |                                                 (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Lubrication::Clone() const
{
  DRT::ELEMENTS::Lubrication* newelement = new DRT::ELEMENTS::Lubrication(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Lubrication element                     (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Lubrication::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Lubrication::Pack(CORE::COMM::PackBuffer& data) const
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

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Lubrication::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // extract internal data
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)         wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Lubrication::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Lubrication::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Lubrication::NumVolume() const
{
  return CORE::FE::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Lubrication::Print(std::ostream& os) const
{
  os << "Lubrication element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(distype_) << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Lubrication::Lines()
{
  return CORE::COMM::GetElementLines<LubricationBoundary, Lubrication>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Lubrication::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<LubricationBoundary, Lubrication>(*this);
}

/*----------------------------------------------------------------------*
 | read element input                                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Lubrication::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, MAT::Factory(material));

  // set discretization type
  SetDisType(CORE::FE::StringToCellType(distype));

  return true;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationBoundary::LubricationBoundary(int id, int owner, int nnode,
    const int* nodeids, CORE::Nodes::Node** nodes, DRT::ELEMENTS::Lubrication* parent,
    const int lsurface)
    : CORE::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationBoundary::LubricationBoundary(
    const DRT::ELEMENTS::LubricationBoundary& old)
    : CORE::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::LubricationBoundary::Clone() const
{
  DRT::ELEMENTS::LubricationBoundary* newelement = new DRT::ELEMENTS::LubricationBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::LubricationBoundary::Shape() const
{
  return CORE::FE::getShapeOfBoundaryElement(num_node(), parent_element()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationBoundary::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("This LubricationBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationBoundary::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This LubricationBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationBoundary::Print(std::ostream& os) const
{
  os << "LubricationBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::LubricationBoundary::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::LubricationBoundary::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::LubricationBoundary::Lines()
{
  FOUR_C_THROW("Lines of LubricationBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::LubricationBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of LubricationBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
