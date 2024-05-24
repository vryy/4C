/*----------------------------------------------------------------------*/
/*! \file
\brief basic thermo element

\level 1
*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "4C_thermo_element.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_element.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::ThermoType DRT::ELEMENTS::ThermoType::instance_;

DRT::ELEMENTS::ThermoType& DRT::ELEMENTS::ThermoType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Thermo* object = new DRT::ELEMENTS::Thermo(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Thermo(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Thermo(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 |                                                           dano 08/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ThermoType::nodal_block_information(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}  // nodal_block_information()


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::ThermoType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ThermoBoundaryType::Create(const int id, const int owner)
{
  // return Teuchos::rcp(new DRT::ELEMENTS::ThermoBoundary(id,owner));
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | setup element                                             dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["THERMO"];

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

  defs["NURBS27"] =
      INPUT::LineDefinition::Builder().AddIntVector("NURBS27", 27).AddNamedInt("MAT").Build();

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

  defs["NURBS4"] =
      INPUT::LineDefinition::Builder().AddIntVector("NURBS4", 4).AddNamedInt("MAT").Build();

  defs["NURBS9"] =
      INPUT::LineDefinition::Builder().AddIntVector("NURBS9", 9).AddNamedInt("MAT").Build();

  defs["LINE2"] =
      INPUT::LineDefinition::Builder().AddIntVector("LINE2", 2).AddNamedInt("MAT").Build();

  defs["LINE3"] =
      INPUT::LineDefinition::Builder().AddIntVector("LINE3", 3).AddNamedInt("MAT").Build();
}  // setup_element_definition()


DRT::ELEMENTS::ThermoBoundaryType DRT::ELEMENTS::ThermoBoundaryType::instance_;

DRT::ELEMENTS::ThermoBoundaryType& DRT::ELEMENTS::ThermoBoundaryType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Thermo::Thermo(int id, int owner)
    : DRT::Element(id, owner), distype_(CORE::FE::CellType::dis_none)
{
  // default: geometrically linear, also including purely thermal probelm
  kintype_ = INPAR::STR::KinemType::linear;
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Thermo::Thermo(const DRT::ELEMENTS::Thermo& old)
    : DRT::Element(old), kintype_(old.kintype_), distype_(old.distype_)
{
  if (old.Shape() == CORE::FE::CellType::nurbs27) SetNurbsElement() = true;
  return;
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance of Thermo and return              dano 09/09 |
 | pointer to it (public)                                               |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Thermo::Clone() const
{
  DRT::ELEMENTS::Thermo* newelement = new DRT::ELEMENTS::Thermo(*this);
  return newelement;
}  // Clone()


/*----------------------------------------------------------------------*
 | return the shape of a Thermo element (public)             dano 09/09 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Thermo::Shape() const { return distype_; }  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Thermo::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype
  AddtoPack(data, kintype_);
  // distype
  AddtoPack(data, distype_);

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Thermo::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));
  // distype
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  if (distype_ == CORE::FE::CellType::nurbs27) SetNurbsElement() = true;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}  // Unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Thermo::Print(std::ostream& os) const
{
  os << "Thermo element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  return;
}  // Print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Thermo::Lines()
{
  return CORE::COMM::GetElementLines<ThermoBoundary, Thermo>(*this);
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of surfaces (public)                           dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Thermo::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<ThermoBoundary, Thermo>(*this);
}  // Surfaces()

/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Thermo::VisNames(std::map<std::string, int>& names)
{
  // see whether we have additional data for visualization in our container
  for (int k = 0; k < numdofpernode_; k++)
  {
    std::ostringstream temp;
    temp << k;
  }  // loop over temperatures

  return;
}  // VisNames()


/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 09/09 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Thermo::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return false;
}  // VisData()

/*----------------------------------------------------------------------------*
 | ENDE DRT::ELEMENTS::Thermo
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ThermoBoundary::ThermoBoundary(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Thermo* parent, const int lsurface)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ThermoBoundary::ThermoBoundary(const DRT::ELEMENTS::ThermoBoundary& old)
    : DRT::FaceElement(old)
{
  return;
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance return pointer to it (public)     dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ThermoBoundary::Clone() const
{
  DRT::ELEMENTS::ThermoBoundary* newelement = new DRT::ELEMENTS::ThermoBoundary(*this);
  return newelement;
}  // Clone()


/*----------------------------------------------------------------------*
 | return shape of this element (public)                     dano 09/09 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::ThermoBoundary::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      if ((parent_element()->Shape() == CORE::FE::CellType::quad8) or
          (parent_element()->Shape() == CORE::FE::CellType::quad9))
        return CORE::FE::CellType::line3;
      else
        return CORE::FE::CellType::tri3;
    case 4:
      return CORE::FE::CellType::quad4;
    case 6:
      return CORE::FE::CellType::tri6;
    case 8:
      return CORE::FE::CellType::quad8;
    case 9:
      if (parent_element()->Shape() == CORE::FE::CellType::hex27)
        return CORE::FE::CellType::quad9;
      else if (parent_element()->Shape() == CORE::FE::CellType::nurbs27)
        return CORE::FE::CellType::nurbs9;
      else
      {
        FOUR_C_THROW(
            "Your parent discretization type is %s. Ccurrently only hex27 and nurbs27 are "
            "implemented.",
            CORE::FE::CellTypeToString(parent_element()->Shape()).c_str());
      }
      break;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ThermoBoundary::Pack(std::vector<char>& data) const
{
  FOUR_C_THROW("This ThermoBoundary element does not support communication");

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ThermoBoundary::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This ThermoBoundary element does not support communication");
  return;
}  // Unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ThermoBoundary::Print(std::ostream& os) const
{
  os << "ThermoBoundary ";
  Element::Print(os);
  return;
}  // Print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ThermoBoundary::Lines()
{
  FOUR_C_THROW("Lines of ThermoBoundary not implemented");
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ThermoBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of ThermoBoundary not implemented");
}  // Surfaces()

FOUR_C_NAMESPACE_CLOSE
