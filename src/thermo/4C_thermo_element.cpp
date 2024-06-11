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
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::ThermoType Discret::ELEMENTS::ThermoType::instance_;

Discret::ELEMENTS::ThermoType& Discret::ELEMENTS::ThermoType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::ThermoType::Create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Thermo* object = new Discret::ELEMENTS::Thermo(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Thermo(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Thermo(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 |                                                           dano 08/12 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ThermoType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}  // nodal_block_information()


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::ThermoType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ThermoBoundaryType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp(new Discret::ELEMENTS::ThermoBoundary(id,owner));
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | setup element                                             dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["THERMO"];

  defs["HEX8"] =
      Input::LineDefinition::Builder().AddIntVector("HEX8", 8).AddNamedInt("MAT").Build();

  defs["HEX20"] =
      Input::LineDefinition::Builder().AddIntVector("HEX20", 20).AddNamedInt("MAT").Build();

  defs["HEX27"] =
      Input::LineDefinition::Builder().AddIntVector("HEX27", 27).AddNamedInt("MAT").Build();

  defs["TET4"] =
      Input::LineDefinition::Builder().AddIntVector("TET4", 4).AddNamedInt("MAT").Build();

  defs["TET10"] =
      Input::LineDefinition::Builder().AddIntVector("TET10", 10).AddNamedInt("MAT").Build();

  defs["WEDGE6"] =
      Input::LineDefinition::Builder().AddIntVector("WEDGE6", 6).AddNamedInt("MAT").Build();

  defs["WEDGE15"] =
      Input::LineDefinition::Builder().AddIntVector("WEDGE15", 15).AddNamedInt("MAT").Build();

  defs["PYRAMID5"] =
      Input::LineDefinition::Builder().AddIntVector("PYRAMID5", 5).AddNamedInt("MAT").Build();

  defs["NURBS27"] =
      Input::LineDefinition::Builder().AddIntVector("NURBS27", 27).AddNamedInt("MAT").Build();

  defs["QUAD4"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddNamedInt("MAT").Build();

  defs["QUAD8"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddNamedInt("MAT").Build();

  defs["QUAD9"] =
      Input::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddNamedInt("MAT").Build();

  defs["TRI3"] =
      Input::LineDefinition::Builder().AddIntVector("TRI3", 3).AddNamedInt("MAT").Build();

  defs["TRI6"] =
      Input::LineDefinition::Builder().AddIntVector("TRI6", 6).AddNamedInt("MAT").Build();

  defs["NURBS4"] =
      Input::LineDefinition::Builder().AddIntVector("NURBS4", 4).AddNamedInt("MAT").Build();

  defs["NURBS9"] =
      Input::LineDefinition::Builder().AddIntVector("NURBS9", 9).AddNamedInt("MAT").Build();

  defs["LINE2"] =
      Input::LineDefinition::Builder().AddIntVector("LINE2", 2).AddNamedInt("MAT").Build();

  defs["LINE3"] =
      Input::LineDefinition::Builder().AddIntVector("LINE3", 3).AddNamedInt("MAT").Build();
}  // setup_element_definition()


Discret::ELEMENTS::ThermoBoundaryType Discret::ELEMENTS::ThermoBoundaryType::instance_;

Discret::ELEMENTS::ThermoBoundaryType& Discret::ELEMENTS::ThermoBoundaryType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Thermo::Thermo(int id, int owner)
    : Core::Elements::Element(id, owner), distype_(Core::FE::CellType::dis_none)
{
  // default: geometrically linear, also including purely thermal probelm
  kintype_ = Inpar::STR::KinemType::linear;
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Thermo::Thermo(const Discret::ELEMENTS::Thermo& old)
    : Core::Elements::Element(old), kintype_(old.kintype_), distype_(old.distype_)
{
  if (old.Shape() == Core::FE::CellType::nurbs27) SetNurbsElement() = true;
  return;
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance of Thermo and return              dano 09/09 |
 | pointer to it (public)                                               |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Thermo::Clone() const
{
  Discret::ELEMENTS::Thermo* newelement = new Discret::ELEMENTS::Thermo(*this);
  return newelement;
}  // Clone()


/*----------------------------------------------------------------------*
 | return the shape of a Thermo element (public)             dano 09/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Thermo::Shape() const { return distype_; }  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Thermo::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype
  add_to_pack(data, kintype_);
  // distype
  add_to_pack(data, distype_);

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Thermo::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<Inpar::STR::KinemType>(ExtractInt(position, data));
  // distype
  distype_ = static_cast<Core::FE::CellType>(ExtractInt(position, data));
  if (distype_ == Core::FE::CellType::nurbs27) SetNurbsElement() = true;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}  // Unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Thermo::Print(std::ostream& os) const
{
  os << "Thermo element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  return;
}  // Print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Thermo::Lines()
{
  return Core::Communication::GetElementLines<ThermoBoundary, Thermo>(*this);
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of surfaces (public)                           dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Thermo::Surfaces()
{
  return Core::Communication::GetElementSurfaces<ThermoBoundary, Thermo>(*this);
}  // Surfaces()

/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Thermo::VisNames(std::map<std::string, int>& names)
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
bool Discret::ELEMENTS::Thermo::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::VisData(name, data)) return true;

  return false;
}  // VisData()

/*----------------------------------------------------------------------------*
 | ENDE Discret::ELEMENTS::Thermo
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ThermoBoundary::ThermoBoundary(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Thermo* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ThermoBoundary::ThermoBoundary(const Discret::ELEMENTS::ThermoBoundary& old)
    : Core::Elements::FaceElement(old)
{
  return;
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance return pointer to it (public)     dano 09/09 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ThermoBoundary::Clone() const
{
  Discret::ELEMENTS::ThermoBoundary* newelement = new Discret::ELEMENTS::ThermoBoundary(*this);
  return newelement;
}  // Clone()


/*----------------------------------------------------------------------*
 | return shape of this element (public)                     dano 09/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::ThermoBoundary::Shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      if ((parent_element()->Shape() == Core::FE::CellType::quad8) or
          (parent_element()->Shape() == Core::FE::CellType::quad9))
        return Core::FE::CellType::line3;
      else
        return Core::FE::CellType::tri3;
    case 4:
      return Core::FE::CellType::quad4;
    case 6:
      return Core::FE::CellType::tri6;
    case 8:
      return Core::FE::CellType::quad8;
    case 9:
      if (parent_element()->Shape() == Core::FE::CellType::hex27)
        return Core::FE::CellType::quad9;
      else if (parent_element()->Shape() == Core::FE::CellType::nurbs27)
        return Core::FE::CellType::nurbs9;
      else
      {
        FOUR_C_THROW(
            "Your parent discretization type is %s. Ccurrently only hex27 and nurbs27 are "
            "implemented.",
            Core::FE::CellTypeToString(parent_element()->Shape()).c_str());
      }
      break;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ThermoBoundary::Pack(std::vector<char>& data) const
{
  FOUR_C_THROW("This ThermoBoundary element does not support communication");

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ThermoBoundary::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This ThermoBoundary element does not support communication");
  return;
}  // Unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ThermoBoundary::Print(std::ostream& os) const
{
  os << "ThermoBoundary ";
  Element::Print(os);
  return;
}  // Print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ThermoBoundary::Lines()
{
  FOUR_C_THROW("Lines of ThermoBoundary not implemented");
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ThermoBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of ThermoBoundary not implemented");
}  // Surfaces()

FOUR_C_NAMESPACE_CLOSE
