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

Thermo::ElementType Thermo::ElementType::instance_;

Thermo::ElementType& Thermo::ElementType::instance() { return instance_; }

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Thermo::ElementType::create(const std::vector<char>& data)
{
  Thermo::Element* object = new Thermo::Element(-1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Thermo::ElementType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(new Thermo::Element(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 09/09 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Thermo::ElementType::create(const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(new Thermo::Element(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 |                                                           dano 08/12 |
 *----------------------------------------------------------------------*/
void Thermo::ElementType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->num_dof_per_node(*(dwele->nodes()[0]));
  dimns = numdf;
  nv = numdf;
}  // nodal_block_information()


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Thermo::ElementType::compute_null_space(
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
Teuchos::RCP<Core::Elements::Element> Thermo::FaceElementType::create(const int id, const int owner)
{
  // return Teuchos::rcp(new FaceElement(id,owner));
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | setup element                                             dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::ElementType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["THERMO"];

  defs["HEX8"] =
      Input::LineDefinition::Builder().add_int_vector("HEX8", 8).add_named_int("MAT").build();

  defs["HEX20"] =
      Input::LineDefinition::Builder().add_int_vector("HEX20", 20).add_named_int("MAT").build();

  defs["HEX27"] =
      Input::LineDefinition::Builder().add_int_vector("HEX27", 27).add_named_int("MAT").build();

  defs["TET4"] =
      Input::LineDefinition::Builder().add_int_vector("TET4", 4).add_named_int("MAT").build();

  defs["TET10"] =
      Input::LineDefinition::Builder().add_int_vector("TET10", 10).add_named_int("MAT").build();

  defs["WEDGE6"] =
      Input::LineDefinition::Builder().add_int_vector("WEDGE6", 6).add_named_int("MAT").build();

  defs["WEDGE15"] =
      Input::LineDefinition::Builder().add_int_vector("WEDGE15", 15).add_named_int("MAT").build();

  defs["PYRAMID5"] =
      Input::LineDefinition::Builder().add_int_vector("PYRAMID5", 5).add_named_int("MAT").build();

  defs["NURBS27"] =
      Input::LineDefinition::Builder().add_int_vector("NURBS27", 27).add_named_int("MAT").build();

  defs["QUAD4"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD4", 4).add_named_int("MAT").build();

  defs["QUAD8"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD8", 8).add_named_int("MAT").build();

  defs["QUAD9"] =
      Input::LineDefinition::Builder().add_int_vector("QUAD9", 9).add_named_int("MAT").build();

  defs["TRI3"] =
      Input::LineDefinition::Builder().add_int_vector("TRI3", 3).add_named_int("MAT").build();

  defs["TRI6"] =
      Input::LineDefinition::Builder().add_int_vector("TRI6", 6).add_named_int("MAT").build();

  defs["NURBS4"] =
      Input::LineDefinition::Builder().add_int_vector("NURBS4", 4).add_named_int("MAT").build();

  defs["NURBS9"] =
      Input::LineDefinition::Builder().add_int_vector("NURBS9", 9).add_named_int("MAT").build();

  defs["LINE2"] =
      Input::LineDefinition::Builder().add_int_vector("LINE2", 2).add_named_int("MAT").build();

  defs["LINE3"] =
      Input::LineDefinition::Builder().add_int_vector("LINE3", 3).add_named_int("MAT").build();
}  // setup_element_definition()


Thermo::FaceElementType Thermo::FaceElementType::instance_;

Thermo::FaceElementType& Thermo::FaceElementType::instance() { return instance_; }

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
Thermo::Element::Element(int id, int owner)
    : Core::Elements::Element(id, owner), distype_(Core::FE::CellType::dis_none)
{
  // default: geometrically linear, also including purely thermal probelm
  kintype_ = Inpar::Solid::KinemType::linear;
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
Thermo::Element::Element(const Element& old)
    : Core::Elements::Element(old), kintype_(old.kintype_), distype_(old.distype_)
{
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance of Thermo and return              dano 09/09 |
 | pointer to it (public)                                               |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Thermo::Element::clone() const
{
  Thermo::Element* newelement = new Thermo::Element(*this);
  return newelement;
}  // clone()


/*----------------------------------------------------------------------*
 | return the shape of a Thermo element (public)             dano 09/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Thermo::Element::shape() const { return distype_; }  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::Element::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Core::Elements::Element::pack(data);
  // kintype
  add_to_pack(data, kintype_);
  // distype
  add_to_pack(data, distype_);

  return;
}  // pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::Element::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::extract_and_assert_id(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::unpack(basedata);
  // kintype_
  kintype_ = static_cast<Inpar::Solid::KinemType>(extract_int(position, data));
  // distype
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}  // unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::Element::print(std::ostream& os) const
{
  os << "Thermo element";
  Core::Elements::Element::print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::cell_type_to_string(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  return;
}  // print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Thermo::Element::lines()
{
  return Core::Communication::get_element_lines<FaceElement, Element>(*this);
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of surfaces (public)                           dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Thermo::Element::surfaces()
{
  return Core::Communication::get_element_surfaces<FaceElement, Element>(*this);
}  // Surfaces()

/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::Element::vis_names(std::map<std::string, int>& names)
{
  // see whether we have additional data for visualization in our container
  for (int k = 0; k < numdofpernode_; k++)
  {
    std::ostringstream temp;
    temp << k;
  }  // loop over temperatures

  return;
}  // vis_names()


/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 09/09 |
 *----------------------------------------------------------------------*/
bool Thermo::Element::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return false;
}  // vis_data()

/*----------------------------------------------------------------------------*
 | ENDE Thermo::Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 09/09 |
 *----------------------------------------------------------------------*/
Thermo::FaceElement::FaceElement(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Element* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}  // ctor


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
Thermo::FaceElement::FaceElement(const FaceElement& old) : Core::Elements::FaceElement(old)
{
  return;
}  // copy-ctor


/*----------------------------------------------------------------------*
 | deep copy this instance return pointer to it (public)     dano 09/09 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Thermo::FaceElement::clone() const
{
  FaceElement* newelement = new FaceElement(*this);
  return newelement;
}  // clone()


/*----------------------------------------------------------------------*
 | return shape of this element (public)                     dano 09/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Thermo::FaceElement::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      if ((parent_element()->shape() == Core::FE::CellType::quad8) or
          (parent_element()->shape() == Core::FE::CellType::quad9))
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
      if (parent_element()->shape() == Core::FE::CellType::hex27)
        return Core::FE::CellType::quad9;
      else if (parent_element()->shape() == Core::FE::CellType::nurbs27)
        return Core::FE::CellType::nurbs9;
      else
      {
        FOUR_C_THROW(
            "Your parent discretization type is %s. Ccurrently only hex27 and nurbs27 are "
            "implemented.",
            Core::FE::cell_type_to_string(parent_element()->shape()).c_str());
      }
      break;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}  // Shape()


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::FaceElement::pack(std::vector<char>& data) const
{
  FOUR_C_THROW("This FaceElement element does not support communication");

  return;
}  // pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::FaceElement::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This FaceElement element does not support communication");
  return;
}  // unpack()



/*----------------------------------------------------------------------*
 | print this element (public)                               dano 09/09 |
 *----------------------------------------------------------------------*/
void Thermo::FaceElement::print(std::ostream& os) const
{
  os << "FaceElement ";
  Element::print(os);
  return;
}  // print()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Thermo::FaceElement::lines()
{
  FOUR_C_THROW("Lines of FaceElement not implemented");
}  // Lines()


/*----------------------------------------------------------------------*
 | get vector of lines (public)                              dano 09/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Thermo::FaceElement::surfaces()
{
  FOUR_C_THROW("Surfaces of FaceElement not implemented");
}  // Surfaces()

FOUR_C_NAMESPACE_CLOSE
