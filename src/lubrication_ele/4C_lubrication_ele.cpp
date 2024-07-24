/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::LubricationType Discret::ELEMENTS::LubricationType::instance_;

Discret::ELEMENTS::LubricationType& Discret::ELEMENTS::LubricationType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::LubricationType::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Lubrication* object = new Discret::ELEMENTS::Lubrication(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "LUBRICATION")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Lubrication(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationType::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Lubrication(id, owner));
  return ele;
}


void Discret::ELEMENTS::LubricationType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->num_dof_per_node(*(dwele->nodes()[0]));
  dimns = numdf;
  nv = numdf;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::LubricationType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void Discret::ELEMENTS::LubricationType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["LUBRICATION"];

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

  defs["LINE2"] =
      Input::LineDefinition::Builder().add_int_vector("LINE2", 2).add_named_int("MAT").build();

  defs["LINE3"] =
      Input::LineDefinition::Builder().add_int_vector("LINE3", 3).add_named_int("MAT").build();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::LubricationBoundaryType Discret::ELEMENTS::LubricationBoundaryType::instance_;

Discret::ELEMENTS::LubricationBoundaryType& Discret::ELEMENTS::LubricationBoundaryType::instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationBoundaryType::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new LubricationBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Lubrication::Lubrication(int id, int owner)
    : Core::Elements::Element(id, owner), distype_(Core::FE::CellType::dis_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Lubrication::Lubrication(const Discret::ELEMENTS::Lubrication& old)
    : Core::Elements::Element(old), distype_(old.distype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Lubrication and return pointer to it        |
 |                                                 (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Lubrication::clone() const
{
  Discret::ELEMENTS::Lubrication* newelement = new Discret::ELEMENTS::Lubrication(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Lubrication element                     (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Lubrication::shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Lubrication::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);

  // add internal data
  add_to_pack(data, distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Lubrication::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

  // extract internal data
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)         wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::num_line() const
{
  return Core::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::num_surface() const
{
  return Core::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::num_volume() const
{
  return Core::FE::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Lubrication::print(std::ostream& os) const
{
  os << "Lubrication element";
  Element::print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(distype_) << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Lubrication::lines()
{
  return Core::Communication::GetElementLines<LubricationBoundary, Lubrication>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Lubrication::surfaces()
{
  return Core::Communication::GetElementSurfaces<LubricationBoundary, Lubrication>(*this);
}

/*----------------------------------------------------------------------*
 | read element input                                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Lubrication::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  // set discretization type
  set_dis_type(Core::FE::StringToCellType(distype));

  return true;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::LubricationBoundary::LubricationBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::ELEMENTS::Lubrication* parent,
    const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::LubricationBoundary::LubricationBoundary(
    const Discret::ELEMENTS::LubricationBoundary& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::LubricationBoundary::clone() const
{
  Discret::ELEMENTS::LubricationBoundary* newelement =
      new Discret::ELEMENTS::LubricationBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::LubricationBoundary::shape() const
{
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_element()->shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationBoundary::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("This LubricationBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationBoundary::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This LubricationBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationBoundary::print(std::ostream& os) const
{
  os << "LubricationBoundary element";
  Element::print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::LubricationBoundary::num_line() const
{
  return Core::FE::getNumberOfElementLines(shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::LubricationBoundary::num_surface() const
{
  return Core::FE::getNumberOfElementSurfaces(shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::LubricationBoundary::lines()
{
  FOUR_C_THROW("Lines of LubricationBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::LubricationBoundary::surfaces()
{
  FOUR_C_THROW("Surfaces of LubricationBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
