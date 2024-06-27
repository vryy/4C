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

Discret::ELEMENTS::LubricationType& Discret::ELEMENTS::LubricationType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::LubricationType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Lubrication* object = new Discret::ELEMENTS::Lubrication(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationType::Create(
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


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Lubrication(id, owner));
  return ele;
}


void Discret::ELEMENTS::LubricationType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::LubricationType::ComputeNullSpace(
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

Discret::ELEMENTS::LubricationBoundaryType& Discret::ELEMENTS::LubricationBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::LubricationBoundaryType::Create(
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
Core::Elements::Element* Discret::ELEMENTS::Lubrication::Clone() const
{
  Discret::ELEMENTS::Lubrication* newelement = new Discret::ELEMENTS::Lubrication(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Lubrication element                     (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Lubrication::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Lubrication::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
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

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

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
int Discret::ELEMENTS::Lubrication::NumLine() const
{
  return Core::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::NumVolume() const
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
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Lubrication::Lines()
{
  return Core::Communication::GetElementLines<LubricationBoundary, Lubrication>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Lubrication::Surfaces()
{
  return Core::Communication::GetElementSurfaces<LubricationBoundary, Lubrication>(*this);
}

/*----------------------------------------------------------------------*
 | read element input                                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Lubrication::ReadElement(
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
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
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
Core::Elements::Element* Discret::ELEMENTS::LubricationBoundary::Clone() const
{
  Discret::ELEMENTS::LubricationBoundary* newelement =
      new Discret::ELEMENTS::LubricationBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::LubricationBoundary::Shape() const
{
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_element()->Shape());
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
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::LubricationBoundary::NumLine() const
{
  return Core::FE::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::LubricationBoundary::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::LubricationBoundary::Lines()
{
  FOUR_C_THROW("Lines of LubricationBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::LubricationBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of LubricationBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
