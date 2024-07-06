/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the main fluid element


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_fluid_ele_tds.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidType Discret::ELEMENTS::FluidType::instance_;

Discret::ELEMENTS::FluidType& Discret::ELEMENTS::FluidType::instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::FluidType::create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Fluid* object = new Discret::ELEMENTS::Fluid(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUID")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::Fluid(id, owner));
  }
  else if (eletype == "FLUID2" || eletype == "FLUID3")
  {
    FOUR_C_THROW("Fluid element types FLUID2 and FLUID3 are no longer in use. Switch to FLUID.");
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidType::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Fluid(id, owner));
}


void Discret::ELEMENTS::FluidType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->num_dof_per_node(*(dwele->nodes()[0]));
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}


Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::FluidType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void Discret::ELEMENTS::FluidType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["FLUID"];

  defsgeneral["HEX8"] = Input::LineDefinition::Builder()
                            .add_int_vector("HEX8", 8)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["HEX20"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX20", 20)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["HEX27"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX27", 27)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["TET4"] = Input::LineDefinition::Builder()
                            .add_int_vector("TET4", 4)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["TET10"] = Input::LineDefinition::Builder()
                             .add_int_vector("TET10", 10)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["WEDGE6"] = Input::LineDefinition::Builder()
                              .add_int_vector("WEDGE6", 6)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["WEDGE15"] = Input::LineDefinition::Builder()
                               .add_int_vector("WEDGE15", 15)
                               .add_named_int("MAT")
                               .add_named_string("NA")
                               .build();

  defsgeneral["PYRAMID5"] = Input::LineDefinition::Builder()
                                .add_int_vector("PYRAMID5", 5)
                                .add_named_int("MAT")
                                .add_named_string("NA")
                                .build();

  defsgeneral["NURBS8"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS8", 8)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["NURBS27"] = Input::LineDefinition::Builder()
                               .add_int_vector("NURBS27", 27)
                               .add_named_int("MAT")
                               .add_named_string("NA")
                               .build();

  // 2D elements
  defsgeneral["QUAD4"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD4", 4)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["QUAD8"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD8", 8)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["QUAD9"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD9", 9)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["TRI3"] = Input::LineDefinition::Builder()
                            .add_int_vector("TRI3", 3)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["TRI6"] = Input::LineDefinition::Builder()
                            .add_int_vector("TRI6", 6)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["NURBS4"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS4", 4)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["NURBS9"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS9", 9)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Fluid::Fluid(int id, int owner)
    : Core::Elements::Element(id, owner), is_ale_(false)
{
  distype_ = Core::FE::CellType::dis_none;
  tds_ = Teuchos::null;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Fluid::Fluid(const Discret::ELEMENTS::Fluid& old)
    : Core::Elements::Element(old), distype_(old.distype_), is_ale_(old.is_ale_)
{
  tds_ = Teuchos::null;
  if (old.tds_ != Teuchos::null)
    FOUR_C_THROW("clone() method for deep copying tds_ not yet implemented!");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Fluid::clone() const
{
  Discret::ELEMENTS::Fluid* newelement = new Discret::ELEMENTS::Fluid(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Fluid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);
  // is_ale_
  add_to_pack(data, is_ale_);
  // Discretisation type
  add_to_pack(data, distype_);

  // time-dependent subgrid scales
  bool is_tds(false);
  if (tds_ != Teuchos::null)
  {
    is_tds = true;
    add_to_pack(data, is_tds);
    tds_->pack(data);
  }
  else
  {
    add_to_pack(data, is_tds);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Fluid::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);
  // is_ale_
  is_ale_ = extract_int(position, data);
  // distype
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));

  // time-dependent subgrid scales
  bool is_tds = extract_int(position, data);
  if (is_tds)
  {
    tds_ = Teuchos::rcp(new FLD::TDSEleData());
    std::vector<char> pbtest;
    extract_from_pack(position, data, pbtest);
    if (pbtest.size() == 0) FOUR_C_THROW("Seems no TDS data available");
    tds_->unpack(pbtest);
  }
  else
    tds_ = Teuchos::null;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Fluid::print(std::ostream& os) const
{
  os << "Fluid ";
  Element::print(os);
  // cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                 ae  02/010|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Fluid::lines()
{
  return Core::Communication::GetElementLines<FluidBoundary, Fluid>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          ehrl  02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Fluid::surfaces()
{
  return Core::Communication::GetElementSurfaces<FluidBoundary, Fluid>(*this);
}


/*----------------------------------------------------------------------*
 |  get face element (public)                               schott 03/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Fluid::create_face_element(
    Core::Elements::Element* parent_slave,  //!< parent slave fluid3 element
    int nnode,                              //!< number of surface nodes
    const int* nodeids,                     //!< node ids of surface element
    Core::Nodes::Node** nodes,              //!< nodes of surface element
    const int lsurface_master,              //!< local surface number w.r.t master parent element
    const int lsurface_slave,               //!< local surface number w.r.t slave parent element
    const std::vector<int>& localtrafomap   //! local trafo map
)
{
  // dynamic cast for slave parent element
  Discret::ELEMENTS::Fluid* slave_pele = dynamic_cast<Discret::ELEMENTS::Fluid*>(parent_slave);


  // insert both parent elements
  return Core::Communication::ElementIntFaceFactory<FluidIntFace, Fluid>(
      -1,               //!< internal face element id
      -1,               //!< owner of internal face element
      nnode,            //!< number of surface nodes
      nodeids,          //!< node ids of surface element
      nodes,            //!< nodes of surface element
      this,             //!< master parent element
      slave_pele,       //!< slave parent element
      lsurface_master,  //!< local surface number w.r.t master parent element
      lsurface_slave,   //!< local surface number w.r.t slave parent element
      localtrafomap     //!< local trafo map
  );
}


/*----------------------------------------------------------------------*
 |  activate time dependent subgrid scales (public)      gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Fluid::activate_tds(
    int nquad, int nsd, double** saccn, double** sveln, double** svelnp)
{
  if (tds_ == Teuchos::null) tds_ = Teuchos::rcp(new FLD::TDSEleData());

  tds_->activate_tds(nquad, nsd, saccn, sveln, svelnp);
}

FOUR_C_NAMESPACE_CLOSE
