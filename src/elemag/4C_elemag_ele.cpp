/*----------------------------------------------------------------------*/
/*! \file
\brief Electromagnetic element implementations


\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_elemag_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_elemag_ele_boundary_calc.hpp"
#include "4C_elemag_ele_intfaces_calc.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_so3_nullspace.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::ElemagType Discret::ELEMENTS::ElemagType::instance_;
Discret::ELEMENTS::ElemagBoundaryType Discret::ELEMENTS::ElemagBoundaryType::instance_;
Discret::ELEMENTS::ElemagIntFaceType Discret::ELEMENTS::ElemagIntFaceType::instance_;

Discret::ELEMENTS::ElemagType& Discret::ELEMENTS::ElemagType::instance() { return instance_; }

Discret::ELEMENTS::ElemagBoundaryType& Discret::ELEMENTS::ElemagBoundaryType::instance()
{
  return instance_;
}

Discret::ELEMENTS::ElemagIntFaceType& Discret::ELEMENTS::ElemagIntFaceType::instance()
{
  return instance_;
}


Core::Communication::ParObject* Discret::ELEMENTS::ElemagType::create(const std::vector<char>& data)
{
  Discret::ELEMENTS::Elemag* object = new Discret::ELEMENTS::Elemag(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ELECTROMAGNETIC")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::Elemag(id, owner));
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagType::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Elemag(id, owner));
}


void Discret::ELEMENTS::ElemagType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  nv = Core::FE::getDimension(dwele->shape()) - 1;
  dimns = nv;
  numdf = dimns;
  return;
}


Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::ElemagType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented right now!");
  return nullspace;
}


void Discret::ELEMENTS::ElemagType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["ELECTROMAGNETIC"];

  // 3D elements
  defs["HEX8"] = Input::LineDefinition::Builder()
                     .add_int_vector("HEX8", 8)
                     .add_named_int("MAT")
                     .add_named_int("DEG")
                     .add_named_int("SPC")
                     .build();

  defs["TET4"] = Input::LineDefinition::Builder()
                     .add_int_vector("TET4", 4)
                     .add_named_int("MAT")
                     .add_named_int("DEG")
                     .add_named_int("SPC")
                     .build();

  // 2D elements
  defs["QUAD4"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD4", 4)
                      .add_named_int("MAT")
                      .add_named_int("DEG")
                      .add_named_int("SPC")
                      .build();

  defs["QUAD9"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD9", 9)
                      .add_named_int("MAT")
                      .add_named_int("DEG")
                      .add_named_int("SPC")
                      .build();

  defs["TRI3"] = Input::LineDefinition::Builder()
                     .add_int_vector("TRI3", 3)
                     .add_named_int("MAT")
                     .add_named_int("DEG")
                     .add_named_int("SPC")
                     .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                       berardocco 02/18|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Elemag::Elemag(int id, int owner)
    : Core::Elements::Element(id, owner), degree_(1), completepol_(true)
{
  distype_ = Core::FE::CellType::dis_none;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                  berardocco 02/18|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Elemag::Elemag(const Discret::ELEMENTS::Elemag& old)
    : Core::Elements::Element(old),
      distype_(old.distype_),
      degree_(old.degree_),
      completepol_(old.completepol_)
{
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Elemag and return pointer to it (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Elemag::clone() const
{
  Discret::ELEMENTS::Elemag* newelement = new Discret::ELEMENTS::Elemag(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Elemag::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);

  // Discretisation type
  add_to_pack(data, distype_);
  int degree = degree_;
  add_to_pack(data, degree);
  add_to_pack(data, completepol_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Elemag::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

  // distype
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));
  int val = 0;
  extract_from_pack(position, data, val);
  FOUR_C_ASSERT(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  extract_from_pack(position, data, val);
  completepol_ = val;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                         berardocco 02/18|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Elemag::print(std::ostream& os) const
{
  os << "Elemag ";
  Element::print(os);
}


bool Discret::ELEMENTS::Elemag::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));
  degree_ = container.get<int>("DEG");

  completepol_ = container.get<int>("SPC");

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  set_dis_type(Core::FE::StringToCellType(distype));

  return true;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)           berardocco 02/18|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Elemag::lines()
{
  return Core::Communication::GetElementLines<ElemagBoundary, Elemag>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                     berardocco 02/18|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Elemag::surfaces()
{
  return Core::Communication::GetElementSurfaces<ElemagBoundary, Elemag>(*this);
}

/*----------------------------------------------------------------------*
 |  get face element (public)                           berardocco 02/18|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Elemag::create_face_element(
    Core::Elements::Element* parent_slave, int nnode, const int* nodeids, Core::Nodes::Node** nodes,
    const int lsurface_master, const int lsurface_slave, const std::vector<int>& localtrafomap)
{
  // dynamic cast for slave parent element
  Discret::ELEMENTS::Elemag* slave_pele = dynamic_cast<Discret::ELEMENTS::Elemag*>(parent_slave);

  // insert both parent elements
  return Core::Communication::ElementIntFaceFactory<ElemagIntFace, Elemag>(-1, -1, nnode, nodeids,
      nodes, this, slave_pele, lsurface_master, lsurface_slave, localtrafomap);
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagBoundaryType::create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                      berardocco 02/18 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagBoundary::ElemagBoundary(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Elemag* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  set_parent_master_element(parent, lsurface);
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                 berardocco 02/18 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagBoundary::ElemagBoundary(const Discret::ELEMENTS::ElemagBoundary& old)
    : Core::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                     berardocco 02/18 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ElemagBoundary::clone() const
{
  Discret::ELEMENTS::ElemagBoundary* newelement = new Discret::ELEMENTS::ElemagBoundary(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::ElemagBoundary::shape() const
{
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_master_element()->shape());
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagBoundary::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);

  // Discretisation type
  // add_to_pack(data,distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagBoundary::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

  // distype
  // distype_ = static_cast<Core::FE::CellType>( extract_int(position,data) );

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                        berardocco 02/18 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagBoundary::print(std::ostream& os) const
{
  os << "ElemagBoundary ";
  Element::print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       berardocco 02/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagBoundary::lines()
{
  FOUR_C_THROW("Lines of ElemagBoundary not implemented");
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       berardocco 02/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagBoundary::surfaces()
{
  FOUR_C_THROW("Surfaces of ElemagBoundary not implemented");
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                      berardocco 02/18 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::ElemagBoundary::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  Discret::ELEMENTS::ElemagBoundaryImplInterface::impl(this)->evaluate(
      this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
  return 0;
}


/*-----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition berardocco 02/18 |
 *-----------------------------------------------------------------------*/
int Discret::ELEMENTS::ElemagBoundary::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("dummy function called");
  return 0;
}

/*------------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) berardocco 02/18 |
 *------------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagBoundary::location_vector(const Core::FE::Discretization& dis,
    LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  // we have to do it this way, just as for weak Dirichlet conditions
  parent_master_element()->location_vector(dis, la, false);
  return;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagIntFaceType::create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                       berardocco 02/18|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagIntFace::ElemagIntFace(int id,  // element id
    int owner,                  // owner (= owner of parent element with smallest gid)
    int nnode,                  // number of nodes
    const int* nodeids,         // node ids
    Core::Nodes::Node** nodes,  // nodes of surface
    Discret::ELEMENTS::Elemag* parent_master,  // master parent element
    Discret::ELEMENTS::Elemag* parent_slave,   // slave parent element
    const int lsurface_master,  // local surface index with respect to master parent element
    const int lsurface_slave,   // local surface index with respect to slave parent element
    const std::vector<int>
        localtrafomap  // get the transformation map between the local coordinate systems of the
                       // face w.r.t the master parent element's face's coordinate system and the
                       // slave element's face's coordinate system
    )
    : Core::Elements::FaceElement(id, owner)
{
  set_parent_master_element(parent_master, lsurface_master);
  set_parent_slave_element(parent_slave, lsurface_slave);

  if (parent_slave != nullptr)
    degree_ = std::max(parent_master->degree(), parent_slave->degree());
  else
    degree_ = parent_master->degree();

  set_local_trafo_map(localtrafomap);

  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                  berardocco 02/18|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagIntFace::ElemagIntFace(const Discret::ELEMENTS::ElemagIntFace& old)
    : Core::Elements::FaceElement(old), degree_(old.degree_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                      berardocco 02/18|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ElemagIntFace::clone() const
{
  Discret::ELEMENTS::ElemagIntFace* newelement = new Discret::ELEMENTS::ElemagIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                     berardocco 02/18 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::ElemagIntFace::shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_master_element()->shape());
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                     berardocco 02/18 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagIntFace::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this ElemagIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                     berardocco 02/18 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagIntFace::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this ElemagIntFace element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  create the patch location vector (public)          berardocco 02/18 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagIntFace::patch_location_vector(
    Core::FE::Discretization& discretization,  // discretization
    std::vector<int>& nds_master,              // nodal dofset w.r.t master parent element
    std::vector<int>& nds_slave,               // nodal dofset w.r.t slave parent element
    std::vector<int>& patchlm,                 // local map for gdof ids for patch of elements
    std::vector<int>& master_lm,               // local map for gdof ids for master element
    std::vector<int>& slave_lm,                // local map for gdof ids for slave element
    std::vector<int>& face_lm,                 // local map for gdof ids for face element
    std::vector<int>& lm_masterToPatch,        // local map between lm_master and lm_patch
    std::vector<int>& lm_slaveToPatch,         // local map between lm_slave and lm_patch
    std::vector<int>& lm_faceToPatch,          // local map between lm_face and lm_patch
    std::vector<int>& lm_masterNodeToPatch,    // local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch      // local map between slave nodes and nodes in patch
)
{
  // create one patch location vector containing all dofs of master, slave and
  // *this ElemagIntFace element only once (no duplicates)

  //-----------------------------------------------------------------------
  const int m_numnode = parent_master_element()->num_node();
  Core::Nodes::Node** m_nodes = parent_master_element()->nodes();

  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    FOUR_C_THROW("wrong number of nodes for master element");
  }

  //-----------------------------------------------------------------------
  const int s_numnode = parent_slave_element()->num_node();
  Core::Nodes::Node** s_nodes = parent_slave_element()->nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    FOUR_C_THROW("wrong number of nodes for slave element");
  }

  //-----------------------------------------------------------------------
  const int f_numnode = num_node();
  Core::Nodes::Node** f_nodes = nodes();

  //-----------------------------------------------------------------------
  // create the patch local map and additional local maps between elements lm and patch lm

  patchlm.clear();

  master_lm.clear();
  slave_lm.clear();
  face_lm.clear();

  lm_masterToPatch.clear();
  lm_slaveToPatch.clear();
  lm_faceToPatch.clear();

  // maps between master/slave nodes and nodes in patch
  lm_masterNodeToPatch.clear();
  lm_slaveNodeToPatch.clear();

  // for each master node, the offset for node's dofs in master_lm
  std::map<int, int> m_node_lm_offset;


  // ---------------------------------------------------
  int dofset = 0;  // assume dofset 0

  int patchnode_count = 0;

  // fill patch lm with master's nodes
  for (int k = 0; k < m_numnode; ++k)
  {
    Core::Nodes::Node* node = m_nodes[k];
    std::vector<int> dof = discretization.dof(dofset, node);

    // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D
    // case, does not return dofset's numnode)
    const int size = discretization.num_dof(dofset, node);
    const int offset = size * nds_master[k];

    FOUR_C_ASSERT(
        dof.size() >= static_cast<unsigned>(offset + size), "illegal physical dofs offset");

    // insert a pair of node-Id and current length of master_lm ( to get the start offset for node's
    // dofs)
    m_node_lm_offset.insert(std::pair<int, int>(node->id(), master_lm.size()));

    for (int j = 0; j < size; ++j)
    {
      int actdof = dof[offset + j];

      // current last index will be the index for next push_back operation
      lm_masterToPatch.push_back((patchlm.size()));

      patchlm.push_back(actdof);
      master_lm.push_back(actdof);
    }

    lm_masterNodeToPatch.push_back(patchnode_count);

    patchnode_count++;
  }

  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  for (int k = 0; k < s_numnode; ++k)
  {
    Core::Nodes::Node* node = s_nodes[k];

    // slave node already contained?
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->id());

    if (m_offset == m_node_lm_offset.end())  // node not included yet
    {
      std::vector<int> dof = discretization.dof(dofset, node);

      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numnode)
      const int size = discretization.num_dof(dofset, node);
      const int offset = size * nds_slave[k];

      FOUR_C_ASSERT(
          dof.size() >= static_cast<unsigned>(offset + size), "illegal physical dofs offset");
      for (int j = 0; j < size; ++j)
      {
        int actdof = dof[offset + j];

        lm_slaveToPatch.push_back(patchlm.size());

        patchlm.push_back(actdof);
        slave_lm.push_back(actdof);
      }

      lm_slaveNodeToPatch.push_back(patchnode_count);

      patchnode_count++;
    }
    else  // node is also a master's node
    {
      const int size = discretization.num_dof(dofset, node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        int actdof = master_lm[offset + j];

        slave_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back(lm_masterToPatch[offset + j]);
      }

      if (offset % size != 0)
        FOUR_C_THROW("there was at least one node with not %d dofs per node", size);
      int patchnode_index = offset / size;

      lm_slaveNodeToPatch.push_back(patchnode_index);
      // no patchnode_count++; (node already contained)
    }
  }

  // ---------------------------------------------------
  // extract face's lm from patch_lm
  for (int k = 0; k < f_numnode; ++k)
  {
    Core::Nodes::Node* node = f_nodes[k];

    // face node must be contained
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->id());

    if (m_offset != m_node_lm_offset.end())  // node not included yet
    {
      const int size = discretization.num_dof(dofset, node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        int actdof = master_lm[offset + j];

        face_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_faceToPatch.push_back(lm_masterToPatch[offset + j]);
      }
    }
    else
      FOUR_C_THROW("face's nodes not contained in masternodes_offset map");
  }

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                        berardocco 02/18 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagIntFace::print(std::ostream& os) const
{
  os << "ElemagIntFace ";
  Element::print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       berardocco 02/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagIntFace::lines()
{
  FOUR_C_THROW("Lines of ElemagIntFace not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                       berardocco 02/18 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagIntFace::surfaces()
{
  FOUR_C_THROW("Surfaces of ElemagIntFace not implemented");
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                      berardocco 02/18 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::ElemagIntFace::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // REMARK: this line ensures that the static Discret::ELEMENTS::ElemagIntFaceImplInterface::Impl
  // is created
  //         this line avoids linker errors
  // Discret::ELEMENTS::ElemagIntFaceImplInterface::Impl(this);

  FOUR_C_THROW("not available");

  return 0;
}


/*------------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition  berardocco 02/18 |
 *------------------------------------------------------------------------*/
int Discret::ELEMENTS::ElemagIntFace::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("not available");

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
