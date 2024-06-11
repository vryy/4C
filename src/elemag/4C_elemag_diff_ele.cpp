/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of electromagnetic diffusion elements

<pre>
\level 2

</pre>
*/
/*----------------------------------------------------------------------*/

#include "4C_elemag_diff_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_elemag_ele_boundary_calc.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::ElemagDiffType Discret::ELEMENTS::ElemagDiffType::instance_;
Discret::ELEMENTS::ElemagDiffBoundaryType Discret::ELEMENTS::ElemagDiffBoundaryType::instance_;
Discret::ELEMENTS::ElemagDiffIntFaceType Discret::ELEMENTS::ElemagDiffIntFaceType::instance_;

Discret::ELEMENTS::ElemagDiffType& Discret::ELEMENTS::ElemagDiffType::Instance()
{
  return instance_;
}

Discret::ELEMENTS::ElemagDiffBoundaryType& Discret::ELEMENTS::ElemagDiffBoundaryType::Instance()
{
  return instance_;
}

Discret::ELEMENTS::ElemagDiffIntFaceType& Discret::ELEMENTS::ElemagDiffIntFaceType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::ElemagDiffType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::ElemagDiff* object = new Discret::ELEMENTS::ElemagDiff(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagDiffType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "ELECTROMAGNETICDIFF")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::ElemagDiff(id, owner));
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagDiffType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::ElemagDiff(id, owner));
}

void Discret::ELEMENTS::ElemagDiffType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = Core::FE::getDimension(dwele->Shape()) - 1;  // 2;  // Bad Luca! Hard coding is not nice!
  dimns = numdf;
  nv = numdf;
  np = 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::ElemagDiffType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented right now!");
  return nullspace;
}

/*----------------------------------------------------------------------*
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["ELECTROMAGNETICDIFF"];

  // 3D elements
  defs["HEX8"] = Input::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedInt("DEG")
                     .AddNamedInt("SPC")
                     .Build();

  defs["TET4"] = Input::LineDefinition::Builder()
                     .AddIntVector("TET4", 4)
                     .AddNamedInt("MAT")
                     .AddNamedInt("DEG")
                     .AddNamedInt("SPC")
                     .Build();

  // 2D elements
  defs["QUAD4"] = Input::LineDefinition::Builder()
                      .AddIntVector("QUAD4", 4)
                      .AddNamedInt("MAT")
                      .AddNamedInt("DEG")
                      .AddNamedInt("SPC")
                      .Build();

  defs["QUAD9"] = Input::LineDefinition::Builder()
                      .AddIntVector("QUAD9", 9)
                      .AddNamedInt("MAT")
                      .AddNamedInt("DEG")
                      .AddNamedInt("SPC")
                      .Build();

  defs["TRI3"] = Input::LineDefinition::Builder()
                     .AddIntVector("TRI3", 3)
                     .AddNamedInt("MAT")
                     .AddNamedInt("DEG")
                     .AddNamedInt("SPC")
                     .Build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                       berardocco 03/19|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiff::ElemagDiff(int id, int owner) : Elemag(id, owner)
{
  distype_ = Core::FE::CellType::dis_none;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    berardocco 03/19|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiff::ElemagDiff(const Discret::ELEMENTS::ElemagDiff& old) : Elemag(old) {}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Elemag and return pointer to it (public)   |
 |                                                        berardocco 03/19|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ElemagDiff::Clone() const
{
  Discret::ELEMENTS::ElemagDiff* newelement = new Discret::ELEMENTS::ElemagDiff(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                         berardocco 03/19|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiff::Print(std::ostream& os) const
{
  os << "ElemagDiff ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)           berardocco 03/19|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagDiff::Lines()
{
  return Core::Communication::GetElementLines<ElemagDiffBoundary, ElemagDiff>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                     berardocco 03/19|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::ElemagDiff::Surfaces()
{
  return Core::Communication::GetElementSurfaces<ElemagDiffBoundary, ElemagDiff>(*this);
}


/*----------------------------------------------------------------------*
 |  get face element (public)                           berardocco 03/19|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagDiff::CreateFaceElement(
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
  Discret::ELEMENTS::ElemagDiff* slave_pele =
      dynamic_cast<Discret::ELEMENTS::ElemagDiff*>(parent_slave);

  // insert both parent elements
  return Core::Communication::ElementIntFaceFactory<ElemagDiffIntFace, ElemagDiff>(
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


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagDiffBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                      berardocco 03/19 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiffBoundary::ElemagDiffBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::ELEMENTS::ElemagDiff* parent,
    const int lsurface)
    : ElemagBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
  //  set_parent_master_element(parent,lsurface);
  //  SetNodeIds(nnode,nodeids);
  //  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                 berardocco 03/19 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiffBoundary::ElemagDiffBoundary(
    const Discret::ELEMENTS::ElemagDiffBoundary& old)
    : ElemagBoundary(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                     berardocco 03/19 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ElemagDiffBoundary::Clone() const
{
  Discret::ELEMENTS::ElemagDiffBoundary* newelement =
      new Discret::ELEMENTS::ElemagDiffBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffBoundary::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);

  // Discretisation type
  // add_to_pack(data,distype_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffBoundary::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);

  // distype
  // distype_ = static_cast<Core::FE::CellType>( ExtractInt(position,data) );

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                        berardocco 03/19 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffBoundary::Print(std::ostream& os) const
{
  os << "ElemagDiffBoundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                      berardocco 03/19 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::ElemagDiffBoundary::Evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  Discret::ELEMENTS::ElemagBoundaryImplInterface::Impl(this)->Evaluate(
      this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
  return 0;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) berardocco 03/19 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffBoundary::LocationVector(const Core::FE::Discretization& dis,
    LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  // we have to do it this way, just as for weak Dirichlet conditions
  ParentMasterElement()->LocationVector(dis, la, false);
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

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::ElemagDiffIntFaceType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                       berardocco 03/19|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiffIntFace::ElemagDiffIntFace(int id,  ///< element id
    int owner,                  ///< owner (= owner of parent element with smallest gid)
    int nnode,                  ///< number of nodes
    const int* nodeids,         ///< node ids
    Core::Nodes::Node** nodes,  ///< nodes of surface
    Discret::ELEMENTS::ElemagDiff* parent_master,  ///< master parent element
    Discret::ELEMENTS::ElemagDiff* parent_slave,   ///< slave parent element
    const int lsurface_master,  ///< local surface index with respect to master parent element
    const int lsurface_slave,   ///< local surface index with respect to slave parent element
    const std::vector<int>
        localtrafomap  ///< get the transformation map between the local coordinate systems of the
                       ///< face w.r.t the master parent element's face's coordinate system and the
                       ///< slave element's face's coordinate system
    )
    : ElemagIntFace(id, owner, nnode, nodeids, nodes, parent_master, parent_slave, lsurface_master,
          lsurface_slave, localtrafomap)
{
  set_parent_master_element(parent_master, lsurface_master);
  set_parent_slave_element(parent_slave, lsurface_slave);
  set_local_trafo_map(localtrafomap);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                  berardocco 03/19|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagDiffIntFace::ElemagDiffIntFace(
    const Discret::ELEMENTS::ElemagDiffIntFace& old)
    : ElemagIntFace(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                      berardocco 03/19|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::ElemagDiffIntFace::Clone() const
{
  Discret::ELEMENTS::ElemagDiffIntFace* newelement =
      new Discret::ELEMENTS::ElemagDiffIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  create the patch location vector (public)          berardocco 03/19 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffIntFace::PatchLocationVector(
    Core::FE::Discretization& discretization,  ///< discretization
    std::vector<int>& nds_master,              ///< nodal dofset w.r.t master parent element
    std::vector<int>& nds_slave,               ///< nodal dofset w.r.t slave parent element
    std::vector<int>& patchlm,                 ///< local map for gdof ids for patch of elements
    std::vector<int>& master_lm,               ///< local map for gdof ids for master element
    std::vector<int>& slave_lm,                ///< local map for gdof ids for slave element
    std::vector<int>& face_lm,                 ///< local map for gdof ids for face element
    std::vector<int>& lm_masterToPatch,        ///< local map between lm_master and lm_patch
    std::vector<int>& lm_slaveToPatch,         ///< local map between lm_slave and lm_patch
    std::vector<int>& lm_faceToPatch,          ///< local map between lm_face and lm_patch
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch    ///< local map between slave nodes and nodes in patch
)
{
  // create one patch location vector containing all dofs of master, slave and
  // *this ElemagDiffIntFace element only once (no duplicates)

  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->num_node();
  Core::Nodes::Node** m_nodes = ParentMasterElement()->Nodes();

  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    FOUR_C_THROW("wrong number of nodes for master element");
  }

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->num_node();
  Core::Nodes::Node** s_nodes = ParentSlaveElement()->Nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    FOUR_C_THROW("wrong number of nodes for slave element");
  }

  //-----------------------------------------------------------------------
  const int f_numnode = num_node();
  Core::Nodes::Node** f_nodes = Nodes();

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
    std::vector<int> dof = discretization.Dof(dofset, node);

    // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D
    // case, does not return dofset's numnode)
    const int size = discretization.NumDof(dofset, node);
    const int offset = size * nds_master[k];

    FOUR_C_ASSERT(
        dof.size() >= static_cast<unsigned>(offset + size), "illegal physical dofs offset");

    // insert a pair of node-Id and current length of master_lm ( to get the start offset for node's
    // dofs)
    m_node_lm_offset.insert(std::pair<int, int>(node->Id(), master_lm.size()));

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
    m_offset = m_node_lm_offset.find(node->Id());

    if (m_offset == m_node_lm_offset.end())  // node not included yet
    {
      std::vector<int> dof = discretization.Dof(dofset, node);

      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numnode)
      const int size = discretization.NumDof(dofset, node);
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
      const int size = discretization.NumDof(dofset, node);

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
    m_offset = m_node_lm_offset.find(node->Id());

    if (m_offset != m_node_lm_offset.end())  // node not included yet
    {
      const int size = discretization.NumDof(dofset, node);

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
 |  print this element (public)                        berardocco 03/19 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ElemagDiffIntFace::Print(std::ostream& os) const
{
  os << "ElemagDiffIntFace ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
