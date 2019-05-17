/*----------------------------------------------------------------------*/
/*!

\brief Fluid internal faces elements

\maintainer Christoph Ager

\level 2

*/
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_discret.H"

#include "fluid_ele.H"
#include "fluid_ele_intfaces_calc.H"

#include <Teuchos_TimeMonitor.hpp>


DRT::ELEMENTS::FluidIntFaceType DRT::ELEMENTS::FluidIntFaceType::instance_;

DRT::ELEMENTS::FluidIntFaceType& DRT::ELEMENTS::FluidIntFaceType::Instance() { return instance_; }


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidIntFaceType::Create(const int id, const int owner)
{
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                           schott 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFace::FluidIntFace(int id,  ///< element id
    int owner,                            ///< owner (= owner of parent element with smallest gid)
    int nnode,                            ///< number of nodes
    const int* nodeids,                   ///< node ids
    DRT::Node** nodes,                    ///< nodes of surface
    DRT::ELEMENTS::Fluid* parent_master,  ///< master parent element
    DRT::ELEMENTS::Fluid* parent_slave,   ///< slave parent element
    const int lsurface_master,  ///< local surface index with respect to master parent element
    const int lsurface_slave,   ///< local surface index with respect to slave parent element
    const std::vector<int>
        localtrafomap  ///< get the transformation map between the local coordinate systems of the
                       ///< face w.r.t the master parent element's face's coordinate system and the
                       ///< slave element's face's coordinate system
    )
    : DRT::FaceElement(id, owner)
{
  SetParentMasterElement(parent_master, lsurface_master);
  SetParentSlaveElement(parent_slave, lsurface_slave);
  SetLocalTrafoMap(localtrafomap);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      schott 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFace::FluidIntFace(const DRT::ELEMENTS::FluidIntFace& old)
    : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                          schott 03/12|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidIntFace::Clone() const
{
  DRT::ELEMENTS::FluidIntFace* newelement = new DRT::ELEMENTS::FluidIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         schott 03/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::FluidIntFace::Shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::Pack(DRT::PackBuffer& data) const
{
  dserror("this FluidIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::Unpack(const std::vector<char>& data)
{
  dserror("this FluidIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                          schott 03/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFace::~FluidIntFace() { return; }


/*----------------------------------------------------------------------*
 |  create the patch location vector (public)              schott 06/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::PatchLocationVector(
    DRT::Discretization& discretization,     ///< discretization
    std::vector<int>& nds_master,            ///< nodal dofset w.r.t master parent element
    std::vector<int>& nds_slave,             ///< nodal dofset w.r.t slave parent element
    std::vector<int>& patchlm,               ///< local map for gdof ids for patch of elements
    std::vector<int>& lm_masterToPatch,      ///< local map between lm_master and lm_patch
    std::vector<int>& lm_slaveToPatch,       ///< local map between lm_slave and lm_patch
    std::vector<int>& lm_faceToPatch,        ///< local map between lm_face and lm_patch
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    Teuchos::RCP<std::map<int, int>>
        pbcconnectivity  ///< connectivity between slave and PBC's master nodes
)
{
  // create one patch location vector containing all dofs of master, slave and
  // *this FluidIntFace element only once (no duplicates)
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: PatchLocationVector");


  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();

  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    throw std::runtime_error("wrong number of nodes for master element");
  }

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    throw std::runtime_error("wrong number of nodes for slave element");
  }

  //-----------------------------------------------------------------------
  const int f_numnode = NumNode();
  DRT::Node** f_nodes = Nodes();



  // for each master node, the offset for node's dofs in master_lm
  std::map<int, int> m_node_lm_offset;

  //----------------------------------------------------
  int curr_patch_lm_size = 0;  // patch_lm.size() (equal to master_lm.size() during the fill of
                               // patch data with master data)

  //----------------------------------------------------
  // check for PBC nodes
  bool has_PBC = (pbcconnectivity != Teuchos::null);


  // ---------------------------------------------------
  const int dofset = 0;  // assume dofset 0

  int patchnode_count = 0;

  // fill patch lm with master's nodes
  for (int k = 0; k < m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof;
    discretization.Dof(dof, node, dofset, nds_master[k]);

    const int size = dof.size();

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // insert a pair of node-Id and current length of master_lm ( to get the start offset for node's
    // dofs)
    m_node_lm_offset.insert(std::pair<int, int>(nid, curr_patch_lm_size));

    for (int j = 0; j < size; ++j)
    {
      lm_masterToPatch.push_back(curr_patch_lm_size);

      int actdof = dof[j];
      patchlm.push_back(actdof);
      curr_patch_lm_size++;
    }

    lm_masterNodeToPatch.push_back(patchnode_count);

    patchnode_count++;
  }



  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  for (int k = 0; k < s_numnode; ++k)
  {
    DRT::Node* node = s_nodes[k];

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // slave node already contained?
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(nid);

    if (m_offset == m_node_lm_offset.end())  // node not included yet
    {
      std::vector<int> dof;  // = discretization.Dof(dofset,node);
      discretization.Dof(dof, node, dofset,
          nds_slave[k]);  // in case of pbcs, the right dofs are stored also for the slave node

      const int size = dof.size();

      for (int j = 0; j < size; ++j)
      {
        lm_slaveToPatch.push_back(curr_patch_lm_size);

        int actdof = dof[j];
        patchlm.push_back(actdof);
        curr_patch_lm_size++;
      }

      lm_slaveNodeToPatch.push_back(patchnode_count);

      patchnode_count++;
    }
    else  // node is also a master's node
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back(lm_masterToPatch[offset + j]);
      }

      if (offset % size != 0)
        dserror("there was at least one node with not %d dofs per node", size);
      int patchnode_index = offset / size;

      lm_slaveNodeToPatch.push_back(patchnode_index);
      // no patchnode_count++; (node already contained)
    }
  }

  // ---------------------------------------------------
  // extract face's lm from patch_lm

  for (int k = 0; k < f_numnode; ++k)
  {
    DRT::Node* node = f_nodes[k];

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // face node must be contained
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(nid);

    if (m_offset != m_node_lm_offset.end())  // node not included yet
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        // copy from lm_masterToPatch
        lm_faceToPatch.push_back(lm_masterToPatch[offset + j]);
      }
    }
    else
      throw std::runtime_error("face's nodes not contained in masternodes_offset map");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create the patch location vector (public)              schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::PatchLocationVector(
    DRT::Discretization& discretization,     ///< discretization
    std::vector<int>& nds_master,            ///< nodal dofset w.r.t master parent element
    std::vector<int>& nds_slave,             ///< nodal dofset w.r.t slave parent element
    std::vector<int>& patchlm,               ///< local map for gdof ids for patch of elements
    std::vector<int>& master_lm,             ///< local map for gdof ids for master element
    std::vector<int>& slave_lm,              ///< local map for gdof ids for slave element
    std::vector<int>& face_lm,               ///< local map for gdof ids for face element
    std::vector<int>& lm_masterToPatch,      ///< local map between lm_master and lm_patch
    std::vector<int>& lm_slaveToPatch,       ///< local map between lm_slave and lm_patch
    std::vector<int>& lm_faceToPatch,        ///< local map between lm_face and lm_patch
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    Teuchos::RCP<std::map<int, int>>
        pbcconnectivity  ///< connectivity between slave and PBC's master nodes
)
{
  // create one patch location vector containing all dofs of master, slave and
  // *this FluidIntFace element only once (no duplicates)
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: PatchLocationVector");

  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();

  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    throw std::runtime_error("wrong number of nodes for master element");
  }

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    throw std::runtime_error("wrong number of nodes for slave element");
  }

  //-----------------------------------------------------------------------
  const int f_numnode = NumNode();
  DRT::Node** f_nodes = Nodes();

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

  //----------------------------------------------------
  int curr_patch_lm_size = 0;   // patch_lm.size()
  int curr_master_lm_size = 0;  // master_lm.size()

  //----------------------------------------------------
  // check for PBC nodes
  bool has_PBC = (pbcconnectivity != Teuchos::null);

  // ---------------------------------------------------
  const int dofset = 0;  // assume dofset 0

  int patchnode_count = 0;


  // fill patch lm with master's nodes
  for (int k = 0; k < m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof;
    discretization.Dof(dof, node, dofset, nds_master[k]);

    const int size = dof.size();

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // insert a pair of node-Id and current length of master_lm ( to get the start offset for node's
    // dofs)

    m_node_lm_offset.insert(std::pair<int, int>(nid, curr_master_lm_size));

    for (int j = 0; j < size; ++j)
    {
      int actdof = dof[j];

      // current last index will be the index for next push_back operation
      lm_masterToPatch.push_back(curr_patch_lm_size);

      patchlm.push_back(actdof);
      curr_patch_lm_size++;

      master_lm.push_back(actdof);
      curr_master_lm_size++;
    }

    lm_masterNodeToPatch.push_back(patchnode_count);

    patchnode_count++;
  }


  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  for (int k = 0; k < s_numnode; ++k)
  {
    DRT::Node* node = s_nodes[k];

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // slave node already contained?
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(nid);

    if (m_offset == m_node_lm_offset.end())  // node not included yet
    {
      std::vector<int> dof;  // = discretization.Dof(dofset,node);
      discretization.Dof(dof, node, dofset, nds_slave[k]);

      const int size = dof.size();

      for (int j = 0; j < size; ++j)
      {
        int actdof = dof[j];

        lm_slaveToPatch.push_back(curr_patch_lm_size);

        patchlm.push_back(actdof);
        curr_patch_lm_size++;

        slave_lm.push_back(actdof);
      }

      lm_slaveNodeToPatch.push_back(patchnode_count);

      patchnode_count++;
    }
    else  // node is also a master's node
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        int actdof = master_lm[offset + j];

        slave_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back(lm_masterToPatch[offset + j]);
      }

      if (offset % size != 0)
        dserror("there was at least one node with not %d dofs per node", size);
      int patchnode_index = offset / size;

      lm_slaveNodeToPatch.push_back(patchnode_index);
      // no patchnode_count++; (node already contained)
    }
  }

  // ---------------------------------------------------
  // extract face's lm from patch_lm
  for (int k = 0; k < f_numnode; ++k)
  {
    DRT::Node* node = f_nodes[k];

    int nid = node->Id();  // node id of node and id of master node if it is a PBC node

    if (has_PBC)  // set the id of the master node if the node is a PBC node
    {
      std::map<int, int>::iterator slave_it = pbcconnectivity->find(
          nid);  // find the slave node id, is there a corresponding pbc master node?

      if (slave_it != pbcconnectivity->end()) nid = slave_it->second;
    }

    // face node must be contained
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(nid);

    if (m_offset != m_node_lm_offset.end())  // node not included yet
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

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
      throw std::runtime_error("face's nodes not contained in masternodes_offset map");
  }


  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                            schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::Print(std::ostream& os) const
{
  os << "FluidIntFace ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           schott 03/12 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidIntFace::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of FluidIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element>> lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           schott 03/12 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::FluidIntFace::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of FluidIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element>> surfaces(0);
  return surfaces;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                          schott 03/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidIntFace::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // REMARK: this line ensures that the static DRT::ELEMENTS::FluidIntFaceImplInterface::Impl is
  // created
  //         this line avoids linker errors
  DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(this);

  dserror("not available");

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition    schott 03/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidIntFace::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not available");

  return 0;
}
