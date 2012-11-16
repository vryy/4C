/*!----------------------------------------------------------------------
\file fluid_ele_intfaces.cpp
\brief

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_discret.H"

#include "fluid_ele.H"
#include "fluid_ele_intfaces_calc.H"

#include <Teuchos_TimeMonitor.hpp>


DRT::ELEMENTS::FluidIntFaceType DRT::ELEMENTS::FluidIntFaceType::instance_;


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidIntFaceType::Create( const int id, const int owner )
{
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                           schott 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFace::FluidIntFace(int id,                                ///< element id
                                          int owner,                             ///< owner (= owner of parent element with smallest gid)
                                          int nnode,                             ///< number of nodes
                                          const int* nodeids,                    ///< node ids
                                          DRT::Node** nodes,                     ///< nodes of surface
                                          DRT::ELEMENTS::Fluid* parent_master,  ///< master parent element
                                          DRT::ELEMENTS::Fluid* parent_slave,   ///< slave parent element
                                          const int lsurface_master              ///< local surface index with respect to master parent element
                                          ) :
DRT::Element(id,owner),
parent_master_(parent_master),
parent_slave_(parent_slave),
lsurface_master_(lsurface_master)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      schott 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFace::FluidIntFace(const DRT::ELEMENTS::FluidIntFace& old) :
DRT::Element(old),
parent_master_(old.parent_master_),
parent_slave_(old.parent_slave_),
lsurface_master_(old.lsurface_master_)
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
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), parent_master_->Shape());
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
DRT::ELEMENTS::FluidIntFace::~FluidIntFace()
{
  return;
}


/*----------------------------------------------------------------------*
 |  create the patch location vector (public)              schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::PatchLocationVector(
    DRT::Discretization & discretization,       ///< discretization
    std::vector<int>&     nds_master,           ///< nodal dofset w.r.t master parent element
    std::vector<int>&     nds_slave,            ///< nodal dofset w.r.t slave parent element
    std::vector<int>&     patchlm,              ///< local map for gdof ids for patch of elements
    std::vector<int>&     master_lm,            ///< local map for gdof ids for master element
    std::vector<int>&     slave_lm,             ///< local map for gdof ids for slave element
    std::vector<int>&     face_lm,              ///< local map for gdof ids for face element
    std::vector<int>&     lm_masterToPatch,     ///< local map between lm_master and lm_patch
    std::vector<int>&     lm_slaveToPatch,      ///< local map between lm_slave and lm_patch
    std::vector<int>&     lm_faceToPatch,       ///< local map between lm_face and lm_patch
    std::vector<int>&     lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&     lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
    )
{
  // create one patch location vector containing all dofs of master, slave and
  // *this FluidIntFace element only once (no duplicates)
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: PatchLocationVector" );

  //-----------------------------------------------------------------------
  const int m_numnode = parent_master_->NumNode();
  DRT::Node** m_nodes = parent_master_->Nodes();

  if ( m_numnode != static_cast<int>( nds_master.size() ) )
  {
    throw std::runtime_error( "wrong number of nodes for master element" );
  }

  //-----------------------------------------------------------------------
  const int s_numnode = parent_slave_->NumNode();
  DRT::Node** s_nodes = parent_slave_->Nodes();

  if ( s_numnode != static_cast<int>( nds_slave.size() ) )
  {
    throw std::runtime_error( "wrong number of nodes for slave element" );
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

  // ---------------------------------------------------
  int dofset = 0; // assume dofset 0

  int patchnode_count = 0;

  // fill patch lm with master's nodes
  for (int k=0; k<m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof = discretization.Dof(dofset,node);

    // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D case, does not return dofset's numnode)
    const int size = NumDofPerNode(dofset,*node);
    const int offset = size*nds_master[k];

#ifdef DEBUG
    if ( dof.size() < static_cast<unsigned>( offset+size ) )
    {
      dserror( "illegal physical dofs offset" );
    }
#endif

    //insert a pair of node-Id and current length of master_lm ( to get the start offset for node's dofs)
    m_node_lm_offset.insert(std::pair<int,int>(node->Id(), master_lm.size()));

    for (int j=0; j< size; ++j)
    {
      int actdof = dof[offset + j];

      // current last index will be the index for next push_back operation
      lm_masterToPatch.push_back( (patchlm.size()) );

      patchlm.push_back(actdof);
      master_lm.push_back(actdof);
    }

    lm_masterNodeToPatch.push_back(patchnode_count);

    patchnode_count++;
  }

  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  for (int k=0; k<s_numnode; ++k)
  {
    DRT::Node* node = s_nodes[k];

    // slave node already contained?
    std::map<int,int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->Id());

    if(m_offset==m_node_lm_offset.end()) // node not included yet
    {
      std::vector<int> dof = discretization.Dof(dofset,node);

      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D case, does not return dofset's numnode)
      const int size = NumDofPerNode(dofset,*node);
      const int offset = size*nds_slave[k];

  #ifdef DEBUG
      if ( dof.size() < static_cast<unsigned>( offset+size ) )
      {
        dserror( "illegal physical dofs offset" );
      }
  #endif
      for (int j=0; j< size; ++j)
      {
        int actdof = dof[offset + j];

        lm_slaveToPatch.push_back( patchlm.size() );

        patchlm.push_back(actdof);
        slave_lm.push_back(actdof);

      }

      lm_slaveNodeToPatch.push_back(patchnode_count);

      patchnode_count++;

    }
    else // node is also a master's node
    {
      const int size = NumDofPerNode(dofset,*node);

      int offset = m_offset->second;

      for (int j=0; j< size; ++j)
      {
        int actdof = master_lm[offset + j];

        slave_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back( lm_masterToPatch[offset + j] );
      }

      if(offset%size != 0) dserror("there was at least one node with not %d dofs per node", size);
      int patchnode_index = offset/size;

      lm_slaveNodeToPatch.push_back(patchnode_index);
      // no patchnode_count++; (node already contained)

    }
  }

  // ---------------------------------------------------
  // extract face's lm from patch_lm
  for (int k=0; k<f_numnode; ++k)
  {
    DRT::Node* node = f_nodes[k];

    // face node must be contained
    std::map<int,int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->Id());

    if(m_offset!=m_node_lm_offset.end()) // node not included yet
    {
      const int size = NumDofPerNode(dofset,*node);

      int offset = m_offset->second;

      for (int j=0; j< size; ++j)
      {
        int actdof = master_lm[offset + j];

        face_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_faceToPatch.push_back( lm_masterToPatch[offset + j] );
      }
    }
    else throw std::runtime_error( "face's nodes not contained in masternodes_offset map" );

  }


  return;
}





/*----------------------------------------------------------------------*
 |  print this element (public)                            schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFace::Print(ostream& os) const
{
  os << "FluidIntFace ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           schott 03/12 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::FluidIntFace::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of FluidIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           schott 03/12 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::FluidIntFace::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of FluidIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                          schott 03/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidIntFace::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // REMARK: this line ensures that the static DRT::ELEMENTS::FluidIntFaceImplInterface::Impl is created
  //         this line avoids linker errors
  DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(this);

  dserror("not available");

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition    schott 03/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidIntFace::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not available");

  return 0;
}
