/*!----------------------------------------------------------------------
\file combust3_elementface.cpp
\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "combust3.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_xfem/dof_management_element.H"
#include "../drt_xfem/field_enriched.H"


DRT::ELEMENTS::Combust3IntFaceType DRT::ELEMENTS::Combust3IntFaceType::instance_;

DRT::ELEMENTS::Combust3IntFaceType& DRT::ELEMENTS::Combust3IntFaceType::Instance()
{
  return instance_;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Combust3IntFaceType::Create( const int id, const int owner )
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                        rasthofer 02/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3IntFace::Combust3IntFace(int id,                                ///< element id
                                          int owner,                             ///< owner (= owner of parent element with smallest gid)
                                          int nnode,                             ///< number of nodes
                                          const int* nodeids,                    ///< node ids
                                          DRT::Node** nodes,                     ///< nodes of surface
                                          DRT::ELEMENTS::Combust3* parent_master,   ///< master parent element
                                          DRT::ELEMENTS::Combust3* parent_slave,    ///< slave parent element
                                          const int lsurface_master,             ///< local surface index with respect to master parent element
                                          const int lsurface_slave,              ///< local surface index with respect to slave parent element
                                          const std::vector<int> localtrafomap   ///< get the transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system and the slave element's face's coordinate system
):
DRT::FaceElement(id,owner)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent_master,lsurface_master);
  SetParentSlaveElement(parent_slave,lsurface_slave);
  SetLocalTrafoMap(localtrafomap);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   rasthofer 02/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3IntFace::Combust3IntFace(const DRT::ELEMENTS::Combust3IntFace& old) :
DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                       rasthofer 02/13|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Combust3IntFace::Clone() const
{
  DRT::ELEMENTS::Combust3IntFace* newelement = new DRT::ELEMENTS::Combust3IntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                      rasthofer 02/13 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Combust3IntFace::Shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      rasthofer 02/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3IntFace::Pack(DRT::PackBuffer& data) const
{
  dserror("this Combust3IntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         rasthofer 02/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3IntFace::Unpack(const std::vector<char>& data)
{
  dserror("this Combust3IntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                       rasthofer 02/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3IntFace::~Combust3IntFace()
{
  return;
}


/*----------------------------------------------------------------------*
 |  create the patch location vector (public)           rasthofer 02/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3IntFace::PatchLocationVector(
    DRT::Discretization & discretization,       ///< discretization
    std::vector<int>&     master_lm,            ///< local map for gdof ids for master element
    std::vector<int>&     slave_lm,             ///< local map for gdof ids for slave element
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&  patch_components_lm,     ///< rearranged local map for gdof ids for patch of elements
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&  patch_components_lmowner ///< owner of patch
    )
{
  // create one patch location vector containing all dofs of master, slave and
  // *this Combust3IntFace element only once (no duplicates)

  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

  //-----------------------------------------------------------------------
  const int f_numnode = NumNode();
  DRT::Node** f_nodes = Nodes();

  // modify the patch owner to the owner of the internal face element
  int owner = this->Owner();

  // get element dof manager
  Teuchos::RCP<XFEM::ElementDofManager> master_eledofmanager = ParentMasterElement()->GetEleDofManager();
  Teuchos::RCP<XFEM::ElementDofManager> slave_eledofmanager = ParentSlaveElement()->GetEleDofManager();
//      std::cout << "Master num dof node  " << master_eledofmanager->NumNodeDof() << " Slave num dof node  " << slave_eledofmanager->NumNodeDof() << std::endl;
//      std::cout << "Master dof field  " << master_eledofmanager->NumDofPerField(XFEM::PHYSICS::Velx) << "Slave dof field  " << slave_eledofmanager->NumDofPerField(XFEM::PHYSICS::Velx) << std::endl;


  const int numfieldpernode = master_eledofmanager->NumFields();
  if (numfieldpernode != slave_eledofmanager->NumFields())
    dserror("Same number of fields expected");

  //-----------------------------------------------------------------------
  // create the patch local map and additional local maps between elements lm and patch lm

  // local maps for patch dofs
  std::vector<int> patchlm;

  // local maps for face dofs
  std::vector<int> face_lm;

  // local maps between master/slave dofs and position in patch dofs (lm_patch)
  std::vector<int> lm_masterToPatch;
  std::vector<int> lm_slaveToPatch;
  std::vector<int> lm_faceToPatch;

  patchlm.clear();

  master_lm.clear();
  slave_lm.clear();
  face_lm.clear();

  lm_masterToPatch.clear();
  lm_slaveToPatch.clear();
  lm_faceToPatch.clear();

  // for each master node, the offset for node's dofs in master_lm
  std::map<int, int> m_node_lm_offset;
  std::map<XFEM::PHYSICS::Field,std::map<int, int> > m_dof_per_field_offset;
  std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator iterator;
  for (iterator = lm_masterDofPerFieldToPatch.begin(); iterator != lm_masterDofPerFieldToPatch.end(); iterator++)
      m_dof_per_field_offset[iterator->first]=std::map<int, int>();

  // ---------------------------------------------------
  int dofset = 0; // assume dofset 0

  // patch node counter
  std::map<XFEM::PHYSICS::Field,int > fieldpatchdof_count;
  std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator iter;
  for (iter = lm_masterDofPerFieldToPatch.begin(); iter != lm_masterDofPerFieldToPatch.end(); iter++)
    fieldpatchdof_count[iter->first] = 0;

  int masterdof_count = 0;

  // fill patch lm with master's nodes
  for (int k=0; k<m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof = discretization.Dof(dofset,node);

    // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D case, does not return dofset's numnode)
    const std::size_t size = master_eledofmanager->NumDofPerNode(node->Id());

    //insert a pair of node-Id and current length of master_lm ( to get the start offset for node's dofs)
    m_node_lm_offset.insert(std::pair<int,int>(node->Id(), master_lm.size()));

    std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator it;
    for (it = lm_masterDofPerFieldToPatch.begin(); it != lm_masterDofPerFieldToPatch.end(); it++)
    {
        m_dof_per_field_offset[it->first].insert(std::pair<int,int>(node->Id(), lm_masterDofPerFieldToPatch[it->first].size()));
    }

//    std::cout << "NodeID  " << node->Id() << std::endl;

    for (int j=0; j< ((int)size); ++j)
    {
      int actdof = dof[j];

      // current last index will be the index for next push_back operation
      lm_masterToPatch.push_back( (patchlm.size()) );

      patchlm.push_back(actdof);
      master_lm.push_back(actdof);

      const XFEM::FieldEnr fieldenr = master_eledofmanager->FieldEnrSetPerDof(masterdof_count);
//      std::cout << "Field  " << physVarToString(fieldenr.getField()) << " dof  " << fieldenr.getEnrichment().toString() << std::endl;

      patch_components_lm[fieldenr.getField()].push_back(actdof);
      patch_components_lmowner[fieldenr.getField()].push_back(owner);

      lm_masterDofPerFieldToPatch[fieldenr.getField()].push_back(fieldpatchdof_count[fieldenr.getField()]);
//      std::cout << " dof per field " << fieldpatchdof_count[fieldenr.getField()] << std::endl;
      fieldpatchdof_count[fieldenr.getField()]++;

      masterdof_count++;
    }
  }

  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  int slavedof_count = 0;

  for (int k=0; k<s_numnode; ++k)
  {
    DRT::Node* node = s_nodes[k];

//    std::cout << "NodeID  " << node->Id() << std::endl;

    // slave node already contained?
    std::map<int,int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->Id());

    if(m_offset==m_node_lm_offset.end()) // node not included yet
    {
      std::vector<int> dof = discretization.Dof(dofset,node);

      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D case, does not return dofset's numnode)
      const int size = slave_eledofmanager->NumDofPerNode(node->Id());

      for (int j=0; j< size; ++j)
      {
          int actdof = dof[j];

        lm_slaveToPatch.push_back( patchlm.size() );

        patchlm.push_back(actdof);
        slave_lm.push_back(actdof);

        const XFEM::FieldEnr fieldenr = slave_eledofmanager->FieldEnrSetPerDof(slavedof_count);

//        std::cout << "slave new dof" << std::endl;
//        std::cout << "Field  " << physVarToString(fieldenr.getField()) << " dof  " << fieldenr.getEnrichment().toString() << std::endl;

        patch_components_lm[fieldenr.getField()].push_back(actdof);
        patch_components_lmowner[fieldenr.getField()].push_back(owner);

        lm_slaveDofPerFieldToPatch[fieldenr.getField()].push_back(fieldpatchdof_count[fieldenr.getField()]);
        fieldpatchdof_count[fieldenr.getField()]++;

        slavedof_count++;
      }

    }
    else // node is also a master's node
    {
      const std::size_t size = slave_eledofmanager->NumDofPerNode(node->Id());

      int offset = m_offset->second;

      std::map<XFEM::PHYSICS::Field, int> offset_per_field;
      std::map<XFEM::PHYSICS::Field,std::map<int, int> >::iterator it = m_dof_per_field_offset.begin();
      for (it = m_dof_per_field_offset.begin(); it != m_dof_per_field_offset.end(); it++)
        offset_per_field[it->first] = (it->second)[node->Id()];

      for (int j=0; j< ((int)size); ++j)
      {
        int actdof = master_lm[offset + j];

        slave_lm.push_back(actdof);

        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back( lm_masterToPatch[offset + j] );

        const XFEM::FieldEnr fieldenr = slave_eledofmanager->FieldEnrSetPerDof(slavedof_count);

//        std::cout << "slave old dof  "<< std::endl;
//        std::cout << "Field  " << physVarToString(fieldenr.getField()) << " offset "<< offset_per_field[fieldenr.getField()] << std::endl;


        lm_slaveDofPerFieldToPatch[fieldenr.getField()].push_back(offset_per_field[fieldenr.getField()]);
        offset_per_field[fieldenr.getField()]++;

        slavedof_count++;
      }

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
      const int size = master_eledofmanager->NumDofPerNode(node->Id());

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

// some debug output
//      const std::vector<XFEM::PHYSICS::Field> fields = master_eledofmanager->GetFields();
//      std::cout << "lm_patch  " << patchlm.size() << std::endl;
//      std::cout << "lm_master  " << master_lm.size() << std::endl;
//      std::cout << "lm_slave  " << slave_lm.size() << std::endl;
//      std::cout << "lm_face  " << face_lm.size() << std::endl;
//      std::cout << "lm_masterToPatch  " << lm_masterToPatch.size() << std::endl;
//      std::cout << "lm_slaveToPatch  " << lm_slaveToPatch.size() << std::endl;
//      std::cout << "lm_faceToPatch  " << lm_faceToPatch.size() << std::endl;
//      std::cout << "lm_masterDofPerFieldToPatch  " << lm_masterDofPerFieldToPatch.size() << std::endl;
//      std::cout << "lm_slaveDofPerFieldToPatch  " << lm_slaveDofPerFieldToPatch.size() << std::endl;
//      for (std::size_t i=0; i< lm_masterDofPerFieldToPatch.size(); i++)
//      {
//          std::cout << "master  " << std::endl;
//         for (std::size_t j=0; j<lm_masterDofPerFieldToPatch[fields[i]].size(); j++)
//         {
//           std::cout << (lm_masterDofPerFieldToPatch[fields[i]])[j] << std::endl;
//         }
//         std::cout << "salve  " << std::endl;
//         for (std::size_t j=0; j<lm_slaveDofPerFieldToPatch[fields[i]].size(); j++)
//         {
//           std::cout << (lm_slaveDofPerFieldToPatch[fields[i]])[j] << std::endl;
//         }
//      }
//
//     std::cout << "patch_components_lm  " << patch_components_lm.size() << std::endl;
//     std::map<XFEM::PHYSICS::Field,std::vector<int> >::iterator myit;
//     for (myit=patch_components_lm.begin(); myit != patch_components_lm.end(); myit++)
//     {
//         std::cout << XFEM::PHYSICS::physVarToString(myit->first) << std::endl;
//         std::cout<< (myit->second).size() << std::endl;
//         for (std::size_t i = 0; i< (myit->second).size(); i++)
//           std::cout<< (myit->second)[i] << std::endl;
//     }

  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                         rasthofer 02/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3IntFace::Print(std::ostream& os) const
{
  os << "Combust3IntFace ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                        rasthofer 02/13 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Combust3IntFace::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of Combust3IntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                        rasthofer 02/13 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Combust3IntFace::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of Combust3IntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}
