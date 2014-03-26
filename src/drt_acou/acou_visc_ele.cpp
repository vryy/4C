/*!----------------------------------------------------------------------
\file acou_visc_ele.cpp
\brief

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "acou_visc_ele.H"
#include "acou_ele_boundary_calc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::AcouViscType DRT::ELEMENTS::AcouViscType::instance_;
DRT::ELEMENTS::AcouViscBoundaryType DRT::ELEMENTS::AcouViscBoundaryType::instance_;
DRT::ELEMENTS::AcouViscIntFaceType DRT::ELEMENTS::AcouViscIntFaceType::instance_;

/*----------------------------------------------------------------------*
 |                                                        schoeder 02/14|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::AcouViscType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::AcouVisc* object = new DRT::ELEMENTS::AcouVisc(-1,-1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schoeder 02/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouViscType::Create(const std::string  eletype,
                                                             const std::string  eledistype,
                                                             const int     id,
                                                             const int     owner)
{
  if (eletype=="ACOUSTICVISC")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::AcouVisc(id,owner));
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schoeder 02/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouViscType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::AcouVisc(id,owner));
}

/*----------------------------------------------------------------------*
 |                                                        schoeder 02/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["ACOUSTICVISC"];

  // 3D elements
  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  // 2D elements
  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

}

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schoeder 02/14|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouVisc::AcouVisc(int id, int owner) :
Acou(id,owner)
{
  distype_= dis_none;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schoeder 02/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouVisc::AcouVisc(const DRT::ELEMENTS::AcouVisc& old) :
Acou(old             ),
distype_    (old.distype_    )
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Acou and return pointer to it (public)   |
 |                                                        schoeder 02/14|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::AcouVisc::Clone() const
{
  DRT::ELEMENTS::AcouVisc* newelement = new DRT::ELEMENTS::AcouVisc(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouVisc::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  // Discretisation type
  AddtoPack(data,distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouVisc::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // distype
  distype_ = static_cast<DiscretizationType>( ExtractInt(position,data) );

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                         schoeder 02/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouVisc::~AcouVisc()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                           schoeder 02/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouVisc::Print(std::ostream& os) const
{
  os << "AcouVisc ";
  Element::Print(os);
  return;
}

bool DRT::ELEMENTS::AcouVisc::ReadElement(const std::string& eletype,
                                      const std::string& distype,
                                      DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(DRT::StringToDistype(distype));

  return true;
}



/*----------------------------------------------------------------------*
 |  get vector of lines              (public)             schoeder 07/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouVisc::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<AcouViscBoundary,AcouVisc>(DRT::UTILS::buildLines,this);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                       schoeder 07/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouVisc::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<AcouViscBoundary,AcouVisc>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1) // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else  // 1D elements
  {
    dserror("Surfaces() does not exist for 1D-element ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)             schoeder 07/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouVisc::Volumes()
{
  if (NumVolume()==1) // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element> > volumes(1);
    volumes[0]= Teuchos::rcp(this, false);
    return volumes;
  }
  else //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Surfaces();
  }
}



/*----------------------------------------------------------------------*
 |  get face element (public)                             schoeder 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouVisc::CreateFaceElement( DRT::Element* parent_slave,           //!< parent slave fluid3 element
                                                                   int nnode,                            //!< number of surface nodes
                                                                   const int* nodeids,                   //!< node ids of surface element
                                                                   DRT::Node** nodes,                    //!< nodes of surface element
                                                                   const int lsurface_master,            //!< local surface number w.r.t master parent element
                                                                   const int lsurface_slave,             //!< local surface number w.r.t slave parent element
                                                                   const std::vector<int> &localtrafomap //! local trafo map
)
{
  // dynamic cast for slave parent element
  DRT::ELEMENTS::AcouVisc * slave_pele = dynamic_cast<DRT::ELEMENTS::AcouVisc *>( parent_slave );

  // insert both parent elements
  return DRT::UTILS::ElementIntFaceFactory<AcouViscIntFace,AcouVisc>( -1,             //!< internal face element id
                                                              -1,             //!< owner of internal face element
                                                              nnode,          //!< number of surface nodes
                                                              nodeids,        //!< node ids of surface element
                                                              nodes,          //!< nodes of surface element
                                                              this,           //!< master parent element
                                                              slave_pele,     //!< slave parent element
                                                              lsurface_master,//!< local surface number w.r.t master parent element
                                                              lsurface_slave, //!< local surface number w.r.t slave parent element
                                                              localtrafomap   //!< local trafo map
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouViscBoundaryType::Create( const int id, const int owner )
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                        schoeder 07/13 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscBoundary::AcouViscBoundary(int id, int owner,
                                          int nnode, const int* nodeids,
                                          DRT::Node** nodes,
                                          DRT::ELEMENTS::AcouVisc* parent,
                                          const int lsurface) :
AcouBoundary(id,owner,nnode,nodeids,nodes,parent,lsurface)
{
//  SetParentMasterElement(parent,lsurface);
//  SetNodeIds(nnode,nodeids);
//  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   schoeder 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscBoundary::AcouViscBoundary(const DRT::ELEMENTS::AcouViscBoundary& old) :
AcouBoundary(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                       schoeder 07/13 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::AcouViscBoundary::Clone() const
{
  DRT::ELEMENTS::AcouViscBoundary* newelement = new DRT::ELEMENTS::AcouViscBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscBoundary::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  // Discretisation type
  //AddtoPack(data,distype_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscBoundary::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // distype
  //distype_ = static_cast<DiscretizationType>( ExtractInt(position,data) );

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                         schoeder 07/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscBoundary::~AcouViscBoundary()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                          schoeder 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscBoundary::Print(std::ostream& os) const
{
  os << "AcouViscBoundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        schoeder 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::AcouViscBoundary::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // TODO
  DRT::ELEMENTS::AcouBoundaryImplInterface::Impl(this)->Absorbing(this,params,discretization,lm,elemat1,elemat2,elevec1,elevec2,elevec3);
  return 0;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) schoeder 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscBoundary::LocationVector(
    const Discretization&    dis,
    LocationArray&           la,
    bool                     doDirichlet,
    const std::string&       condstring,
    Teuchos::ParameterList&  params
    ) const
{
  // we have to do it this way, just as for weak Dirichlet conditions
  ParentMasterElement()->LocationVector(dis,la,false);
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouViscIntFaceType::Create( const int id, const int owner )
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscIntFace::AcouViscIntFace(int id,                                ///< element id
                                          int owner,                             ///< owner (= owner of parent element with smallest gid)
                                          int nnode,                             ///< number of nodes
                                          const int* nodeids,                    ///< node ids
                                          DRT::Node** nodes,                     ///< nodes of surface
                                          DRT::ELEMENTS::AcouVisc* parent_master,   ///< master parent element
                                          DRT::ELEMENTS::AcouVisc* parent_slave,    ///< slave parent element
                                          const int lsurface_master,             ///< local surface index with respect to master parent element
                                          const int lsurface_slave,              ///< local surface index with respect to slave parent element
                                          const std::vector<int> localtrafomap   ///< get the transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system and the slave element's face's coordinate system
):
AcouIntFace(id,owner,nnode,nodeids,nodes,parent_master,parent_slave,lsurface_master,lsurface_slave,localtrafomap),
localtrafomap_(localtrafomap)
{

  SetParentMasterElement(parent_master,lsurface_master);
  SetParentSlaveElement(parent_slave,lsurface_slave);
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscIntFace::AcouViscIntFace(const DRT::ELEMENTS::AcouViscIntFace& old) :
AcouIntFace(old),
localtrafomap_(old.localtrafomap_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                        schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::AcouViscIntFace::Clone() const
{
  DRT::ELEMENTS::AcouViscIntFace* newelement = new DRT::ELEMENTS::AcouViscIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::AcouViscIntFace::Shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouViscIntFace::~AcouViscIntFace()
{
  return;
}

/*----------------------------------------------------------------------*
 |  create the patch location vector (public)            schoeder 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscIntFace::PatchLocationVector(
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
  // *this AcouViscIntFace element only once (no duplicates)

  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();

  if ( m_numnode != static_cast<int>( nds_master.size() ) )
  {
    throw std::runtime_error( "wrong number of nodes for master element" );
  }

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

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
    const int size = discretization.NumDof(dofset,node);
    const int offset = size*nds_master[k];

    dsassert ( dof.size() >= static_cast<unsigned>( offset+size ), "illegal physical dofs offset" );

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
      const int size = discretization.NumDof(dofset,node);
      const int offset = size*nds_slave[k];

      dsassert ( dof.size() >= static_cast<unsigned>( offset+size ), "illegal physical dofs offset" );
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
      const int size = discretization.NumDof(dofset,node);

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
      const int size = discretization.NumDof(dofset,node);

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
 |  print this element (public)                          schoeder 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouViscIntFace::Print(std::ostream& os) const
{
  os << "AcouViscIntFace ";
  Element::Print(os);
  return;
}



