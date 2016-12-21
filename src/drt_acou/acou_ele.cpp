/*!----------------------------------------------------------------------
\file acou_ele.cpp
\brief acoustic element implementations

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#include "acou_ele.H"
#include "acou_ele_boundary_calc.H"
#include "acou_ele_intfaces_calc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"


DRT::ELEMENTS::AcouType DRT::ELEMENTS::AcouType::instance_;
DRT::ELEMENTS::AcouBoundaryType DRT::ELEMENTS::AcouBoundaryType::instance_;
DRT::ELEMENTS::AcouIntFaceType DRT::ELEMENTS::AcouIntFaceType::instance_;

DRT::ELEMENTS::AcouType& DRT::ELEMENTS::AcouType::Instance()
{
  return instance_;
}

DRT::ELEMENTS::AcouBoundaryType& DRT::ELEMENTS::AcouBoundaryType::Instance()
{
  return instance_;
}

DRT::ELEMENTS::AcouIntFaceType& DRT::ELEMENTS::AcouIntFaceType::Instance()
{
  return instance_;
}


DRT::ParObject* DRT::ELEMENTS::AcouType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Acou* object = new DRT::ELEMENTS::Acou(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouType::Create(const std::string  eletype,
                                                           const std::string  eledistype,
                                                           const int     id,
                                                           const int     owner)
{
  if ( eletype=="ACOUSTIC" )
  {
    return Teuchos::rcp(new DRT::ELEMENTS::Acou(id,owner));
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::Acou(id,owner));
}


void DRT::ELEMENTS::AcouType::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  nv = DRT::UTILS::getDimension(dwele->Shape());
  np = 1;
  dimns = nv+np;
  numdf = dimns;
  return;
}


void DRT::ELEMENTS::AcouType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  // DRT::UTILS::ComputeFluidDNullSpace( dis, ns, x0, numdf, dimns );
  return;
}


void DRT::ELEMENTS::AcouType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["ACOUSTIC"];

  // 3D elements
  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedInt("DEG")
    .AddNamedInt("SPC")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedInt("DEG")
    .AddNamedInt("SPC")
    ;

  // 2D elements
  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedInt("DEG")
    .AddNamedInt("SPC")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedInt("DEG")
    .AddNamedInt("SPC")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedInt("DEG")
    .AddNamedInt("SPC")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                         schoeder 07/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Acou::Acou(int id, int owner) :
DRT::Element(id,owner),
degree_(1),
completepol_(true)
{
  distype_= dis_none;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schoeder 07/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Acou::Acou(const DRT::ELEMENTS::Acou& old) :
DRT::Element(old             ),
distype_    (old.distype_    ),
degree_(old.degree_),
completepol_(old.completepol_)
{
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Acou and return pointer to it (public)   |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Acou::Clone() const
{
  DRT::ELEMENTS::Acou* newelement = new DRT::ELEMENTS::Acou(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Acou::Pack(DRT::PackBuffer& data) const
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
  int degree = degree_;
  AddtoPack(data, degree);
  AddtoPack(data, completepol_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Acou::Unpack(const std::vector<char>& data)
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
  int val = 0;
  ExtractfromPack(position,data,val);
  dsassert(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  ExtractfromPack(position,data,val);
  completepol_ = val;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                         schoeder 07/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Acou::~Acou()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                           schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Acou::Print(std::ostream& os) const
{
  os << "Acou ";
  Element::Print(os);
  return;
}


bool DRT::ELEMENTS::Acou::ReadElement(const std::string& eletype,
                                      const std::string& distype,
                                      DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);
  int degree;
  linedef->ExtractInt("DEG", degree);
  degree_ = degree;

  linedef->ExtractInt("SPC", degree);
  completepol_ = degree;

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(DRT::StringToDistype(distype));

  return true;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)             schoeder 07/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Acou::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<AcouBoundary,Acou>(DRT::UTILS::buildLines,this);
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
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Acou::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<AcouBoundary,Acou>(DRT::UTILS::buildSurfaces,this);
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
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Acou::Volumes()
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
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Acou::CreateFaceElement( DRT::Element* parent_slave,           //!< parent slave fluid3 element
                                                                   int nnode,                            //!< number of surface nodes
                                                                   const int* nodeids,                   //!< node ids of surface element
                                                                   DRT::Node** nodes,                    //!< nodes of surface element
                                                                   const int lsurface_master,            //!< local surface number w.r.t master parent element
                                                                   const int lsurface_slave,             //!< local surface number w.r.t slave parent element
                                                                   const std::vector<int> &localtrafomap //! local trafo map
)
{
  // dynamic cast for slave parent element
  DRT::ELEMENTS::Acou * slave_pele = dynamic_cast<DRT::ELEMENTS::Acou *>( parent_slave );

  // insert both parent elements
  return DRT::UTILS::ElementIntFaceFactory<AcouIntFace,Acou>( -1,             //!< internal face element id
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouBoundaryType::Create( const int id, const int owner )
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                        schoeder 07/13 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouBoundary::AcouBoundary(int id, int owner,
                                          int nnode, const int* nodeids,
                                          DRT::Node** nodes,
                                          DRT::ELEMENTS::Acou* parent,
                                          const int lsurface) :
DRT::FaceElement(id,owner)
{
  SetParentMasterElement(parent,lsurface);
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   schoeder 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouBoundary::AcouBoundary(const DRT::ELEMENTS::AcouBoundary& old) :
DRT::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                       schoeder 07/13 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::AcouBoundary::Clone() const
{
  DRT::ELEMENTS::AcouBoundary* newelement = new DRT::ELEMENTS::AcouBoundary(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::AcouBoundary::Shape() const
{
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                        schoeder 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouBoundary::Pack(DRT::PackBuffer& data) const
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
void DRT::ELEMENTS::AcouBoundary::Unpack(const std::vector<char>& data)
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
DRT::ELEMENTS::AcouBoundary::~AcouBoundary()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                          schoeder 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouBoundary::Print(std::ostream& os) const
{
  os << "AcouBoundary ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         schoeder 07/13 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouBoundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of AcouBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         schoeder 07/13 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouBoundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of AcouBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        schoeder 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::AcouBoundary::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::AcouBoundaryImplInterface::Impl(this)->Evaluate(this,params,discretization,lm,elemat1,elemat2,elevec1,elevec2,elevec3);
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition  schoeder 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::AcouBoundary::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  dserror("dummy function called");
  return 0;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) schoeder 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouBoundary::LocationVector(
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::AcouIntFaceType::Create( const int id, const int owner )
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                         schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouIntFace::AcouIntFace(int id,                                ///< element id
                                          int owner,                             ///< owner (= owner of parent element with smallest gid)
                                          int nnode,                             ///< number of nodes
                                          const int* nodeids,                    ///< node ids
                                          DRT::Node** nodes,                     ///< nodes of surface
                                          DRT::ELEMENTS::Acou* parent_master,   ///< master parent element
                                          DRT::ELEMENTS::Acou* parent_slave,    ///< slave parent element
                                          const int lsurface_master,             ///< local surface index with respect to master parent element
                                          const int lsurface_slave,              ///< local surface index with respect to slave parent element
                                          const std::vector<int> localtrafomap   ///< get the transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system and the slave element's face's coordinate system
):
DRT::FaceElement(id,owner)
{
  SetParentMasterElement(parent_master,lsurface_master);
  SetParentSlaveElement(parent_slave,lsurface_slave);

  if(parent_slave != NULL)
    degree_ = std::max(parent_master->Degree(),parent_slave->Degree());
  else
    degree_ = parent_master->Degree();

  SetLocalTrafoMap(localtrafomap);

  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouIntFace::AcouIntFace(const DRT::ELEMENTS::AcouIntFace& old) :
DRT::FaceElement(old),
degree_(old.degree_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                        schoeder 01/14|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::AcouIntFace::Clone() const
{
  DRT::ELEMENTS::AcouIntFace* newelement = new DRT::ELEMENTS::AcouIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::AcouIntFace::Shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouIntFace::Pack(DRT::PackBuffer& data) const
{
  dserror("this AcouIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouIntFace::Unpack(const std::vector<char>& data)
{
  dserror("this AcouIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouIntFace::~AcouIntFace()
{
  return;
}


/*----------------------------------------------------------------------*
 |  create the patch location vector (public)            schoeder 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouIntFace::PatchLocationVector(
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
  // *this AcouIntFace element only once (no duplicates)

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
void DRT::ELEMENTS::AcouIntFace::Print(std::ostream& os) const
{
  os << "AcouIntFace ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouIntFace::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of AcouIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::AcouIntFace::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of AcouIntFace not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::AcouIntFace::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // REMARK: this line ensures that the static DRT::ELEMENTS::AcouIntFaceImplInterface::Impl is created
  //         this line avoids linker errors
  DRT::ELEMENTS::AcouIntFaceImplInterface::Impl(this);

  dserror("not available");

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition  schoeder 01/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::AcouIntFace::EvaluateNeumann(
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


