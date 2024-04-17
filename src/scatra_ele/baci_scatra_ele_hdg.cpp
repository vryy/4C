/*----------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG element

\level 3

*----------------------------------------------------------------------*/

#include "baci_scatra_ele_hdg.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_discretization_fem_general_utils_polynomial.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret_faces.hpp"
#include "baci_lib_discret_hdg.hpp"
#include "baci_mat_list.hpp"
#include "baci_mat_myocard.hpp"
#include "baci_scatra_ele_action.hpp"
#include "baci_scatra_ele_factory.hpp"
#include "baci_scatra_ele_hdg_boundary_calc.hpp"
#include "baci_scatra_ele_hdg_intfaces_calc.hpp"
#include "baci_scatra_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN


// initialize static variable
DRT::ELEMENTS::ScaTraHDGType DRT::ELEMENTS::ScaTraHDGType::instance_;
DRT::ELEMENTS::ScaTraHDGBoundaryType DRT::ELEMENTS::ScaTraHDGBoundaryType::instance_;
DRT::ELEMENTS::ScaTraHDGIntFaceType DRT::ELEMENTS::ScaTraHDGIntFaceType::instance_;


DRT::ELEMENTS::ScaTraHDGType& DRT::ELEMENTS::ScaTraHDGType::Instance() { return instance_; }

DRT::ELEMENTS::ScaTraHDGBoundaryType& DRT::ELEMENTS::ScaTraHDGBoundaryType::Instance()
{
  return instance_;
}

DRT::ELEMENTS::ScaTraHDGIntFaceType& DRT::ELEMENTS::ScaTraHDGIntFaceType::Instance()
{
  return instance_;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::ScaTraHDGType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::ScaTraHDG* object = new DRT::ELEMENTS::ScaTraHDG(-1, -1);
  object->Unpack(data);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ScaTraHDGType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRANSPHDG")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::ScaTraHDG(id, owner));
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ScaTraHDGType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::ScaTraHDG(id, owner));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 1;  // Only one scalar (so far) is the unknown that is solved for
  dimns = numdf;
  nv = numdf;

  if (GLOBAL::Problem::Instance(0)->GetProblemType() == GLOBAL::ProblemType::elch)
  {
    if (nv > 1)  // only when we have more than 1 dof per node!
    {
      nv -= 1;  // ion concentrations
      np = 1;   // electric potential
    }
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::ScaTraHDGType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  dserror("method ComputeNullSpace not implemented right now!");
  return nullspace;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGType ::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_scatra;
  TransportType::SetupElementDefinition(definitions_scatra);

  std::map<std::string, INPUT::LineDefinition>& defs_scatra = definitions_scatra["TRANSP"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["TRANSPHDG"];

  for (const auto& [key, scatra_line_def] : defs_scatra)
  {
    defs[key] = INPUT::LineDefinition::Builder(scatra_line_def)
                    .AddNamedInt("DEG")
                    .AddOptionalNamedInt("SPC")
                    .Build();
  }
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                         hoermann 09/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDG::ScaTraHDG(int id, int owner)
    : Transport(id, owner),
      diff1_(0.0),
      ndofs_(0),
      onfdofs_(0),
      onfdofs_old_(0),
      degree_(1),
      degree_old_(0),
      completepol_(true),
      padpatele_(true),
      matinit_(false)
{
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    hoermann 09/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDG::ScaTraHDG(const DRT::ELEMENTS::ScaTraHDG& old)
    : Transport(old),
      diff1_(0.0),
      ndofs_(0),
      onfdofs_(0),
      onfdofs_old_(0),
      degree_(old.degree_),
      degree_old_(old.degree_old_),
      completepol_(old.completepol_),
      padpatele_(true),
      matinit_(false)
{
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of ScaTra and return pointer to it (public) |
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ScaTraHDG::Clone() const
{
  DRT::ELEMENTS::ScaTraHDG* newelement = new DRT::ELEMENTS::ScaTraHDG(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  Pack data (public)                                   hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDG::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Transport::Pack(data);

  int degree = degree_;
  AddtoPack(data, degree);
  degree = completepol_;
  AddtoPack(data, degree);
  degree = degree_old_;
  AddtoPack(data, degree);
}



/*----------------------------------------------------------------------*
 |  Unpack data (public)                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDG::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  Transport::ExtractfromPack(position, data, basedata);
  Transport::Unpack(basedata);

  int val = 0;
  ExtractfromPack(position, data, val);
  dsassert(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  ExtractfromPack(position, data, val);
  completepol_ = val;
  ExtractfromPack(position, data, val);
  degree_old_ = val;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  PackMaterial data (public)                           hoermann 12/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDG::PackMaterial(CORE::COMM::PackBuffer& data) const
{
  // add material
  if (Material() != Teuchos::null)
  {
    // pack only first material
    Material()->Pack(data);
  }
  else
    dserror("No material defined to pack!");
}

/*----------------------------------------------------------------------*
 |  UnPackMaterial data (public)                         hoermann 12/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDG::UnpackMaterial(const std::vector<char>& data) const
{
  Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType() == INPAR::MAT::m_myocard)
  {
    // Note: We need to do a dynamic_cast here
    Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(mat);
    actmat->UnpackMaterial(data);
  }
  else
    dserror("No material defined to unpack!");
}

/*----------------------------------------------------------------------*
 |  init the element                                     hoermann 12/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDG::Initialize()
{
  Teuchos::RCP<MAT::Material> mat = Material();
  // for now, we only need to do something in case of reactions (for the initialization of functions
  // in case of reactions by function)
  if (mat->MaterialType() == INPAR::MAT::m_myocard)
  {
    int gp;
    // Note: We need to do a dynamic_cast here
    Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(mat);
    int deg = 0;
    if (degree_old_ == 1)
      deg = 4 * degree_old_;
    else
      deg = 3 * degree_old_;
    if (this->Shape() == CORE::FE::CellType::tet4 or this->Shape() == CORE::FE::CellType::tet10)
    {
      switch (deg)
      {
        case 0:
          gp = 1;
          break;
        case 3:
          gp = 5;
          break;
        case 4:
          gp = 11;
          break;
        case 6:
          gp = 24;
          break;
        case 9:
          gp = 125;
          break;
        case 12:
          gp = 343;
          break;
        case 15:
          gp = 729;
          break;
        default:
          dserror(
              "Integration rule for TET elements only until polynomial order 5 for TET defined. "
              "You specified a degree of %d ",
              degree_old_);
          gp = 0;
          break;
      }
    }
    else
    {
      Teuchos::RCP<CORE::FE::GaussPoints> quadrature_(
          CORE::FE::GaussPointCache::Instance().Create(this->Shape(), deg));
      gp = quadrature_->NumPoints();
    }
    if (actmat->Parameter() != nullptr and
        !actmat->MyocardMat())  // in case we are not in post-process mode
    {
      actmat->SetGP(gp);
      actmat->Initialize();
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  Read element from input (public)                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::ScaTraHDG::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  bool success = Transport::ReadElement(eletype, distype, linedef);
  int degree;
  linedef->ExtractInt("DEG", degree);
  degree_ = degree;
  degree_old_ = degree_;

  if (linedef->HaveNamed("SPC"))
  {
    linedef->ExtractInt("SPC", degree);
    completepol_ = degree;
  }
  else
    completepol_ = false;

  return success;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)             hoermann 09/15|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDG::Lines()
{
  return CORE::COMM::GetElementLines<ScaTraHDGBoundary, ScaTraHDG>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                       hoermann 09/15|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDG::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<ScaTraHDGBoundary>(*this);
}


/*----------------------------------------------------------------------*
 |  get face element (public)                             hoermann 09/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ScaTraHDG::CreateFaceElement(
    DRT::Element* parent_slave,            //!< parent slave fluid3 element
    int nnode,                             //!< number of surface nodes
    const int* nodeids,                    //!< node ids of surface element
    DRT::Node** nodes,                     //!< nodes of surface element
    const int lsurface_master,             //!< local surface number w.r.t master parent element
    const int lsurface_slave,              //!< local surface number w.r.t slave parent element
    const std::vector<int>& localtrafomap  //! local trafo map
)
{
  // dynamic cast for slave parent element
  DRT::ELEMENTS::ScaTraHDG* slave_pele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(parent_slave);


  // insert both parent elements
  return CORE::COMM::ElementIntFaceFactory<ScaTraHDGIntFace, ScaTraHDG>(
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


/*---------------------------------------------------------------------*
|  evaluate the element (public)                         hoermann 09/15|
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDG::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  // get the action required
  SCATRA::Action act;
  if (params.get<bool>("hdg_action", false))
  {
    switch (Teuchos::getIntegralValue<DRT::HDGAction>(params, "action"))
    {
      case DRT::HDGAction::project_dirich_field:
        act = SCATRA::Action::project_dirich_field;
        break;
      default:
        dserror("HDG Action type not supported");
    }
  }
  else
  {
    act = Teuchos::getIntegralValue<SCATRA::Action>(params, "action");
  }

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();

  // switch between different physical types as used below
  switch (act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case SCATRA::Action::calc_mat_and_rhs:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImplHDG(
          Shape(), ImplType(), numdofpernode, numscal, discretization.Name())
          ->Evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    break;

    case SCATRA::Action::interpolate_hdg_to_node:
    case SCATRA::Action::update_interior_variables:
    case SCATRA::Action::project_dirich_field:
    case SCATRA::Action::project_material_field:
    case SCATRA::Action::project_neumann_field:
    case SCATRA::Action::set_initial_field:
    case SCATRA::Action::time_update_material:
    case SCATRA::Action::get_material_internal_state:
    case SCATRA::Action::set_material_internal_state:
    case SCATRA::Action::calc_mat_initial:
    case SCATRA::Action::project_field:
    case SCATRA::Action::calc_padaptivity:
    case SCATRA::Action::calc_error:

    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImplHDG(
          Shape(), ImplType(), numdofpernode, numscal, discretization.Name())
          ->EvaluateService(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }

    case SCATRA::Action::calc_initial_time_deriv:
    case SCATRA::Action::set_general_scatra_parameter:
    case SCATRA::Action::set_nodeset_parameter:
    case SCATRA::Action::set_time_parameter:
    case SCATRA::Action::set_turbulence_scatra_parameter:
      break;

    default:
      dserror("Unknown type of action '%i' for ScaTraHDG", act);
      break;
  }  // switch(action)

  return 0;
}  // DRT::ELEMENTS::ScaTra::Evaluate


/*----------------------------------------------------------------------*
 |  print this element (public)                           hoermann 09/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDG::Print(std::ostream& os) const
{
  os << "ScaTraHDG ";
  Element::Print(os);
}



//===========================================================================



Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ScaTraHDGBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                        hoermann 09/15 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGBoundary::ScaTraHDGBoundary(int id, int owner, int nnode,
    const int* nodeids, DRT::Node** nodes, DRT::Element* parent, const int lsurface)
    : DRT::FaceElement(id, owner)
{
  SetParentMasterElement(parent, lsurface);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   hoermann 09/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGBoundary::ScaTraHDGBoundary(const DRT::ELEMENTS::ScaTraHDGBoundary& old)
    : DRT::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ScaTraHDGBoundary::Clone() const
{
  DRT::ELEMENTS::ScaTraHDGBoundary* newelement = new DRT::ELEMENTS::ScaTraHDGBoundary(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                        hoermann 09/15|
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::ScaTraHDGBoundary::Shape() const
{
  return CORE::FE::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                        hoermann 09/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGBoundary::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  // Discretisation type
  // AddtoPack(data,distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        hoermann 09/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGBoundary::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // distype
  // distype_ = static_cast<CORE::FE::CellType>( ExtractInt(position,data) );

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                          hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGBoundary::Print(std::ostream& os) const
{
  os << "ScaTraHDGBoundary ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         hoermann 09/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDGBoundary::Lines()
{
  dserror("Lines of ScaTraHDGBoundary not implemented");
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         hoermann 09/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDGBoundary::Surfaces()
{
  dserror("Surfaces of ScaTraHDGBoundary not implemented");
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        hoermann 09/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDGBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition  hoermann 09/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDGBoundary::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)

{
  // add Neumann boundary condition to parameter list
  params.set<DRT::Condition*>("condition", &condition);

  // build location array from location vector
  //(this a little ugly. one could fix this by introducing a EvaluateNeumann() method
  // with LocationArray as input in the DRT::Element ...)
  LocationArray la(1);
  la[0].lm_ = lm;

  DRT::ELEMENTS::ScaTraHDGBoundaryImplInterface::Impl(this)->EvaluateNeumann(
      this, params, discretization, la, *elemat1, elevec1);

  return 0;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGBoundary::LocationVector(const Discretization& dis, LocationArray& la,
    bool doDirichlet, const std::string& condstring, Teuchos::ParameterList& params) const
{
  // we have to do it this way
  ParentMasterElement()->LocationVector(dis, la, false);
  return;
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element (public) hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGBoundary::LocationVector(const Discretization& dis,
    std::vector<int>& lm, std::vector<int>& lmowner, std::vector<int>& lmstride) const
{
  // we have to do it this way
  ParentMasterElement()->LocationVector(dis, lm, lmowner, lmstride);
  return;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ScaTraHDGIntFaceType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                         hoermann 09/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGIntFace::ScaTraHDGIntFace(int id,  ///< element id
    int owner,           ///< owner (= owner of parent element with smallest gid)
    int nnode,           ///< number of nodes
    const int* nodeids,  ///< node ids
    DRT::Node** nodes,   ///< nodes of surface
    DRT::ELEMENTS::ScaTraHDG* parent_master,  ///< master parent element
    DRT::ELEMENTS::ScaTraHDG* parent_slave,   ///< slave parent element
    const int lsurface_master,  ///< local surface index with respect to master parent element
    const int lsurface_slave,   ///< local surface index with respect to slave parent element
    const std::vector<int>
        localtrafomap  ///< get the transformation map between the local coordinate systems of the
                       ///< face w.r.t the master parent element's face's coordinate system and the
                       ///< slave element's face's coordinate system
    )
    : DRT::FaceElement(id, owner), degree_(0), degree_old_(0)
{
  SetParentMasterElement(parent_master, lsurface_master);
  SetParentSlaveElement(parent_slave, lsurface_slave);

  if (parent_slave != nullptr)
  {
    degree_ = std::max(parent_master->Degree(), parent_slave->Degree());
    degree_old_ = std::max(parent_master->DegreeOld(), parent_slave->DegreeOld());
  }
  else
  {
    degree_ = parent_master->Degree();
    degree_old_ = parent_master->DegreeOld();
  }

  SetLocalTrafoMap(localtrafomap);

  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    hoermann 09/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGIntFace::ScaTraHDGIntFace(const DRT::ELEMENTS::ScaTraHDGIntFace& old)
    : DRT::FaceElement(old), degree_(old.degree_), degree_old_(old.degree_old_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                        hoermann 09/15|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ScaTraHDGIntFace::Clone() const
{
  DRT::ELEMENTS::ScaTraHDGIntFace* newelement = new DRT::ELEMENTS::ScaTraHDGIntFace(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::ScaTraHDGIntFace::Shape() const
{
  // could be called for master parent or slave parent element, doesn't matter
  return CORE::FE::getShapeOfBoundaryElement(NumNode(), ParentMasterElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGIntFace::Pack(CORE::COMM::PackBuffer& data) const
{
  dserror("this ScaTraHDGIntFace element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGIntFace::Unpack(const std::vector<char>& data)
{
  dserror("this ScaTraHDGIntFace element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  create the patch location vector (public)            hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGIntFace::PatchLocationVector(
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
    std::vector<int>& lm_slaveNodeToPatch    ///< local map between slave nodes and nodes in patch
)
{
  // create one patch location vector containing all dofs of master, slave and
  // *this ScaTraHDGIntFace element only once (no duplicates)

  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();

  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    dserror("wrong number of nodes for master element");
  }

  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    dserror("wrong number of nodes for slave element");
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
  int dofset = 0;  // assume dofset 0

  int patchnode_count = 0;

  // fill patch lm with master's nodes
  for (int k = 0; k < m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof = discretization.Dof(dofset, node);

    // get maximum of numdof per node with the help of master and/or slave element (returns 4 in 3D
    // case, does not return dofset's numnode)
    const int size = discretization.NumDof(dofset, node);
    const int offset = size * nds_master[k];

    dsassert(dof.size() >= static_cast<unsigned>(offset + size), "illegal physical dofs offset");

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
    DRT::Node* node = s_nodes[k];

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

      dsassert(dof.size() >= static_cast<unsigned>(offset + size), "illegal physical dofs offset");
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
      dserror("face's nodes not contained in masternodes_offset map");
  }

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                          hoermann 09/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraHDGIntFace::Print(std::ostream& os) const
{
  os << "ScaTraHDGIntFace ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         hoermann 09/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDGIntFace::Lines()
{
  dserror("Lines of ScaTraHDGIntFace not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         hoermann 09/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::ScaTraHDGIntFace::Surfaces()
{
  dserror("Surfaces of ScaTraHDGIntFace not implemented");
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        hoermann 09/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDGIntFace::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // REMARK: this line ensures that the static DRT::ELEMENTS::ScaTraHDGIntFaceImplInterface::Impl is
  // created
  //         this line avoids linker errors
  DRT::ELEMENTS::ScaTraHDGIntFaceImplInterface::Impl(this);

  dserror("not available");

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition  hoermann 09/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ScaTraHDGIntFace::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  dserror("not available");

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
