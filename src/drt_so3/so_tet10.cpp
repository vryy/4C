/*!----------------------------------------------------------------------
\file so_tet10.cpp

\brief ToDo Add meaningful comment.

\level 1

\maintainer Jonas Biehler

*----------------------------------------------------------------------*/

#include "so_tet10.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_mat/so3_material.H"

// inverse design object
#include "inversedesign.H"
#include "prestress.H"
// remove later


DRT::ELEMENTS::So_tet10Type DRT::ELEMENTS::So_tet10Type::instance_;

DRT::ELEMENTS::So_tet10Type& DRT::ELEMENTS::So_tet10Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_tet10Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_tet10* object = new DRT::ELEMENTS::So_tet10(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT10")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_tet10(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_tet10(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_tet10Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_tet10Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::So_tet10Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT10"];

  defs["TET10"]
      .AddIntVector("TET10", 10)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3)
      .AddOptionalNamedDouble("STRENGTH")
      .AddOptionalNamedDouble("HU")
      .AddOptionalNamedDouble("lambda")
      .AddOptionalNamedDouble("GROWTHTRIG");
}


/*----------------------------------------------------------------------***
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::So_tet10(int id, int owner)
    : So_base(id, owner), data_(), pstype_(INPAR::STR::prestress_none), pstime_(0.0), time_(0.0)
{
  invJ_.resize(NUMGPT_SOTET10, LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>(true));
  detJ_.resize(NUMGPT_SOTET10, 0.0);
  invJ_mass_.resize(NUMGPT_MASS_SOTET10, LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>(true));
  detJ_mass_.resize(NUMGPT_MASS_SOTET10, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
  }

  if (pstype_ == INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET10, NUMGPT_SOTET10));

  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::So_tet10(const DRT::ELEMENTS::So_tet10& old)
    : So_base(old),
      data_(old.data_),
      detJ_(old.detJ_),
      detJ_mass_(old.detJ_mass_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
// try out later detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

  invJ_mass_.resize(old.invJ_mass_.size());
  for (int i = 0; i < (int)invJ_mass_.size(); ++i)
  {
    invJ_mass_[i] = old.invJ_mass_[i];
  }

  if (pstype_ == INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_tet10::Clone() const
{
  DRT::ELEMENTS::So_tet10* newelement = new DRT::ELEMENTS::So_tet10(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_tet10::Shape() const { return tet10; }

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);
  ;
  // data
  AddtoPack(data, data_);
  // detJ_
  AddtoPack(data, detJ_);
  AddtoPack(data, detJ_mass_);

  // invJ
  const int size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  const int size_mass = (int)invJ_mass_.size();
  AddtoPack(data, size_mass);
  for (int i = 0; i < size_mass; ++i) AddtoPack(data, invJ_mass_[i]);

  // prestress_
  AddtoPack(data, pstype_);
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (pstype_ == INPAR::STR::prestress_mulf)
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }

  return;
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  So_base::Unpack(basedata);
  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  // detJ_
  ExtractfromPack(position, data, detJ_);
  ExtractfromPack(position, data, detJ_mass_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  int size_mass = 0;
  ExtractfromPack(position, data, size_mass);
  invJ_mass_.resize(size_mass, LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>(true));
  for (int i = 0; i < size_mass; ++i) ExtractfromPack(position, data, invJ_mass_[i]);

  // prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (pstype_ == INPAR::STR::prestress_mulf)
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET10, NUMGPT_SOTET10));
    prestress_->Unpack(tmpprestress);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------***
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::~So_tet10() { return; }


/*----------------------------------------------------------------------***
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Print(std::ostream& os) const
{
  os << "So_tet10 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}



/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      popp 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_expol(
    LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  static LINALG::Matrix<NUMNOD_SOTET10, NUMGPT_SOTET10> expol;
  static bool isfilled;
  LINALG::Matrix<NUMNOD_SOTET10, MAT::NUM_STRESS_3D> nodalstresses;

  if (isfilled == true)
  {
    nodalstresses.Multiply(expol, stresses);
  }
  else
  {
    // get gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_tet_4point);

    // transformation of nodes of virtual tet4 element defined by GP
    const double palpha = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
    const double pbeta = (5.0 - sqrt(5.0)) / 20.0;
    // we need to setup a mapping which allows us to evaluate the position of the nodes
    // in the parameter space of a fictitious element defined by the GP

    // loop over all GP
    for (int ip = 0; ip < NUMGPT_SOTET10; ++ip)
    {
      // gaussian coordinates
      const double e1 = intpoints.qxg[ip][0];
      const double e2 = intpoints.qxg[ip][1];
      const double e3 = intpoints.qxg[ip][2];

      // coordinates of node in the fictitious GP element
      double e1expol;
      double e2expol;
      double e3expol;

      e1expol = (e1 - pbeta * (1 + palpha - pbeta)) / ((palpha - pbeta) * (palpha - pbeta));
      e2expol = (e2 - pbeta * (1 + palpha - pbeta)) / ((palpha - pbeta) * (palpha - pbeta));
      e3expol = (e3 - pbeta * (1 + palpha - pbeta)) / ((palpha - pbeta) * (palpha - pbeta));

      // shape functions for the extrapolated coordinates
      // (yes, we REALLY mean DiscretizationType tet4 here, because we are using 4 GP for stiffness
      // matrix integration )
      LINALG::Matrix<NUMNOD_SOTET10, 1> funct;
      DRT::UTILS::shape_function_3D(funct, e1expol, e2expol, e3expol, tet4);

      // extrapolation matrix
      // Evaluate shape functions (of fictious GP element) at the corner nodes of the actual element
      for (int i = 0; i < NUMGPT_SOTET10; ++i)
      {
        expol(ip, i) = funct(i);
      }
    }


    // linear combination
    // interpolate stresses at edge/element center nodes from linear combination
    // of corner nodes
    for (int i = 0; i < NUMGPT_SOTET10; ++i)
    {
      expol(4, i) = (0.5 * expol(0, i) + 0.5 * expol(1, i));
      expol(5, i) = (0.5 * expol(1, i) + 0.5 * expol(2, i));
      expol(6, i) = (0.5 * expol(0, i) + 0.5 * expol(2, i));
      expol(7, i) = (0.5 * expol(0, i) + 0.5 * expol(3, i));
      expol(8, i) = (0.5 * expol(1, i) + 0.5 * expol(3, i));
      expol(9, i) = (0.5 * expol(2, i) + 0.5 * expol(3, i));
    }
    // do extrapolation
    nodalstresses.Multiply(expol, stresses);

    isfilled = true;
  }

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    const int lid = expolstresses.Map().LID(NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = 1.0 / Nodes()[i]->NumElement();
      for (int j = 0; j < MAT::NUM_STRESS_3D; ++j)
        (*(expolstresses(j)))[lid] += nodalstresses(i, j) * invmyadjele;
    }
  }
  return;
}



/*====================================================================*/
/* 10-node tetrahedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (ksi1, ksi2, ksi3) of nodes
 * of a common tetrahedron [0,1]x[0,1]x[0,1]
 *  10-node hexahedron: node 0,1,...,9
 *
 * -----------------------
 *- this is the numbering used in GiD & EXODUS!!
 *      3-
 *      |\ ---
 *      |  \    --9
 *      |    \      ---
 *      |      \        -2
 *      |        \       /\
 *      |          \   /   \
 *      7            8      \
 *      |          /   \     \
 *      |        6       \    5
 *      |      /           \   \
 *      |    /               \  \
 *      |  /                   \ \
 *      |/                       \\
 *      0------------4-------------1
 */
/*====================================================================*/

/*----------------------------------------------------------------------***
 |  get vector of volumes (length 1) (public)                           |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_tet10::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}


/*----------------------------------------------------------------------**#
|  get vector of surfaces (public)                                     |
|  surface normals always point outward                                |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_tet10::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_tet10::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      DRT::UTILS::buildLines, this);
}
/*----------------------------------------------------------------------*
 |  get location of element center                              jb 08/11|
 *----------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::So_tet10::ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOTET10, 1> funct;
  // Centroid of a tet with (0,1)(0,1)(0,1) is (0.25, 0.25, 0.25)
  DRT::UTILS::shape_function_3D(funct, 0.25, 0.25, 0.25, distype);
  LINALG::Matrix<1, NUMDIM_SOTET10> midpoint;
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                 st 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);
  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                          st 01/10|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOTET10, this->Id());
}
