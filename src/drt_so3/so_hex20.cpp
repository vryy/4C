/*----------------------------------------------------------------------*/
/*! \file

\brief 3D quadratic serendipity element

\level 1


*----------------------------------------------------------------------*/

#include "so_hex20.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/prestress_service.H"
#include "../drt_mat/so3_material.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

DRT::ELEMENTS::So_hex20Type DRT::ELEMENTS::So_hex20Type::instance_;

DRT::ELEMENTS::So_hex20Type& DRT::ELEMENTS::So_hex20Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_hex20Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_hex20* object = new DRT::ELEMENTS::So_hex20(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex20Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH20")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex20(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex20Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex20(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_hex20Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_hex20Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::So_hex20Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH20"];

  defs["HEX20"]
      .AddIntVector("HEX20", 20)
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


/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex20::So_hex20(int id, int owner)
    : So_base(id, owner), data_(), pstype_(INPAR::STR::PreStress::none), pstime_(0.0), time_(0.0)
{
  invJ_.resize(NUMGPT_SOH20, LINALG::Matrix<NUMDIM_SOH20, NUMDIM_SOH20>(true));
  detJ_.resize(NUMGPT_SOH20, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = ::UTILS::PRESTRESS::GetType();
    pstime_ = ::UTILS::PRESTRESS::GetTime();
  }
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH20, NUMGPT_SOH20));

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex20::So_hex20(const DRT::ELEMENTS::So_hex20& old)
    : So_base(old),
      data_(old.data_),
      detJ_(old.detJ_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH20 x NUMDIM_SOH20?
    // invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
    invJ_[i] = old.invJ_[i];
  }

  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex20::Clone() const
{
  DRT::ELEMENTS::So_hex20* newelement = new DRT::ELEMENTS::So_hex20(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_hex20::Shape() const { return hex20; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);
  // data_
  AddtoPack(data, data_);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);
  // Pack prestress
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::Unpack(const std::vector<char>& data)
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
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOH20, NUMDIM_SOH20>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // Extract prestress
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH20, NUMGPT_SOH20));
    prestress_->Unpack(tmpprestress);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex20::~So_hex20() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::Print(std::ostream& os) const
{
  os << "So_hex20 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      popp 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_expol(
    LINALG::Matrix<NUMGPT_SOH20, MAT::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  static LINALG::Matrix<NUMNOD_SOH20, NUMGPT_SOH20> expol;
  static bool isfilled;
  LINALG::Matrix<NUMNOD_SOH20, MAT::NUM_STRESS_3D> nodalstresses;

  if (isfilled == true)
  {
    nodalstresses.Multiply(expol, stresses);
  }
  else
  {
    // get gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_27point);

    // loop over all nodes
    for (int ip = 0; ip < NUMNOD_SOH20; ++ip)
    {
      // gaussian coordinates
      const double e1 = intpoints.qxg[ip][0];
      const double e2 = intpoints.qxg[ip][1];
      const double e3 = intpoints.qxg[ip][2];

      // coordinates of node in the fictitious GP element
      double e1expol;
      double e2expol;
      double e3expol;

      if (e1 != 0)
        e1expol = 1 / e1;
      else
        e1expol = 0;
      if (e2 != 0)
        e2expol = 1 / e2;
      else
        e2expol = 0;
      if (e3 != 0)
        e3expol = 1 / e3;
      else
        e3expol = 0;

      // shape functions for the extrapolated coordinates
      // (yes, we REALLY mean DiscretizationType hex27 here!)
      LINALG::Matrix<NUMGPT_SOH20, 1> funct;
      DRT::UTILS::shape_function_3D(funct, e1expol, e2expol, e3expol, hex27);

      // extrapolation matrix
      for (int i = 0; i < NUMGPT_SOH20; ++i) expol(ip, i) = funct(i);
    }

    // do extrapolation
    nodalstresses.Multiply(expol, stresses);

    isfilled = true;
  }

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_SOH20; ++i)
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


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                           |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex20::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex20::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex20::Lines()
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
 |  Return names of visualization data (public)                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);
  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                                  |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex20::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOH20, this->Id());
}
