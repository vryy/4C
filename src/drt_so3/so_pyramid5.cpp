/*----------------------------------------------------------------------*/
/*! \file

\brief pyramid shaped solid element

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_pyramid5.H"
#include "so_pyramid5fbar.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/so3_material.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

DRT::ELEMENTS::So_pyramid5Type DRT::ELEMENTS::So_pyramid5Type::instance_;

DRT::ELEMENTS::So_pyramid5Type& DRT::ELEMENTS::So_pyramid5Type::Instance() { return instance_; }


DRT::ParObject* DRT::ELEMENTS::So_pyramid5Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_pyramid5* object = new DRT::ELEMENTS::So_pyramid5(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_pyramid5Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDP5")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_pyramid5(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_pyramid5Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_pyramid5(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_pyramid5Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_pyramid5Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::So_pyramid5Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDP5"];

  defs["PYRAMID5"]
      .AddIntVector("PYRAMID5", 5)
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
DRT::ELEMENTS::So_pyramid5::So_pyramid5(int id, int owner)
    : So_base(id, owner), data_(), pstype_(INPAR::STR::PreStress::none), pstime_(0.0), time_(0.0)
{
  kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  invJ_.resize(NUMGPT_SOP5, LINALG::Matrix<NUMDIM_SOP5, NUMDIM_SOP5>(true));
  detJ_.resize(NUMGPT_SOP5, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
  }
  if (pstype_ == INPAR::STR::PreStress::mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOP5, NUMGPT_SOP5));

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_pyramid5::So_pyramid5(const DRT::ELEMENTS::So_pyramid5& old)
    : So_base(old),
      kintype_(old.kintype_),
      data_(old.data_),
      detJ_(old.detJ_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOP5 x NUMDIM_SOP5?
    invJ_[i] = old.invJ_[i];
  }

  if (pstype_ == INPAR::STR::PreStress::mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_pyramid5::Clone() const
{
  DRT::ELEMENTS::So_pyramid5* newelement = new DRT::ELEMENTS::So_pyramid5(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_pyramid5::Shape() const { return pyramid5; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_pyramid5::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data, kintype_);
  // data_
  AddtoPack(data, data_);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  // prestress_
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (pstype_ == INPAR::STR::PreStress::mulf)
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_pyramid5::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOP5, NUMDIM_SOP5>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (pstype_ == INPAR::STR::PreStress::mulf)
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
    {
      int numgpt = NUMGPT_SOP5;
      // see whether I am actually a So_pyramid5fbar element
      DRT::ELEMENTS::So_pyramid5fbar* me = dynamic_cast<DRT::ELEMENTS::So_pyramid5fbar*>(this);
      if (me) numgpt += 1;  // one more history entry for centroid data in pyramid5fbar
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOP5, numgpt));
    }
    prestress_->Unpack(tmpprestress);
    // end
  }

  // invdesign_
  else if (pstype_ == INPAR::STR::PreStress::id)
  {
    std::vector<char> tmpinvdesign(0);
    ExtractfromPack(position, data, tmpinvdesign);
    if (invdesign_ == Teuchos::null)
      invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOP5, NUMGPT_SOP5));
    invdesign_->Unpack(tmpinvdesign);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_pyramid5::~So_pyramid5() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_pyramid5::Print(std::ostream& os) const
{
  os << "So_pyramid5 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes     seitz 03/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_pyramid5::sop5_expol(
    LINALG::Matrix<NUMGPT_SOP5, MAT::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  // static variables, that are the same for every element
  static LINALG::Matrix<NUMNOD_SOP5, NUMGPT_SOP5> expol;
  static bool isfilled;

  if (isfilled == false)
  {
    // the 8 gauss-points create a HEX-element inside the pyramid.
    // expol holds the shapefunction-values of the HEX-element, evaluated at the pyramid-nodes.
    expol(0, 0) = 2.408235313815748;
    expol(0, 1) = -0.6452847075210328;
    expol(0, 2) = 0.1729035162684118;
    expol(0, 3) = -0.6452847075210328;
    expol(0, 4) = -0.542209910031327;
    expol(0, 5) = 0.1452847075210439;
    expol(0, 6) = -0.03892892005285509;
    expol(0, 7) = 0.1452847075210439;
    expol(1, 0) = -0.6452847075210328;
    expol(1, 1) = 2.408235313815748;
    expol(1, 2) = -0.6452847075210328;
    expol(1, 3) = 0.1729035162684118;
    expol(1, 4) = 0.1452847075210439;
    expol(1, 5) = -0.542209910031327;
    expol(1, 6) = 0.1452847075210439;
    expol(1, 7) = -0.03892892005285509;
    expol(2, 0) = 0.1729035162684118;
    expol(2, 1) = -0.6452847075210328;
    expol(2, 2) = 2.408235313815748;
    expol(2, 3) = -0.6452847075210328;
    expol(2, 4) = -0.03892892005285509;
    expol(2, 5) = 0.1452847075210439;
    expol(2, 6) = -0.542209910031327;
    expol(2, 7) = 0.1452847075210439;
    expol(3, 0) = -0.6452847075210328;
    expol(3, 1) = 0.1729035162684118;
    expol(3, 2) = -0.6452847075210328;
    expol(3, 3) = 2.408235313815748;
    expol(3, 4) = 0.1452847075210439;
    expol(3, 5) = -0.03892892005285509;
    expol(3, 6) = 0.1452847075210439;
    expol(3, 7) = -0.542209910031327;
    expol(4, 0) = -0.2702847075210531;
    expol(4, 1) = -0.2702847075210531;
    expol(4, 2) = -0.2702847075210531;
    expol(4, 3) = -0.2702847075210531;
    expol(4, 4) = 0.520284707521053;
    expol(4, 5) = 0.520284707521053;
    expol(4, 6) = 0.520284707521053;
    expol(4, 7) = 0.520284707521053;

    isfilled = true;
  }

  LINALG::Matrix<NUMNOD_SOP5, MAT::NUM_STRESS_3D> nodalstresses;
  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_SOP5; ++i)
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
/* 5-node pyramid node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (r,s,t) of nodes
 * of biunit pyramid [-1,1]x[-1,1]x[0,1]
 * 5-node pyramid: node 1,2,3,4,5


 *                /(5)\
 *              / //\\ \
 *            // //  \\ \\
 *          //  //    \\  \\
 *        //   //  t   \\   \\
 *      //    //   |    \\    \\
 *    //     //    |     \\     \\
 *  (4)-----//------------\\-----(3)
 *  ||     //      |       \\     ||
 *  ||    //       |        \\    ||
 *  ||   //        o---------\\---------s
 *  ||  //         |          \\  ||
 *  || //          |           \\ ||
 *  ||//           |            \\||
 *  (1)============|=============(2)
 *                 |
 *                 r
 *
 */
/*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                           |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_pyramid5::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_pyramid5::Surfaces()
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
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_pyramid5::Lines()
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
void DRT::ELEMENTS::So_pyramid5::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);
  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                                  |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_pyramid5::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOP5, this->Id());
}
