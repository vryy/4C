/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "so_hex8fbar.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"


DRT::ELEMENTS::So_hex8Type DRT::ELEMENTS::So_hex8Type::instance_;

DRT::ELEMENTS::So_hex8Type& DRT::ELEMENTS::So_hex8Type::Instance() { return instance_; }

namespace
{
  const std::string name = DRT::ELEMENTS::So_hex8Type::Instance().Name();
}

DRT::ParObject* DRT::ELEMENTS::So_hex8Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_hex8* object = new DRT::ELEMENTS::So_hex8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex8(id, owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex8(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_hex8Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_hex8Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::So_hex8Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8"];

  defs["HEX8"]
      .AddIntVector("HEX8", 8)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
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

// initialization of static gauss point rule for the so_hex8 element
const DRT::UTILS::IntPointsAndWeights<NUMDIM_SOH8> DRT::ELEMENTS::So_hex8::gp_rule_(
    DRT::UTILS::IntPointsAndWeights<NUMDIM_SOH8>(
        static_cast<const enum DRT::UTILS::GaussRule3D>(GP_RULE_SOH8::rule)));

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(int id, int owner)
    : So_base(id, owner),
      data_(),
      analyticalmaterialtangent_(true),
      pstype_(INPAR::STR::PreStress::none),
      pstime_(0.0),
      time_(0.0),
      old_step_length_(0.0)
{
  eastype_ = soh8_easnone;
  neas_ = 0;
  invJ_.resize(NUMGPT_SOH8, LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>(true));
  detJ_.resize(NUMGPT_SOH8, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
    if (DRT::INPUT::IntegralValue<int>(sdyn, "MATERIALTANGENT")) analyticalmaterialtangent_ = false;
  }
  if (pstype_ == INPAR::STR::PreStress::mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8, NUMGPT_SOH8));

  if (pstype_ == INPAR::STR::PreStress::id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOH8, NUMGPT_SOH8));

  if (DRT::Problem::Instance()->GetProblemType() == prb_struct_ale)
  {
    if (kintype_ == INPAR::STR::kinem_linear)
      dserror("Structure-Ale approach only for nonlinear kinematics !!!");

    structale_ = true;
  }
  else
    structale_ = false;


  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(const DRT::ELEMENTS::So_hex8& old)
    : So_base(old),
      eastype_(old.eastype_),
      neas_(old.neas_),
      data_(old.data_),
      detJ_(old.detJ_),
      analyticalmaterialtangent_(old.analyticalmaterialtangent_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_),
      old_step_length_(old.old_step_length_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH8 x NUMDIM_SOH8?
    // invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
    invJ_[i] = old.invJ_[i];
  }

  if (pstype_ == INPAR::STR::PreStress::mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  if (pstype_ == INPAR::STR::PreStress::id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(*(old.invdesign_)));

  if (DRT::Problem::Instance()->GetProblemType() == prb_struct_ale)
  {
    if (kintype_ == INPAR::STR::kinem_linear)
      dserror("Structure-Ale approach only for nonlinear kinematics !!!");

    structale_ = true;
  }
  else
    structale_ = false;

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex8::Clone() const
{
  DRT::ELEMENTS::So_hex8* newelement = new DRT::ELEMENTS::So_hex8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_hex8::Shape() const { return hex8; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);
  // eastype_
  AddtoPack(data, eastype_);
  // neas_
  AddtoPack(data, neas_);
  // analyticalmaterialtangent_
  AddtoPack(data, analyticalmaterialtangent_);
  // data_
  AddtoPack(data, data_);
  // line search
  AddtoPack(data, old_step_length_);
  // prestress_
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (pstype_ == INPAR::STR::PreStress::mulf)
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }
  // invdesign_
  else if (pstype_ == INPAR::STR::PreStress::id)
  {
    DRT::ParObject::AddtoPack(data, *invdesign_);
  }

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Unpack(const std::vector<char>& data)
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
  // eastype_
  eastype_ = static_cast<EASType>(ExtractInt(position, data));
  // neas_
  ExtractfromPack(position, data, neas_);
  // analyticalmaterialtangent_
  analyticalmaterialtangent_ = ExtractInt(position, data);
  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);
  // line search
  ExtractfromPack(position, data, old_step_length_);
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
      int numgpt = NUMGPT_SOH8;
      // see whether I am actually a So_hex8fbar element
      DRT::ELEMENTS::So_hex8fbar* me = dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(this);
      if (me) numgpt += 1;  // one more history entry for centroid data in hex8fbar
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8, numgpt));
    }
    prestress_->Unpack(tmpprestress);
  }
  // invdesign_
  else if (pstype_ == INPAR::STR::PreStress::id)
  {
    std::vector<char> tmpinvdesign(0);
    ExtractfromPack(position, data, tmpinvdesign);
    if (invdesign_ == Teuchos::null)
      invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOH8, NUMGPT_SOH8));
    invdesign_->Unpack(tmpinvdesign);
  }

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::~So_hex8() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Print(std::ostream& os) const
{
  os << "So_hex8 ";
  Element::Print(os);
  // std::cout << std::endl;
  // std::cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_expol(
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  // static variables, that are the same for every element
  static LINALG::Matrix<NUMNOD_SOH8, NUMGPT_SOH8> expol;
  static bool isfilled;

  if (isfilled == false)
  {
    double sq3 = sqrt(3.0);
    expol(0, 0) = 1.25 + 0.75 * sq3;
    expol(0, 1) = -0.25 - 0.25 * sq3;
    expol(0, 2) = -0.25 + 0.25 * sq3;
    expol(0, 3) = -0.25 - 0.25 * sq3;
    expol(0, 4) = -0.25 - 0.25 * sq3;
    expol(0, 5) = -0.25 + 0.25 * sq3;
    expol(0, 6) = 1.25 - 0.75 * sq3;
    expol(0, 7) = -0.25 + 0.25 * sq3;
    expol(1, 1) = 1.25 + 0.75 * sq3;
    expol(1, 2) = -0.25 - 0.25 * sq3;
    expol(1, 3) = -0.25 + 0.25 * sq3;
    expol(1, 4) = -0.25 + 0.25 * sq3;
    expol(1, 5) = -0.25 - 0.25 * sq3;
    expol(1, 6) = -0.25 + 0.25 * sq3;
    expol(1, 7) = 1.25 - 0.75 * sq3;
    expol(2, 2) = 1.25 + 0.75 * sq3;
    expol(2, 3) = -0.25 - 0.25 * sq3;
    expol(2, 4) = 1.25 - 0.75 * sq3;
    expol(2, 5) = -0.25 + 0.25 * sq3;
    expol(2, 6) = -0.25 - 0.25 * sq3;
    expol(2, 7) = -0.25 + 0.25 * sq3;
    expol(3, 3) = 1.25 + 0.75 * sq3;
    expol(3, 4) = -0.25 + 0.25 * sq3;
    expol(3, 5) = 1.25 - 0.75 * sq3;
    expol(3, 6) = -0.25 + 0.25 * sq3;
    expol(3, 7) = -0.25 - 0.25 * sq3;
    expol(4, 4) = 1.25 + 0.75 * sq3;
    expol(4, 5) = -0.25 - 0.25 * sq3;
    expol(4, 6) = -0.25 + 0.25 * sq3;
    expol(4, 7) = -0.25 - 0.25 * sq3;
    expol(5, 5) = 1.25 + 0.75 * sq3;
    expol(5, 6) = -0.25 - 0.25 * sq3;
    expol(5, 7) = -0.25 + 0.25 * sq3;
    expol(6, 6) = 1.25 + 0.75 * sq3;
    expol(6, 7) = -0.25 - 0.25 * sq3;
    expol(7, 7) = 1.25 + 0.75 * sq3;

    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        expol(i, j) = expol(j, i);
      }
    }
    isfilled = true;
  }

  LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D> nodalstresses;
  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_SOH8; ++i)
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
/* 8-node hexhedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (r,s,t) of nodes
 * of biunit cube [-1,1]x[-1,1]x[-1,1]
 *  8-node hexahedron: node 0,1,...,7
 *                      t
 *                      |
 *             4========|================7
 *           //|        |               /||
 *          // |        |              //||
 *         //  |        |             // ||
 *        //   |        |            //  ||
 *       //    |        |           //   ||
 *      //     |        |          //    ||
 *     //      |        |         //     ||
 *     5=========================6       ||
 *    ||       |        |        ||      ||
 *    ||       |        o--------||---------s
 *    ||       |       /         ||      ||
 *    ||       0------/----------||------3
 *    ||      /      /           ||     //
 *    ||     /      /            ||    //
 *    ||    /      /             ||   //
 *    ||   /      /              ||  //
 *    ||  /      /               || //
 *    || /      r                ||//
 *    ||/                        ||/
 *     1=========================2
 *
 */
/*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex8::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             maf 04/07|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex8::Surfaces()
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
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So_hex8::Lines()
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
std::vector<double> DRT::ELEMENTS::So_hex8::ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOH8, 1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  LINALG::Matrix<1, NUMDIM_SOH8> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}
/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOH8, this->Id());
}
