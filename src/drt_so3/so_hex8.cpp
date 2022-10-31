/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element

\level 1


*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include <Epetra_SerialDenseMatrix.h>
#include "so_element_service.H"
#include "so_hex8fbar.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/prestress_service.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_utils_nullspace.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

#include "../drt_fiber/drt_fiber_node.H"
#include "../drt_fiber/drt_fiber_utils.H"
#include "so_utils.H"
#include "../drt_fiber/nodal_fiber_holder.H"


DRT::ELEMENTS::So_hex8Type DRT::ELEMENTS::So_hex8Type::instance_;

DRT::ELEMENTS::So_hex8Type& DRT::ELEMENTS::So_hex8Type::Instance() { return instance_; }

namespace
{
  const std::string name = DRT::ELEMENTS::So_hex8Type::Instance().Name();
}

DRT::ParObject* DRT::ELEMENTS::So_hex8Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So_hex8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
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

Teuchos::SerialDenseMatrix<int, double> DRT::ELEMENTS::So_hex8Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return LINALG::ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::So_hex8Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

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
        static_cast<enum DRT::UTILS::GaussRule3D>(GP_RULE_SOH8::rule)));

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

    pstype_ = ::UTILS::PRESTRESS::GetType();
    pstime_ = ::UTILS::PRESTRESS::GetPrestressTime();
    if (DRT::INPUT::IntegralValue<int>(sdyn, "MATERIALTANGENT")) analyticalmaterialtangent_ = false;
  }
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8, NUMGPT_SOH8));

  if (::UTILS::PRESTRESS::IsInverseDesign(pstype_))
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

  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  if (::UTILS::PRESTRESS::IsInverseDesign(pstype_))
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
  auto* newelement = new DRT::ELEMENTS::So_hex8(*this);
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
  // Pack prestress type
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }
  // invdesign_
  else if (::UTILS::PRESTRESS::IsInverseDesign(pstype_))
  {
    DRT::ParObject::AddtoPack(data, *invdesign_);
  }

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
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
  // Extract prestress
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
    {
      int numgpt = NUMGPT_SOH8;
      // see whether I am actually a So_hex8fbar element
      auto* me = dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(this);
      if (me) numgpt += 1;  // one more history entry for centroid data in hex8fbar
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8, numgpt));
    }
    prestress_->Unpack(tmpprestress);
  }
  // invdesign_
  else if (::UTILS::PRESTRESS::IsInverseDesign(pstype_))
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
  LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D> nodalstresses;
  static LINALG::Matrix<NUMNOD_SOH8, NUMGPT_SOH8> extrapolationMatrix(
      GaussPointsToNodesExtrapolation());
  nodalstresses.Multiply(extrapolationMatrix, stresses);

  DRT::ELEMENTS::AssembleExtrapolatedNodalValues<LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>>(
      expolstresses, nodalstresses, this);
}

void DRT::ELEMENTS::So_hex8::ExtrapolateGPQuantityToNodesAndAssemble(
    const LINALG::SerialDenseMatrix& gp_quantity, Epetra_MultiVector& global_quantity,
    bool nodal_average)
{
  LINALG::SerialDenseMatrix nodal_quantity(NUMNOD_SOH8, gp_quantity.N());
  // static const LINALG::SerialDenseMatrix extrapolation_matrix(Epetra_DataAccess::View,
  //    GaussPointsToNodesExtrapolation().A(), NUMNOD_SOH8, NUMNOD_SOH8, NUMGPT_SOH8);

  nodal_quantity.Multiply('N', 'N', 1.0, GaussPointsToNodesExtrapolation(), gp_quantity, 0.0);

  DRT::ELEMENTS::AssembleExtrapolatedNodalValues<LINALG::SerialDenseMatrix>(
      global_quantity, nodal_quantity, this, nodal_average);
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

// Compute nodal fibers and call post setup routine of the materials
void DRT::ELEMENTS::So_hex8::MaterialPostSetup(Teuchos::ParameterList& params)
{
  if (DRT::FIBER::UTILS::HaveNodalFibers<hex8>(Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();

    // add fibers to the ParameterList
    // ParameterList does not allow to store a std::vector, so we have to add every gp fiber
    // with a separate key. To keep it clean, It is added to a sublist.
    DRT::FIBER::NodalFiberHolder fiberHolder;

    // Do the interpolation
    DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::hex8>(
        Nodes(), shapefcts, fiberHolder);

    params.set("fiberholder", fiberHolder);
  }

  // Call super post setup
  So_base::MaterialPostSetup(params);

  // Cleanup ParameterList to not carry all fibers the whole simulation
  // do not throw an error if key does not exist.
  params.remove("fiberholder", false);
}