/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief 3D nonlinear Reissner beam element

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_beam3_reissner.hpp"

#include "4C_beam3_spatial_discretization_utils.hpp"
#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_beaminteraction_periodic_boundingbox.hpp"
#include "4C_discretization_fem_general_largerotations.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3rType DRT::ELEMENTS::Beam3rType::instance_;

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3rType& DRT::ELEMENTS::Beam3rType::Instance() { return instance_; }

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::Beam3rType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Beam3r* object = new DRT::ELEMENTS::Beam3r(-1, -1);
  object->Unpack(data);
  return object;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3rType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "BEAM3R")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3rType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(id, owner));
  return ele;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3rType::nodal_block_information(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  DRT::ELEMENTS::Beam3r* currele = dynamic_cast<DRT::ELEMENTS::Beam3r*>(dwele);
  if (!currele) FOUR_C_THROW("cast to Beam3r* failed");

  if (currele->hermite_centerline_interpolation() or currele->num_node() > 2)
  {
    FOUR_C_THROW(
        "method nodal_block_information not implemented for element type beam3r in case of Hermite "
        "interpolation or higher order Lagrange interpolation!");
  }
  else
  {
    numdf = 6;
    dimns = 6;
    nv = 6;
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Beam3rType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented for element type beam3r!");
  return nullspace;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3rType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["BEAM3R"];

  // note: LINE2 refers to linear Lagrange interpolation of centerline AND triad field
  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .add_named_double_vector("TRIADS", 6)
                      .AddOptionalTag("FAD")
                      .Build();

  // note: LINE3 refers to quadratic Lagrange interpolation of centerline AND triad field
  defs["LINE3"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE3", 3)
                      .AddNamedInt("MAT")
                      .add_named_double_vector("TRIADS", 9)
                      .AddOptionalTag("FAD")
                      .Build();

  // note: LINE4 refers to cubic Lagrange interpolation of centerline AND triad field
  defs["LINE4"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE4", 4)
                      .AddNamedInt("MAT")
                      .add_named_double_vector("TRIADS", 12)
                      .AddOptionalTag("FAD")
                      .Build();

  // note: LINE5 refers to quartic Lagrange interpolation of centerline AND triad field
  defs["LINE5"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE5", 5)
                      .AddNamedInt("MAT")
                      .add_named_double_vector("TRIADS", 15)
                      .AddOptionalTag("FAD")
                      .Build();

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LINE2 refers to linear Lagrange interpolation of the triad field*/
  defs["HERM2LINE2"] = INPUT::LineDefinition::Builder()
                           .AddIntVector("HERM2LINE2", 2)
                           .AddNamedInt("MAT")
                           .add_named_double_vector("TRIADS", 6)
                           .AddOptionalTag("FAD")
                           .Build();

  /* note: HERM2 refers to cubic order Hermite interpolation of centerline (2 nodes)
   *       LINE3 refers to quadratic Lagrange interpolation of the triad field*/
  defs["HERM2LINE3"] = INPUT::LineDefinition::Builder()
                           .AddIntVector("HERM2LINE3", 3)
                           .AddNamedInt("MAT")
                           .add_named_double_vector("TRIADS", 9)
                           .AddOptionalTag("FAD")
                           .Build();

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LINE4 refers to cubic Lagrange interpolation of the triad field*/
  defs["HERM2LINE4"] = INPUT::LineDefinition::Builder()
                           .AddIntVector("HERM2LINE4", 4)
                           .AddNamedInt("MAT")
                           .add_named_double_vector("TRIADS", 12)
                           .AddOptionalTag("FAD")
                           .Build();

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LINE5 refers to quartic Lagrange interpolation of the triad field*/
  defs["HERM2LINE5"] = INPUT::LineDefinition::Builder()
                           .AddIntVector("HERM2LINE5", 5)
                           .AddNamedInt("MAT")
                           .add_named_double_vector("TRIADS", 15)
                           .AddOptionalTag("FAD")
                           .Build();
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3rType::Initialize(DRT::Discretization& dis)
{
  // setting up geometric variables for beam3r elements
  for (int num = 0; num < dis.NumMyColElements(); ++num)
  {
    /* in case that current element is not a beam3r element there is nothing to do and we go back
     * to the head of the loop*/
    if (dis.lColElement(num)->ElementType() != *this) continue;

    // if we get so far current element is a beam3r element and we get a pointer at it
    DRT::ELEMENTS::Beam3r* currele = dynamic_cast<DRT::ELEMENTS::Beam3r*>(dis.lColElement(num));
    if (!currele) FOUR_C_THROW("cast to Beam3r* failed");

    // reference node position
    std::vector<double> xrefe;
    std::vector<double> rotrefe;

    /* the triad field is discretized with Lagrange polynomials of order num_node()-1;
     * the centerline is either discretized in the same way or with 3rd order Hermite polynomials;
     * in case of Hermite interpolation of the centerline, always the two boundary nodes are used
     * for centerline interpolation*/
    const bool centerline_hermite = currele->hermite_centerline_interpolation();

    // nnodetriad: number of nodes used for interpolation of triad field
    // nnodecl: number of nodes used for interpolation of centerline
    // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
    const int nnodetriad = currele->num_node();
    int nnodecl = nnodetriad;
    if (centerline_hermite) nnodecl = 2;

    // resize xrefe and rotrefe for the number of (external) DOFs we need to store
    xrefe.resize(3 * nnodecl);
    rotrefe.resize(3 * nnodetriad);

    // getting element's nodal coordinates and treating them as reference configuration
    /* note: in case of Hermite interpolation of centerline, the reference config of tangent DOFs
     *       is computed from the reference triads, i.e. rotrefe*/
    for (int node = 0; node < nnodetriad; node++)
      for (int dim = 0; dim < 3; dim++)
        rotrefe[node * 3 + dim] = currele->InitialNodalRotVecs()[node](dim);

    // the next section is needed in case of periodic boundary conditions and a shifted
    // configuration (i.e. elements cut by the periodic boundary) in the input file
    Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> periodic_boundingbox =
        Teuchos::rcp(new CORE::GEO::MESHFREE::BoundingBox());
    periodic_boundingbox->Init();  // no Setup() call needed here

    std::vector<double> disp_shift;
    int numdof = currele->NumDofPerNode(*(currele->Nodes()[0]));
    disp_shift.resize(numdof * nnodecl);
    for (unsigned int i = 0; i < disp_shift.size(); ++i) disp_shift[i] = 0.0;
    if (periodic_boundingbox->HavePBC())
      currele->UnShiftNodePosition(disp_shift, *periodic_boundingbox);

    for (int node = 0; node < nnodecl; ++node)
    {
      if (currele->Nodes()[node] == nullptr)
        FOUR_C_THROW("beam3r: Cannot get nodes in order to compute reference configuration");

      for (unsigned int dim = 0; dim < 3; ++dim)
        xrefe[node * 3 + dim] = currele->Nodes()[node]->X()[dim] + disp_shift[node * numdof + dim];
    }

    // set_up_reference_geometry is a templated function
    switch (nnodetriad)
    {
      case 2:
      {
        if (!centerline_hermite)
          currele->set_up_reference_geometry<2, 2, 1>(xrefe, rotrefe);
        else
          currele->set_up_reference_geometry<2, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 3:
      {
        if (!centerline_hermite)
          currele->set_up_reference_geometry<3, 3, 1>(xrefe, rotrefe);
        else
          currele->set_up_reference_geometry<3, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 4:
      {
        if (!centerline_hermite)
          currele->set_up_reference_geometry<4, 4, 1>(xrefe, rotrefe);
        else
          currele->set_up_reference_geometry<4, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 5:
      {
        if (!centerline_hermite)
          currele->set_up_reference_geometry<5, 5, 1>(xrefe, rotrefe);
        else
          currele->set_up_reference_geometry<5, 2, 2>(xrefe, rotrefe);
        break;
      }
      default:
        FOUR_C_THROW("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        break;
    }
  }

  return 0;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3r::Beam3r(int id, int owner)
    : DRT::ELEMENTS::Beam3Base(id, owner),
      stiff_ptc_(),
      use_fad_(false),
      isinit_(false),
      jacobi_gp_elastf_(0),
      jacobi_gp_elastm_(0),
      jacobi_gp_mass_(0),
      jacobi_gp_dampstoch_(0),
      jacobi_gp_neumannline_(0),
      eint_(0.0),
      ekin_(0.0),
      ekintorsion_(0.0),
      ekinbending_(0.0),
      ekintrans_(0.0),
      l_(true),
      p_(true)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3r::Beam3r(const DRT::ELEMENTS::Beam3r& old)
    : DRT::ELEMENTS::Beam3Base(old),
      use_fad_(old.use_fad_),
      isinit_(old.isinit_),
      reflength_(old.reflength_),
      theta0node_(old.theta0node_),
      tcurrnode_(old.tcurrnode_),
      kref_gp_(old.kref_gp_),
      gammaref_gp_(old.gammaref_gp_),
      jacobi_gp_elastf_(old.jacobi_gp_elastf_),
      jacobi_gp_elastm_(old.jacobi_gp_elastm_),
      jacobi_gp_mass_(old.jacobi_gp_mass_),
      jacobi_gp_dampstoch_(old.jacobi_gp_dampstoch_),
      jacobi_gp_neumannline_(old.jacobi_gp_neumannline_),
      qconvnode_(old.qconvnode_),
      qnewnode_(old.qnewnode_),
      qconv_gp_mass_(old.qconv_gp_mass_),
      qnew_gp_mass_(old.qnew_gp_mass_),
      wconv_gp_mass_(old.wconv_gp_mass_),
      wnew_gp_mass_(old.wnew_gp_mass_),
      aconv_gp_mass_(old.aconv_gp_mass_),
      anew_gp_mass_(old.anew_gp_mass_),
      amodconv_gp_mass_(old.amodconv_gp_mass_),
      amodnew_gp_mass_(old.amodnew_gp_mass_),
      rttconv_gp_mass_(old.rttconv_gp_mass_),
      rttnew_gp_mass_(old.rttnew_gp_mass_),
      rttmodconv_gp_mass_(old.rttmodconv_gp_mass_),
      rttmodnew_gp_mass_(old.rttmodnew_gp_mass_),
      rtconv_gp_mass_(old.rtconv_gp_mass_),
      rtnew_gp_mass_(old.rtnew_gp_mass_),
      rconv_gp_mass_(old.rconv_gp_mass_),
      rnew_gp_mass_(old.rnew_gp_mass_),
      qconv_gp_dampstoch_(old.qconv_gp_dampstoch_),
      qnew_gp_dampstoch_(old.qnew_gp_dampstoch_),
      eint_(old.eint_),
      ekin_(old.ekin_),
      ekintorsion_(old.ekintorsion_),
      ekinbending_(old.ekinbending_),
      ekintrans_(old.ekintrans_),
      l_(old.l_),
      p_(old.p_),
      axial_strain_gp_elastf_(old.axial_strain_gp_elastf_),
      shear_strain_2_gp_elastf_(old.shear_strain_2_gp_elastf_),
      shear_strain_3_gp_elastf_(old.shear_strain_3_gp_elastf_),
      twist_gp_elastm_(old.twist_gp_elastm_),
      curvature_2_gp_elastm_(old.curvature_2_gp_elastm_),
      curvature_3_gp_elastm_(old.curvature_3_gp_elastm_),
      material_axial_force_gp_elastf_(old.material_axial_force_gp_elastf_),
      material_shear_force_2_gp_elastf_(old.material_shear_force_2_gp_elastf_),
      material_shear_force_3_gp_elastf_(old.material_shear_force_3_gp_elastf_),
      material_torque_gp_elastm_(old.material_torque_gp_elastm_),
      material_bending_moment_2_gp_elastm_(old.material_bending_moment_2_gp_elastm_),
      material_bending_moment_3_gp_elastm_(old.material_bending_moment_3_gp_elastm_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3r and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3r::Clone() const
{
  DRT::ELEMENTS::Beam3r* newelement = new DRT::ELEMENTS::Beam3r(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::Print(std::ostream& os) const
{
  os << "beam3r ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Beam3r::Shape() const
{
  int numnodes = num_node();
  switch (numnodes)
  {
    case 2:
      return CORE::FE::CellType::line2;
      break;
    case 3:
      return CORE::FE::CellType::line3;
      break;
    case 4:
      return CORE::FE::CellType::line4;
      break;
    case 5:
      return CORE::FE::CellType::line5;
      break;
    default:
      FOUR_C_THROW("Only Line2, Line3, Line4 and Line5 elements are implemented.");
      break;
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Beam3Base::Pack(data);

  // add all class variables of beam3r element
  AddtoPack(data, use_fad_);
  AddtoPack(data, centerline_hermite_);
  AddtoPack(data, isinit_);
  AddtoPack(data, reflength_);
  AddtoPack<3, 1>(data, theta0node_);
  AddtoPack<3, 1>(data, Tref_);
  AddtoPack<3, 1>(data, tcurrnode_);
  AddtoPack<3, 1>(data, kref_gp_);
  AddtoPack<3, 1>(data, gammaref_gp_);
  AddtoPack(data, jacobi_gp_elastf_);
  AddtoPack(data, jacobi_gp_elastm_);
  AddtoPack(data, jacobi_gp_mass_);
  AddtoPack(data, jacobi_gp_dampstoch_);
  AddtoPack(data, jacobi_gp_neumannline_);
  AddtoPack<4, 1>(data, qconvnode_);
  AddtoPack<4, 1>(data, qnewnode_);
  AddtoPack<4, 1>(data, qconv_gp_mass_);
  AddtoPack<4, 1>(data, qnew_gp_mass_);
  AddtoPack<3, 1>(data, wconv_gp_mass_);
  AddtoPack<3, 1>(data, wnew_gp_mass_);
  AddtoPack<3, 1>(data, aconv_gp_mass_);
  AddtoPack<3, 1>(data, anew_gp_mass_);
  AddtoPack<3, 1>(data, amodnew_gp_mass_);
  AddtoPack<3, 1>(data, amodconv_gp_mass_);
  AddtoPack<3, 1>(data, rttconv_gp_mass_);
  AddtoPack<3, 1>(data, rttnew_gp_mass_);
  AddtoPack<3, 1>(data, rttmodconv_gp_mass_);
  AddtoPack<3, 1>(data, rttmodnew_gp_mass_);
  AddtoPack<3, 1>(data, rtconv_gp_mass_);
  AddtoPack<3, 1>(data, rtnew_gp_mass_);
  AddtoPack<3, 1>(data, rconv_gp_mass_);
  AddtoPack<3, 1>(data, rnew_gp_mass_);
  AddtoPack<4, 1>(data, qconv_gp_dampstoch_);
  AddtoPack<4, 1>(data, qnew_gp_dampstoch_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Beam3Base::Unpack(basedata);

  // extract all class variables of beam3r element
  use_fad_ = ExtractInt(position, data);
  centerline_hermite_ = ExtractInt(position, data);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack(position, data, reflength_);
  ExtractfromPack<3, 1>(position, data, theta0node_);
  ExtractfromPack<3, 1>(position, data, Tref_);
  ExtractfromPack<3, 1>(position, data, tcurrnode_);
  ExtractfromPack<3, 1>(position, data, kref_gp_);
  ExtractfromPack<3, 1>(position, data, gammaref_gp_);
  ExtractfromPack(position, data, jacobi_gp_elastf_);
  ExtractfromPack(position, data, jacobi_gp_elastm_);
  ExtractfromPack(position, data, jacobi_gp_mass_);
  ExtractfromPack(position, data, jacobi_gp_dampstoch_);
  ExtractfromPack(position, data, jacobi_gp_neumannline_);
  ExtractfromPack<4, 1>(position, data, qconvnode_);
  ExtractfromPack<4, 1>(position, data, qnewnode_);
  ExtractfromPack<4, 1>(position, data, qconv_gp_mass_);
  ExtractfromPack<4, 1>(position, data, qnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, wconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, wnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, aconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, anew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, amodconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, amodnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rttconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rttnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rttmodconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rttmodnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rtconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rtnew_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rconv_gp_mass_);
  ExtractfromPack<3, 1>(position, data, rnew_gp_mass_);
  ExtractfromPack<4, 1>(position, data, qconv_gp_dampstoch_);
  ExtractfromPack<4, 1>(position, data, qnew_gp_dampstoch_);

  // NOT communicated
  eint_ = 0.0;
  ekin_ = 0.0;
  ekintorsion_ = 0.0;
  ekinbending_ = 0.0;
  ekintrans_ = 0.0;
  l_.Clear();
  p_.Clear();
  kmax_ = 0.0;
  axial_strain_gp_elastf_.clear();
  shear_strain_2_gp_elastf_.clear();
  shear_strain_3_gp_elastf_.clear();
  twist_gp_elastm_.clear();
  curvature_2_gp_elastm_.clear();
  curvature_3_gp_elastm_.clear();
  material_axial_force_gp_elastf_.clear();
  material_shear_force_2_gp_elastf_.clear();
  material_shear_force_3_gp_elastf_.clear();
  material_torque_gp_elastm_.clear();
  material_bending_moment_2_gp_elastm_.clear();
  material_bending_moment_3_gp_elastm_.clear();
  spatial_x_force_gp_elastf_.clear();
  spatial_y_force_2_gp_elastf_.clear();
  spatial_z_force_3_gp_elastf_.clear();
  spatial_x_moment_gp_elastm_.clear();
  spatial_y_moment_2_gp_elastm_.clear();
  spatial_z_moment_3_gp_elastm_.clear();

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Beam3r::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------*
 | determine Gauss rule from purpose and interpolation scheme grill 03/16|
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule1D DRT::ELEMENTS::Beam3r::MyGaussRule(const IntegrationPurpose intpurpose) const
{
  const CORE::FE::CellType distype = this->Shape();

  switch (intpurpose)
  {
    // anti-locking: reduced integration of elastic residual contributions from forces (-> 'Gamma
    // terms')
    case res_elastic_force:
    {
      switch (distype)
      {
        case CORE::FE::CellType::line2:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_1point;
          else
            return CORE::FE::GaussRule1D::line_lobatto3point;
        }
        case CORE::FE::CellType::line3:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_2point;
          else
            return CORE::FE::GaussRule1D::line_lobatto3point;
        }
        case CORE::FE::CellType::line4:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_3point;
          else
            return CORE::FE::GaussRule1D::line_lobatto3point;
        }
        case CORE::FE::CellType::line5:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_4point;
          else
            return CORE::FE::GaussRule1D::line_lobatto3point;
        }
        default:
        {
          FOUR_C_THROW("unknown discretization type!");
          break;
        }
      }
      break;
    }

    /* reduced integration of elastic residual contributions from moments (-> 'curvature terms')
     * NOT required for anti-locking, but for historic reasons we keep this in case of Lagrange
     * interpolation of centerline 'full' integration in case of Hermite centerline interpolation */
    case res_elastic_moment:
    {
      switch (distype)
      {
        case CORE::FE::CellType::line2:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_1point;
          else
            return CORE::FE::GaussRule1D::line_2point;
        }
        case CORE::FE::CellType::line3:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_2point;
          else
            return CORE::FE::GaussRule1D::line_3point;
        }
        case CORE::FE::CellType::line4:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_3point;
          else
            return CORE::FE::GaussRule1D::line_4point;
        }
        case CORE::FE::CellType::line5:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_4point;
          else
            return CORE::FE::GaussRule1D::line_5point;
        }
        default:
        {
          FOUR_C_THROW("unknown discretization type!");
          break;
        }
      }
      break;
    }

    // 'full' integration of inertia contributions
    case res_inertia:
    {
      switch (distype)
      {
        case CORE::FE::CellType::line2:
        {
          return CORE::FE::GaussRule1D::line_2point;
        }
        case CORE::FE::CellType::line3:
        {
          return CORE::FE::GaussRule1D::line_3point;
        }
        case CORE::FE::CellType::line4:
        {
          return CORE::FE::GaussRule1D::line_4point;
        }
        case CORE::FE::CellType::line5:
        {
          return CORE::FE::GaussRule1D::line_5point;
        }
        default:
        {
          FOUR_C_THROW("unknown discretization type!");
          break;
        }
      }
      break;
    }

    // 'full' integration of damping and stochastic contributions
    case res_damp_stoch:
    {
      return CORE::FE::GaussRule1D::line_4point;
    }

    /* 'full' integration of Neumann line loads
     * higher order Gauss quadrature scheme may prove useful in case of abnormal convergence
     * behaviour due to 'complex' line loads*/
    case neumann_lineload:
    {
      switch (distype)
      {
        case CORE::FE::CellType::line2:
        {
          if (!centerline_hermite_)
            return CORE::FE::GaussRule1D::line_1point;
          else
            return CORE::FE::GaussRule1D::line_2point;
        }
        case CORE::FE::CellType::line3:
        {
          return CORE::FE::GaussRule1D::line_2point;
        }
        case CORE::FE::CellType::line4:
        {
          return CORE::FE::GaussRule1D::line_3point;
        }
        case CORE::FE::CellType::line5:
        {
          return CORE::FE::GaussRule1D::line_4point;
        }
        default:
        {
          FOUR_C_THROW("unknown discretization type!");
          break;
        }
      }
      break;
    }

    default:
    {
      FOUR_C_THROW("beam3r: unknown purpose for numerical quadrature!");
      break;
    }
  }

  return CORE::FE::GaussRule1D::undefined;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::set_up_reference_geometry(
    const std::vector<double>& xrefe, const std::vector<double>& rotrefe)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  /* in case of Hermite interpolation of the centerline, always the two boundary nodes (ID0 & ID1)
   * are used for centerline interpolation; the triad field may be interpolated with Lagrange
   * polynomials of order 1-4 (linear-quartic), i.e. nnodetriad=2...5*/
  if (centerline_hermite_ and nnodecl != 2)
    FOUR_C_THROW("Only 3rd order Hermite interpolation of beam centerline implemented!");

  /* this method initializes geometric variables of the element; the initialization can usually be
   * applied to elements only once; therefore after the first initialization the flag isinit_ is set
   * to true and from then on this method does not take any action when called again unless it is
   * called on purpose with the additional parameter secondinit. If this parameter is passed into
   * the method and is true the element is initialized another time with xrefe;
   * note: the isinit_ flag is important for avoiding re-initialization upon restart. However, it
   * should be possible to conduct a
   * second initialization in principle (e.g. for periodic boundary conditions*/

  if (!isinit_)
  {
    isinit_ = true;

    // check input data
    if (xrefe.size() != 3 * nnodecl)
      FOUR_C_THROW(
          "size mismatch in given position vector for stress-free reference geometry of beam3r:"
          " expected %d and got %d entries!",
          3 * nnodecl, xrefe.size());

    if (rotrefe.size() != 3 * nnodetriad)
      FOUR_C_THROW(
          "size mismatch in given rotation vector for stress-free reference geometry of beam3r:"
          " expected %d and got %d entries!",
          3 * nnodetriad, rotrefe.size());



    /********************************** Initialize/resize general variables
     ********************************
     *****************************************************************************************************/

    // create object of triad interpolation scheme
    Teuchos::RCP<LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, double>>
        triad_interpolation_scheme_ptr = Teuchos::rcp(
            new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, double>());

    // Get DiscretizationType
    CORE::FE::CellType distype = Shape();

    /* Note: index i refers to the i-th shape function (i = 0 ... nnode*vpernode-1)
     * the vectors store individual shape functions, NOT an assembled matrix of shape functions) */
    /* vector whose numgp-th element is a 1xnnode-matrix with all Lagrange shape functions evaluated
     * at the numgp-th GP these shape functions are used for the interpolation of the triad field*/
    std::vector<CORE::LINALG::Matrix<1, nnodetriad, double>> I_i;
    // same for derivatives
    std::vector<CORE::LINALG::Matrix<1, nnodetriad, double>> I_i_xi;

    /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite)
     * shape functions evaluated at the numgp-th GP these shape functions are used for the
     * interpolation of the beam centerline*/
    std::vector<CORE::LINALG::Matrix<1, vpernode * nnodecl, double>> H_i;
    // same for the derivatives
    std::vector<CORE::LINALG::Matrix<1, vpernode * nnodecl, double>> H_i_xi;

    // beside the nodal reference positions from xrefe, this vector also holds the reference
    // tangents in case of Hermite interpolation of the beam centerline
    CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1> pos_ref_centerline;

    // initial curve in physical space and derivative with respect to curve parameter xi \in [-1;1]
    // on element level
    CORE::LINALG::Matrix<3, 1> r0;
    CORE::LINALG::Matrix<3, 1> dr0dxi;

    // dummy 3x1 vector
    CORE::LINALG::Matrix<3, 1> dummy(true);


    /********************************** Compute nodal quantities
     *******************************************
     *****************************************************************************************************/

    /********************* store given nodal triads as quaternions in class variable
     * *********************/
    qnewnode_.resize(nnodetriad);
    qconvnode_.resize(nnodetriad);

    // nodal triads in stress-free configuration
    for (unsigned int node = 0; node < nnodetriad; node++)
    {
      CORE::LINALG::Matrix<3, 1> rotvec(&rotrefe[3 * node]);
      CORE::LARGEROTATIONS::angletoquaternion(rotvec, qnewnode_[node]);
    }

    qconvnode_ = qnewnode_;

    std::vector<CORE::LINALG::Matrix<4, 1, double>> Qnewnode;

    for (unsigned int inode = 0; inode < nnodetriad; ++inode)
      Qnewnode.push_back(CORE::LINALG::Matrix<4, 1, double>(qnewnode_[inode], true));

    // reset triad interpolation with nodal quaternions
    triad_interpolation_scheme_ptr->Reset(Qnewnode);

    CORE::LINALG::Matrix<3, 3> Gref;
    Tref_.resize(nnodecl);

    for (unsigned int node = 0; node < nnodecl; node++)
    {
      /* Calculate the (initial reference triads) = (initial material triads) at the nodes out of
       * the angles theta0node_. So far the initial value for the relative angle is set to zero,
       * i.e. material coordinate system and reference system in the reference configuration
       * coincidence (only at the nodes)*/
      Gref.Clear();
      CORE::LARGEROTATIONS::quaterniontotriad(qnewnode_[node], Gref);
      // store initial nodal tangents in class variable
      for (int i = 0; i < 3; i++) (Tref_[node])(i) = (Gref)(i, 0);

      // fill disp_refe_centerline with reference nodal centerline positions and tangents
      for (int dim = 0; dim < 3; ++dim)
      {
        pos_ref_centerline(3 * vpernode * node + dim) = xrefe[3 * node + dim];
        if (centerline_hermite_)
          pos_ref_centerline(3 * vpernode * node + 3 + dim) = (Tref_[node])(dim);
      }
    }

    reflength_ = CalcReflength<nnodecl, vpernode>(pos_ref_centerline);

    /************************ Compute quantities required for elasticity
     ***********************************
     *****************************************************************************************************/

    // interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11),
    // Jelenic 1999
    CORE::LINALG::Matrix<3, 1> Psi_l;
    /* derivative of interpolated local relative rotation \Psi^l with respect to arc-length
     * parameter at a certain Gauss point according to (3.11), Jelenic 1999*/
    CORE::LINALG::Matrix<3, 1> Psi_l_s;
    // triad at GP
    CORE::LINALG::Matrix<3, 3> Lambda;

    //*********************** preparation for residual contributions from forces
    //***************************

    // Get the applied integration scheme
    CORE::FE::IntegrationPoints1D gausspoints_elast_force(MyGaussRule(res_elastic_force));

    jacobi_gp_elastf_.resize(gausspoints_elast_force.nquad);
    gammaref_gp_.resize(gausspoints_elast_force.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_elast_force.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_elast_force, H_i_xi, distype, this->RefLength());


    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(axial_strain_gp_elastf_.begin(), axial_strain_gp_elastf_.end(), 0.0);
    shear_strain_2_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(shear_strain_2_gp_elastf_.begin(), shear_strain_2_gp_elastf_.end(), 0.0);
    shear_strain_3_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(shear_strain_3_gp_elastf_.begin(), shear_strain_3_gp_elastf_.end(), 0.0);

    material_axial_force_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(material_axial_force_gp_elastf_.begin(), material_axial_force_gp_elastf_.end(), 0.0);
    material_shear_force_2_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(
        material_shear_force_2_gp_elastf_.begin(), material_shear_force_2_gp_elastf_.end(), 0.0);
    material_shear_force_3_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(
        material_shear_force_3_gp_elastf_.begin(), material_shear_force_3_gp_elastf_.end(), 0.0);

    spatial_x_force_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_x_force_gp_elastf_.begin(), spatial_x_force_gp_elastf_.end(), 0.0);
    spatial_y_force_2_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_y_force_2_gp_elastf_.begin(), spatial_y_force_2_gp_elastf_.end(), 0.0);
    spatial_z_force_3_gp_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_z_force_3_gp_elastf_.begin(), spatial_z_force_3_gp_elastf_.end(), 0.0);

    dummy.Clear();

    // Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
    for (int numgp = 0; numgp < gausspoints_elast_force.nquad; ++numgp)
    {
      Calc_r_xi<nnodecl, vpernode, double>(pos_ref_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point for under-integration
      jacobi_gp_elastf_[numgp] = dr0dxi.Norm2();

      // we need dr0ds for computestrain, just reuse dr0dxi from above for simplicity
      dr0dxi.Scale(1.0 / jacobi_gp_elastf_[numgp]);

      triad_interpolation_scheme_ptr->get_interpolated_triad_at_xi(
          Lambda, gausspoints_elast_force.qxg[numgp][0]);

      /* compute material strain Gamma according to Jelenic 1999, eq. (2.12) for reference
       * configuration, i.e. call this function with gammaref=zerovector*/
      compute_gamma<double>(dr0dxi, Lambda, dummy, gammaref_gp_[numgp]);
    }

    //*********************** preparation for residual contributions from moments
    //***************************

    // Get the applied integration scheme
    CORE::FE::IntegrationPoints1D gausspoints_elast_moment(MyGaussRule(res_elastic_moment));

    jacobi_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    kref_gp_.resize(gausspoints_elast_moment.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    I_i.resize(gausspoints_elast_moment.nquad);
    I_i_xi.resize(gausspoints_elast_moment.nquad);
    H_i_xi.resize(gausspoints_elast_moment.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAllGPs<nnodetriad, 1>(
        gausspoints_elast_moment, I_i, I_i_xi, distype);
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_elast_moment, H_i_xi, distype, this->RefLength());

    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    twist_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(twist_gp_elastm_.begin(), twist_gp_elastm_.end(), 0.0);
    curvature_2_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(curvature_2_gp_elastm_.begin(), curvature_2_gp_elastm_.end(), 0.0);
    curvature_3_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(curvature_3_gp_elastm_.begin(), curvature_3_gp_elastm_.end(), 0.0);

    material_torque_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_torque_gp_elastm_.begin(), material_torque_gp_elastm_.end(), 0.0);
    material_bending_moment_2_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_bending_moment_2_gp_elastm_.begin(),
        material_bending_moment_2_gp_elastm_.end(), 0.0);
    material_bending_moment_3_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_bending_moment_3_gp_elastm_.begin(),
        material_bending_moment_3_gp_elastm_.end(), 0.0);

    spatial_x_moment_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_x_moment_gp_elastm_.begin(), spatial_x_moment_gp_elastm_.end(), 0.0);
    spatial_y_moment_2_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_y_moment_2_gp_elastm_.begin(), spatial_y_moment_2_gp_elastm_.end(), 0.0);
    spatial_z_moment_3_gp_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_z_moment_3_gp_elastm_.begin(), spatial_z_moment_3_gp_elastm_.end(), 0.0);


    dummy.Clear();

    // Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
    for (int numgp = 0; numgp < gausspoints_elast_moment.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(pos_ref_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobi_gp_elastm_[numgp] = dr0dxi.Norm2();

      // we need dr0ds for computestrain, just reuse dr0dxi from above for simplicity
      dr0dxi.Scale(1.0 / jacobi_gp_elastm_[numgp]);

      triad_interpolation_scheme_ptr->get_interpolated_local_rotation_vector(Psi_l, I_i[numgp]);

      triad_interpolation_scheme_ptr->get_interpolated_local_rotation_vector_derivative(
          Psi_l_s, I_i_xi[numgp], jacobi_gp_elastm_[numgp]);

      /* compute material curvature K according to Jelenic 1999, eq. (2.12) for reference
       * configuration, i.e. call this function with kapparef=zerovector*/
      compute_k<double>(Psi_l, Psi_l_s, dummy, kref_gp_[numgp]);
    }


    /******************************* Compute quantities required for inertia
     *******************************
     *****************************************************************************************************/

    // Get the applied integration scheme
    CORE::FE::GaussRule1D gaussrule_inertia = MyGaussRule(res_inertia);
    CORE::FE::IntegrationPoints1D gausspoints_inertia(gaussrule_inertia);

    // these quantities will later be used mainly for calculation of inertia terms -> named 'mass'
    jacobi_gp_mass_.resize(gausspoints_inertia.nquad);
    qconv_gp_mass_.resize(gausspoints_inertia.nquad);
    qnew_gp_mass_.resize(gausspoints_inertia.nquad);
    wconv_gp_mass_.resize(gausspoints_inertia.nquad);
    wnew_gp_mass_.resize(gausspoints_inertia.nquad);
    aconv_gp_mass_.resize(gausspoints_inertia.nquad);
    anew_gp_mass_.resize(gausspoints_inertia.nquad);
    rttconv_gp_mass_.resize(gausspoints_inertia.nquad);
    rttnew_gp_mass_.resize(gausspoints_inertia.nquad);
    rttmodconv_gp_mass_.resize(gausspoints_inertia.nquad);
    rttmodnew_gp_mass_.resize(gausspoints_inertia.nquad);
    rtconv_gp_mass_.resize(gausspoints_inertia.nquad);
    rtnew_gp_mass_.resize(gausspoints_inertia.nquad);
    rconv_gp_mass_.resize(gausspoints_inertia.nquad);
    rnew_gp_mass_.resize(gausspoints_inertia.nquad);
    amodconv_gp_mass_.resize(gausspoints_inertia.nquad);
    amodnew_gp_mass_.resize(gausspoints_inertia.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i.resize(gausspoints_inertia.nquad);
    H_i_xi.resize(gausspoints_inertia.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_inertia, H_i, H_i_xi, distype, this->RefLength());

    // Loop through all GPs for exact integration and compute initial jacobi determinant
    for (int numgp = 0; numgp < gausspoints_inertia.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(pos_ref_centerline, H_i_xi[numgp], dr0dxi);
      Calc_r<nnodecl, vpernode, double>(pos_ref_centerline, H_i[numgp], r0);

      // Store Jacobi determinant at this Gauss point
      jacobi_gp_mass_[numgp] = dr0dxi.Norm2();

      triad_interpolation_scheme_ptr->get_interpolated_quaternion_at_xi(
          qnew_gp_mass_[numgp], gausspoints_inertia.qxg[numgp][0]);

      // copy QnewGPmass_ to QconvGPmass_
      qconv_gp_mass_[numgp] = qnew_gp_mass_[numgp];

      wconv_gp_mass_[numgp].Clear();
      wnew_gp_mass_[numgp].Clear();
      aconv_gp_mass_[numgp].Clear();
      anew_gp_mass_[numgp].Clear();
      amodconv_gp_mass_[numgp].Clear();
      amodnew_gp_mass_[numgp].Clear();
      rttconv_gp_mass_[numgp].Clear();
      rttnew_gp_mass_[numgp].Clear();
      rttmodconv_gp_mass_[numgp].Clear();
      rttmodnew_gp_mass_[numgp].Clear();
      rtconv_gp_mass_[numgp].Clear();
      rtnew_gp_mass_[numgp].Clear();
      rconv_gp_mass_[numgp] = r0;
      rnew_gp_mass_[numgp] = r0;
    }


    /********************* Compute quantities required for damping/stochastic forces
     **********************
     *****************************************************************************************************/

    // compute Jacobi determinant at GPs for integration of damping/stochastic forces

    // Get the applied integration scheme
    CORE::FE::GaussRule1D gaussrule_damp_stoch =
        MyGaussRule(res_damp_stoch);  // TODO reuse/copy quantities if same integration scheme has
                                      // been applied above
    CORE::FE::IntegrationPoints1D gausspoints_damp_stoch(gaussrule_damp_stoch);

    // these quantities will later be used mainly for calculation of damping/stochastic terms ->
    // named 'dampstoch'
    qconv_gp_dampstoch_.resize(gausspoints_damp_stoch.nquad);
    qnew_gp_dampstoch_.resize(gausspoints_damp_stoch.nquad);
    jacobi_gp_dampstoch_.resize(gausspoints_damp_stoch.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_damp_stoch.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_damp_stoch, H_i_xi, distype, this->RefLength());

    // Loop through all GPs
    for (int numgp = 0; numgp < gausspoints_damp_stoch.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(pos_ref_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobi_gp_dampstoch_[numgp] = dr0dxi.Norm2();

      triad_interpolation_scheme_ptr->get_interpolated_quaternion_at_xi(
          qnew_gp_dampstoch_[numgp], gausspoints_damp_stoch.qxg[numgp][0]);

      // copy QnewGPdampstoch_ to QconvGPdampstoch_
      qconv_gp_dampstoch_[numgp] = qnew_gp_dampstoch_[numgp];
    }


    /********************* Compute quantities required for integration of Neumann lineloads
     ***************
     *****************************************************************************************************/

    // Get the applied integration scheme
    CORE::FE::GaussRule1D gaussrule_neumann = MyGaussRule(neumann_lineload);
    CORE::FE::IntegrationPoints1D gausspoints_neumann(gaussrule_neumann);

    // these quantities will later be used for calculation of Neumann lineloads
    jacobi_gp_neumannline_.resize(gausspoints_neumann.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_neumann.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_neumann, H_i_xi, distype, this->RefLength());

    // Loop through all GPs
    for (int numgp = 0; numgp < gausspoints_neumann.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(pos_ref_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobi_gp_neumannline_[numgp] = dr0dxi.Norm2();
    }
  }

  return;
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::GetPosAtXi(
    CORE::LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int numnodalvalues = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->num_node();

  std::vector<double> disp_centerline(3 * numnodalvalues * nnodecl, 0.0);

  /* we assume that either the full disp vector of this element or
   * disp_centerline (without rotational DoFs) is passed in this function call */
  if (disp.size() == 3 * numnodalvalues * nnodecl)
    disp_centerline = disp;
  else if (disp.size() == 3 * numnodalvalues * nnodecl + 3 * nnodetriad)
    extract_centerline_dof_values_from_element_state_vector(disp, disp_centerline);
  else
    FOUR_C_THROW(
        "size mismatch: expected either %d values for disp_centerline or "
        "%d values for full disp state vector of this element and got %d",
        3 * numnodalvalues * nnodecl, 3 * numnodalvalues * nnodecl + 3 * nnodetriad, disp.size());

  switch (nnodecl)
  {
    case 2:
    {
      if (this->hermite_centerline_interpolation())
      {
        CORE::LINALG::Matrix<12, 1> disp_totlag_centerline_fixedsize(disp_centerline.data());
        add_ref_values_disp_centerline<2, 2, double>(disp_totlag_centerline_fixedsize);
        Beam3Base::GetPosAtXi<2, 2, double>(pos, xi, disp_totlag_centerline_fixedsize);
      }
      else
      {
        CORE::LINALG::Matrix<6, 1> disp_totlag_centerline_fixedsize(disp_centerline.data());
        add_ref_values_disp_centerline<2, 1, double>(disp_totlag_centerline_fixedsize);
        Beam3Base::GetPosAtXi<2, 1, double>(pos, xi, disp_totlag_centerline_fixedsize);
      }
      break;
    }
    case 3:
    {
      CORE::LINALG::Matrix<9, 1> disp_totlag_centerline_fixedsize(disp_centerline.data());
      add_ref_values_disp_centerline<3, 1, double>(disp_totlag_centerline_fixedsize);
      Beam3Base::GetPosAtXi<3, 1, double>(pos, xi, disp_totlag_centerline_fixedsize);
      break;
    }
    case 4:
    {
      CORE::LINALG::Matrix<12, 1> disp_totlag_centerline_fixedsize(disp_centerline.data());
      add_ref_values_disp_centerline<4, 1, double>(disp_totlag_centerline_fixedsize);
      Beam3Base::GetPosAtXi<4, 1, double>(pos, xi, disp_totlag_centerline_fixedsize);
      break;
    }
    case 5:
    {
      CORE::LINALG::Matrix<15, 1> disp_totlag_centerline_fixedsize(disp_centerline.data());
      add_ref_values_disp_centerline<5, 1, double>(disp_totlag_centerline_fixedsize);
      Beam3Base::GetPosAtXi<5, 1, double>(pos, xi, disp_totlag_centerline_fixedsize);
      break;
    }
    default:
      FOUR_C_THROW("no valid number for number of centerline nodes");
  }

  return;
}

double DRT::ELEMENTS::Beam3r::GetJacobiFacAtXi(const double& xi) const
{
  double jacfac = 0.0;

  switch (this->NumCenterlineNodes())
  {
    case 2:
    {
      if (this->hermite_centerline_interpolation())
        jacfac = this->GetJacobiFacAtXi<2, 2>(xi);
      else
        jacfac = this->GetJacobiFacAtXi<2, 1>(xi);
      break;
    }
    case 3:
    {
      jacfac = this->GetJacobiFacAtXi<3, 1>(xi);
      break;
    }
    case 4:
    {
      jacfac = this->GetJacobiFacAtXi<4, 1>(xi);
      break;
    }
    case 5:
    {
      jacfac = this->GetJacobiFacAtXi<5, 1>(xi);
      break;
    }
    default:
      FOUR_C_THROW("no valid number for number of centerline nodes");
  }

  return jacfac;
}

void DRT::ELEMENTS::Beam3r::GetTriadAtXi(
    CORE::LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int numnodalvalues = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->num_node();

  std::vector<CORE::LINALG::Matrix<3, 1, double>> nodal_rotvecs(nnodetriad);

  /* we assume that either the full disp vector of this element or only
   * values for nodal rotation vectors are passed in this function call */
  if (disp.size() == 3 * nnodetriad)
  {
    for (unsigned int node = 0; node < nnodetriad; ++node)
      for (unsigned int i = 0; i < 3; ++i) nodal_rotvecs[node](i) = disp[node * 3 + i];
  }
  else if (disp.size() == 3 * numnodalvalues * nnodecl + 3 * nnodetriad)
  {
    extract_rot_vec_dof_values(disp, nodal_rotvecs);
  }
  else
  {
    FOUR_C_THROW(
        "size mismatch: expected either %d values for psi (rotation vecs) or "
        "%d values for for full disp state vector of this element and got %d",
        3 * nnodetriad, 3 * numnodalvalues * nnodecl + 3 * nnodetriad, disp.size());
  }

  // nodal triads
  std::vector<CORE::LINALG::Matrix<4, 1, double>> Qnode(nnodetriad);

  switch (nnodetriad)
  {
    case 2:
    {
      get_nodal_triads_from_disp_theta<2, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<2, double>(triad, xi, Qnode);
      break;
    }
    case 3:
    {
      get_nodal_triads_from_disp_theta<3, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<3, double>(triad, xi, Qnode);
      break;
    }
    case 4:
    {
      get_nodal_triads_from_disp_theta<4, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<4, double>(triad, xi, Qnode);
      break;
    }
    case 5:
    {
      get_nodal_triads_from_disp_theta<5, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<5, double>(triad, xi, Qnode);
      break;
    }
    default:
      FOUR_C_THROW("%d is no valid number of nodes for beam3r triad interpolation", nnodetriad);
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::get_generalized_interpolation_matrix_variations_at_xi(
    CORE::LINALG::SerialDenseMatrix& Ivar, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int vpernode = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->num_node();

  // safety check
  if (static_cast<unsigned int>(Ivar.numRows()) != 6 or
      static_cast<unsigned int>(Ivar.numCols()) != 3 * vpernode * nnodecl + 3 * nnodetriad)
    FOUR_C_THROW("size mismatch! expected %dx%d matrix and got %dx%d", 6,
        3 * vpernode * nnodecl + 3 * nnodetriad, Ivar.numRows(), Ivar.numCols());

  switch (nnodetriad)
  {
    case 2:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 12, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<2, 2, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        CORE::LINALG::Matrix<6, 18, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<2, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 3:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 18, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<3, 3, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        CORE::LINALG::Matrix<6, 21, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<3, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 4:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 24, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<4, 4, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        CORE::LINALG::Matrix<6, 24, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<4, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 5:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 30, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<5, 5, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        CORE::LINALG::Matrix<6, 27, double> Ivar_fixedsize(&Ivar(0, 0), true);
        get_generalized_interpolation_matrix_variations_at_xi<5, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    default:
      FOUR_C_THROW("Beam3r: no valid number of nodes specified");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::get_generalized_interpolation_matrix_variations_at_xi(
    CORE::LINALG::Matrix<6, 3 * vpernode * nnodecl + 3 * nnodetriad, double>& Ivar,
    const double& xi) const
{
  const unsigned int dofperclnode = 3 * vpernode;
  const unsigned int dofpertriadnode = 3;
  const unsigned int dofpercombinode = dofperclnode + dofpertriadnode;

  // these shape functions are used for the interpolation of the triad field
  // (so far always Lagrange polynomials of order 1...5)
  CORE::LINALG::Matrix<1, nnodetriad, double> I_i;
  // these shape functions are used for the interpolation of the beam centerline
  // (either cubic Hermite or Lagrange polynomials of order 1...5)
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> H_i;

  DRT::UTILS::BEAM::EvaluateShapeFunctionsAtXi<nnodetriad, 1>(xi, I_i, this->Shape());
  DRT::UTILS::BEAM::EvaluateShapeFunctionsAtXi<nnodecl, vpernode>(
      xi, H_i, this->Shape(), this->RefLength());

  Ivar.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
  {
    for (unsigned int inode = 0; inode < nnodecl; ++inode)
    {
      Ivar(idim, dofpercombinode * inode + idim) = H_i(vpernode * inode);
      Ivar(3 + idim, dofpercombinode * inode + 3 + idim) = I_i(inode);
      if (vpernode == 2) Ivar(idim, dofpercombinode * inode + 6 + idim) = H_i(vpernode * inode + 1);
    }
    // this loop is only entered in case of nnodetriad>nnodecl
    for (unsigned int inode = nnodecl; inode < nnodetriad; ++inode)
    {
      Ivar(3 + idim, dofperclnode * nnodecl + dofpertriadnode * inode + idim) = I_i(inode);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::get_generalized_interpolation_matrix_increments_at_xi(
    CORE::LINALG::SerialDenseMatrix& Iinc, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int vpernode = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->num_node();

  // safety check
  if (static_cast<unsigned int>(Iinc.numRows()) != 6 or
      static_cast<unsigned int>(Iinc.numCols()) != 3 * vpernode * nnodecl + 3 * nnodetriad)
    FOUR_C_THROW("size mismatch! expected %dx%d matrix and got %dx%d", 6,
        3 * vpernode * nnodecl + 3 * nnodetriad, Iinc.numRows(), Iinc.numCols());

  switch (nnodetriad)
  {
    case 2:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 12, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<2, 2, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        CORE::LINALG::Matrix<6, 18, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<2, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 3:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 18, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<3, 3, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        CORE::LINALG::Matrix<6, 21, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<3, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 4:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 24, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<4, 4, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        CORE::LINALG::Matrix<6, 24, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<4, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 5:
    {
      if (vpernode == 1)
      {
        CORE::LINALG::Matrix<6, 30, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<5, 5, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        CORE::LINALG::Matrix<6, 27, double> Iinc_fixedsize(&Iinc(0, 0), true);
        get_generalized_interpolation_matrix_increments_at_xi<5, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    default:
      FOUR_C_THROW("Beam3r: no valid number of nodes specified");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::get_generalized_interpolation_matrix_increments_at_xi(
    CORE::LINALG::Matrix<6, 3 * vpernode * nnodecl + 3 * nnodetriad, double>& Iinc,
    const double& xi, const std::vector<double>& disp) const
{
  const unsigned int dofperclnode = 3 * vpernode;
  const unsigned int dofpertriadnode = 3;
  const unsigned int dofpercombinode = dofperclnode + dofpertriadnode;

  // these shape functions are used for the interpolation of the beam centerline
  // (either cubic Hermite or Lagrange polynomials of order 1...5)
  CORE::LINALG::Matrix<1, vpernode * nnodecl, double> H_i;

  DRT::UTILS::BEAM::EvaluateShapeFunctionsAtXi<nnodecl, vpernode>(
      xi, H_i, this->Shape(), this->RefLength());

  // nodal triads in form of quaternions
  std::vector<CORE::LINALG::Matrix<4, 1, double>> Qnode(nnodetriad);

  get_nodal_triads_from_full_disp_vec_or_from_disp_theta<nnodetriad, double>(disp, Qnode);

  // vector with nnodetriad elements, who represent the 3x3-matrix-shaped interpolation
  // function \tilde{I}^nnode at a certain point xi according to (3.18), Jelenic 1999
  std::vector<CORE::LINALG::Matrix<3, 3, double>> Itilde(nnodetriad);
  compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads<nnodetriad, double>(
      Qnode, xi, Itilde);

  Iinc.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
  {
    for (unsigned int inode = 0; inode < nnodecl; ++inode)
    {
      Iinc(idim, dofpercombinode * inode + idim) = H_i(vpernode * inode);

      for (unsigned int jdim = 0; jdim < 3; ++jdim)
        Iinc(3 + idim, dofpercombinode * inode + 3 + jdim) = Itilde[inode](idim, jdim);

      if (vpernode == 2) Iinc(idim, dofpercombinode * inode + 6 + idim) = H_i(vpernode * inode + 1);
    }
    // this loop is only entered in case of nnodetriad>nnodecl
    for (unsigned int inode = nnodecl; inode < nnodetriad; ++inode)
    {
      for (unsigned int jdim = 0; jdim < 3; ++jdim)
        Iinc(3 + idim, dofperclnode * nnodecl + dofpertriadnode * inode + jdim) =
            Itilde[inode](idim, jdim);
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 | update (total) displacement vector and set nodal triads (as quaternions) grill 03/16|
 *------------------------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads(const std::vector<double>& disp,
    CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& disp_totlag_centerline,
    std::vector<CORE::LINALG::Matrix<4, 1, T>>& Q_i)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  // get current values of translational nodal DOFs in total Lagrangean manner (initial value +
  // disp) rotational DOFs need different handling, depending on whether FAD is used or not (see
  // comment below)
  extract_centerline_dof_values_from_element_state_vector<nnodecl, vpernode, T>(
      disp, disp_totlag_centerline);
  add_ref_values_disp_centerline<nnodecl, vpernode, T>(disp_totlag_centerline);

  // get current displacement values of rotational DOFs (i.e. relative rotation with respect to
  // reference config)
  std::vector<CORE::LINALG::Matrix<3, 1, double>> disptheta;
  disptheta.resize(nnodetriad);
  extract_rot_vec_dof_values<nnodetriad, nnodecl, vpernode, double>(disp, disptheta);

  // Compute current nodal triads
  get_nodal_triads_from_disp_theta<nnodetriad, T>(disptheta, Q_i);

  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // copy quaternions of nodal triads to class variable
    for (unsigned int i = 0; i < 4; ++i)
      qnewnode_[node](i) = CORE::FADUTILS::CastToDouble(Q_i[node](i));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables(
    CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1, FAD>& disp_totlag_centerline,
    std::vector<CORE::LINALG::Matrix<4, 1, FAD>>& Q_i) const
{
  const int dofperclnode = 3 * vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode + dofpertriadnode;

  // set differentiation variables for FAD: translational DOFs
  for (int dim = 0; dim < 3; ++dim)
  {
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      disp_totlag_centerline(dofperclnode * node + dim)
          .diff(
              dofpercombinode * node + dim, dofperclnode * nnodecl + dofpertriadnode * nnodetriad);

      // have Hermite interpolation? then set tangent DOFs as well
      if (vpernode == 2)
        disp_totlag_centerline(dofperclnode * node + 3 + dim)
            .diff(dofpercombinode * node + 6 + dim,
                dofperclnode * nnodecl + dofpertriadnode * nnodetriad);
    }
  }

  // rotation vector theta at a specific node in a total Lagrangean manner (with respect to global
  // reference coordinate system)
  std::vector<CORE::LINALG::Matrix<3, 1, FAD>> theta_totlag_i(nnodetriad);

  // compute nodal quaternions based on multiplicative increments of rotational DOFs
  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // compute physical total angle theta_totlag
    CORE::LARGEROTATIONS::quaterniontoangle(Q_i[node], theta_totlag_i[node]);
  }

  // set differentiation variables for FAD: rotational DOFs
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    for (unsigned int node = 0; node < nnodecl; ++node)
      theta_totlag_i[node](dim).diff(
          dofpercombinode * node + 3 + dim, dofperclnode * nnodecl + dofpertriadnode * nnodetriad);

    for (unsigned int node = nnodecl; node < nnodetriad; ++node)
      theta_totlag_i[node](dim).diff(dofperclnode * nnodecl + dofpertriadnode * node + dim,
          dofperclnode * nnodecl + dofpertriadnode * nnodetriad);
  }

  /* Attention: although the nodal quaternions Q_i have already been computed correctly, we need the
   * following step in order to track the dependency of subsequently calculated quantities via FAD
   */
  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    Q_i[node].PutScalar(0.0);
    CORE::LARGEROTATIONS::angletoquaternion(theta_totlag_i[node], Q_i[node]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector(
    const std::vector<double>& dofvec,
    CORE::LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& dofvec_centerline,
    bool add_reference_values) const
{
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperclnode = 3 * vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode + dofpertriadnode;

  if (dofvec.size() != dofperclnode * nnodecl + dofpertriadnode * this->num_node())
    FOUR_C_THROW("size mismatch: expected %d values for element state vector and got %d",
        dofperclnode * nnodecl + dofpertriadnode * this->num_node(), dofvec.size());

  // get current values for DOFs relevant for centerline interpolation
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_centerline(3 * vpernode * node + dim) = dofvec[dofpercombinode * node + dim];

      // have Hermite interpolation? then update tangent DOFs as well
      if (vpernode == 2)
        dofvec_centerline(3 * vpernode * node + 3 + dim) = dofvec[dofpercombinode * node + 6 + dim];
    }
  }

  if (add_reference_values) add_ref_values_disp_centerline<nnodecl, vpernode, T>(dofvec_centerline);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector(
    const std::vector<double>& dofvec, std::vector<double>& dofvec_centerline,
    bool add_reference_values) const
{
  const unsigned int vpernode = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();

  dofvec_centerline.resize(3 * vpernode * nnodecl, 0.0);

  switch (nnodecl)
  {
    case 2:
    {
      if (vpernode == 2)
      {
        // we use the method for CORE::LINALG fixed size matrix and create it as a view on the STL
        // vector
        CORE::LINALG::Matrix<12, 1> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);
        this->extract_centerline_dof_values_from_element_state_vector<2, 2, double>(
            dofvec, dofvec_centerline_fixedsize, add_reference_values);
      }
      else
      {
        CORE::LINALG::Matrix<6, 1> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);
        this->extract_centerline_dof_values_from_element_state_vector<2, 1, double>(
            dofvec, dofvec_centerline_fixedsize, add_reference_values);
      }
      break;
    }
    case 3:
    {
      CORE::LINALG::Matrix<9, 1> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);
      this->extract_centerline_dof_values_from_element_state_vector<3, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    case 4:
    {
      CORE::LINALG::Matrix<12, 1> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);
      this->extract_centerline_dof_values_from_element_state_vector<4, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    case 5:
    {
      CORE::LINALG::Matrix<15, 1> dofvec_centerline_fixedsize(dofvec_centerline.data(), true);
      this->extract_centerline_dof_values_from_element_state_vector<5, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    default:
      FOUR_C_THROW("no valid number for number of centerline nodes");
  }
}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values(const std::vector<double>& dofvec,
    std::vector<CORE::LINALG::Matrix<3, 1, T>>& dofvec_rotvec) const
{
  // nnodetriad: number of nodes used for triad interpolation
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperclnode = 3 * vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode + dofpertriadnode;

  if (dofvec.size() != dofperclnode * nnodecl + dofpertriadnode * nnodetriad)
    FOUR_C_THROW("size mismatch: expected %d values for element state vector and got %d",
        dofperclnode * nnodecl + dofpertriadnode * nnodetriad, dofvec.size());

  // get current values for DOFs relevant for triad interpolation
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_rotvec[node](dim) = dofvec[dofpercombinode * node + 3 + dim];
    }
    for (unsigned int node = nnodecl; node < nnodetriad; ++node)
    {
      dofvec_rotvec[node](dim) = dofvec[dofperclnode * nnodecl + dofpertriadnode * node + dim];
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values(const std::vector<double>& dofvec,
    std::vector<CORE::LINALG::Matrix<3, 1, double>>& dofvec_rotvec) const
{
  switch (this->num_node())
  {
    case 2:
    {
      if (this->hermite_centerline_interpolation())
      {
        this->extract_rot_vec_dof_values<2, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->extract_rot_vec_dof_values<2, 2, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 3:
    {
      if (this->hermite_centerline_interpolation())
      {
        this->extract_rot_vec_dof_values<3, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->extract_rot_vec_dof_values<3, 3, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 4:
    {
      if (this->hermite_centerline_interpolation())
      {
        this->extract_rot_vec_dof_values<4, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->extract_rot_vec_dof_values<4, 4, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 5:
    {
      if (this->hermite_centerline_interpolation())
      {
        this->extract_rot_vec_dof_values<5, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->extract_rot_vec_dof_values<5, 5, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    default:
      FOUR_C_THROW("no valid number for number of centerline nodes");
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_disp_theta(
    const std::vector<CORE::LINALG::Matrix<3, 1, double>>& disptheta,
    std::vector<CORE::LINALG::Matrix<4, 1, T>>& Qnode) const
{
  // initial nodal rotation vector in quaternion form
  CORE::LINALG::Matrix<4, 1> Q0;
  // rotational displacement at a certain node in quaternion form
  CORE::LINALG::Matrix<4, 1> deltaQ;

  // Compute nodal triads in quaternion form
  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // get initial nodal rotation vectors and transform to quaternions
    CORE::LARGEROTATIONS::angletoquaternion(theta0node_[node], Q0);

    // rotate initial triads by relative rotation vector from displacement vector (via quaternion
    // product)
    CORE::LARGEROTATIONS::angletoquaternion(disptheta[node], deltaQ);
    CORE::LARGEROTATIONS::quaternionproduct(Q0, deltaQ, Qnode[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and
    // intricate calculations
    Qnode[node].Scale(1.0 / CORE::FADUTILS::VectorNorm(Qnode[node]));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_full_disp_vec_or_from_disp_theta(
    const std::vector<T>& dispvec, std::vector<CORE::LINALG::Matrix<4, 1, T>>& Qnode) const
{
  const unsigned int vpernode = this->hermite_centerline_interpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();

  std::vector<CORE::LINALG::Matrix<3, 1, double>> nodal_rotvecs(nnodetriad);

  /* we assume that either the full disp vector of this element or only
   * values for nodal rotation vectors are passed in this function call */
  if (dispvec.size() == 3 * nnodetriad)
  {
    for (unsigned int node = 0; node < nnodetriad; ++node)
      for (unsigned int i = 0; i < 3; ++i) nodal_rotvecs[node](i) = dispvec[node * 3 + i];
  }
  else if (dispvec.size() == 3 * vpernode * nnodecl + 3 * nnodetriad)
  {
    extract_rot_vec_dof_values(dispvec, nodal_rotvecs);
  }
  else
  {
    FOUR_C_THROW(
        "size mismatch: expected either %d values for psi (rotation vecs) or "
        "%d values for for full disp state vector of this element and got %d",
        3 * nnodetriad, 3 * vpernode * nnodecl + 3 * nnodetriad, dispvec.size());
  }

  get_nodal_triads_from_disp_theta<nnodetriad, double>(nodal_rotvecs, Qnode);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::
    compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads(
        const std::vector<CORE::LINALG::Matrix<4, 1, T>>& Qnode, const double xi,
        std::vector<CORE::LINALG::Matrix<3, 3, T>>& Itilde) const
{
  // create object of triad interpolation scheme
  Teuchos::RCP<LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, T>>
      triad_interpolation_scheme_ptr =
          Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, T>());

  // reset triad interpolation scheme with nodal quaternions
  triad_interpolation_scheme_ptr->Reset(Qnode);

  triad_interpolation_scheme_ptr->get_nodal_generalized_rotation_interpolation_matrices_at_xi(
      Itilde, xi);
}

// explicit template instantations (some compilers do not export symboles defined above)
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<2, 2, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<2, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<3, 3, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<3, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<4, 4, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<4, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<5, 5, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::set_up_reference_geometry<5, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<2, 2, 1, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<6, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<3, 3, 1, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<9, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<4, 4, 1, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<5, 5, 1, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<15, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<2, 2, 2, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<3, 2, 2, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<4, 2, 2, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<5, 2, 2, double>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<2, 2, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<6, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<3, 3, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<9, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<4, 4, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<5, 5, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<15, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<2, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<3, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<4, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void
DRT::ELEMENTS::Beam3r::update_disp_tot_lag_and_nodal_triads<5, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, CORE::LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, Sacado::Fad::DFad<double>>>&);
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<2, 2, 1>(
    CORE::LINALG::Matrix<6, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<3, 3, 1>(
    CORE::LINALG::Matrix<9, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<4, 4, 1>(
    CORE::LINALG::Matrix<12, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<5, 5, 1>(
    CORE::LINALG::Matrix<15, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<2, 2, 2>(
    CORE::LINALG::Matrix<12, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<3, 2, 2>(
    CORE::LINALG::Matrix<12, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<4, 2, 2>(
    CORE::LINALG::Matrix<12, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::set_automatic_differentiation_variables<5, 2, 2>(
    CORE::LINALG::Matrix<12, 1, FAD>&, std::vector<CORE::LINALG::Matrix<4, 1, FAD>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector<2, 1,
    double>(const std::vector<double>&, CORE::LINALG::Matrix<6, 1, double>&, bool) const;
template void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector<3, 1,
    double>(const std::vector<double>&, CORE::LINALG::Matrix<9, 1, double>&, bool) const;
template void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector<4, 1,
    double>(const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&, bool) const;
template void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector<5, 1,
    double>(const std::vector<double>&, CORE::LINALG::Matrix<15, 1, double>&, bool) const;
template void DRT::ELEMENTS::Beam3r::extract_centerline_dof_values_from_element_state_vector<2, 2,
    double>(const std::vector<double>&, CORE::LINALG::Matrix<12, 1, double>&, bool) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<2, 2, 1, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<2, 2, 2, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<3, 3, 1, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<3, 2, 2, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<4, 4, 1, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<4, 2, 2, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<5, 5, 1, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::extract_rot_vec_dof_values<5, 2, 2, double>(
    const std::vector<double>&, std::vector<CORE::LINALG::Matrix<3, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_disp_theta<2, double>(
    const std::vector<CORE::LINALG::Matrix<3, 1, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_disp_theta<3, double>(
    const std::vector<CORE::LINALG::Matrix<3, 1, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_disp_theta<4, double>(
    const std::vector<CORE::LINALG::Matrix<3, 1, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_disp_theta<5, double>(
    const std::vector<CORE::LINALG::Matrix<3, 1, double>>&,
    std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_full_disp_vec_or_from_disp_theta<2,
    double>(const std::vector<double>&, std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_full_disp_vec_or_from_disp_theta<3,
    double>(const std::vector<double>&, std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_full_disp_vec_or_from_disp_theta<4,
    double>(const std::vector<double>&, std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void DRT::ELEMENTS::Beam3r::get_nodal_triads_from_full_disp_vec_or_from_disp_theta<5,
    double>(const std::vector<double>&, std::vector<CORE::LINALG::Matrix<4, 1, double>>&) const;
template void
DRT::ELEMENTS::Beam3r::compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads<2,
    double>(const std::vector<CORE::LINALG::Matrix<4, 1, double>>&, const double,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&) const;
template void
DRT::ELEMENTS::Beam3r::compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads<3,
    double>(const std::vector<CORE::LINALG::Matrix<4, 1, double>>&, const double,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&) const;
template void
DRT::ELEMENTS::Beam3r::compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads<4,
    double>(const std::vector<CORE::LINALG::Matrix<4, 1, double>>&, const double,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&) const;
template void
DRT::ELEMENTS::Beam3r::compute_generalized_nodal_rotation_interpolation_matrix_from_nodal_triads<5,
    double>(const std::vector<CORE::LINALG::Matrix<4, 1, double>>&, const double,
    std::vector<CORE::LINALG::Matrix<3, 3, double>>&) const;

FOUR_C_NAMESPACE_CLOSE
