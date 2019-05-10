/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam3r.cpp

\brief 3D nonlinear Reissner beam element

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3r.H"

#include "triad_interpolation_local_rotation_vectors.H"
#include "beam_spatial_discretization_utils.H"

// Todo @grill: check for obsolete header inclusions
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3rType DRT::ELEMENTS::Beam3rType::instance_;

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3rType& DRT::ELEMENTS::Beam3rType::Instance() { return instance_; }

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Beam3rType::Create(const std::vector<char>& data)
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
void DRT::ELEMENTS::Beam3rType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  DRT::ELEMENTS::Beam3r* currele = dynamic_cast<DRT::ELEMENTS::Beam3r*>(dwele);
  if (!currele) dserror("cast to Beam3r* failed");

  if (currele->HermiteCenterlineInterpolation() or currele->NumNode() > 2)
  {
    dserror(
        "method NodalBlockInformation not implemented for element type beam3r in case of Hermite "
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
void DRT::ELEMENTS::Beam3rType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  dserror("method ComputeNullSpace not implemented for element type beam3r!");
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3rType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["BEAM3R"];

  // note: LIN2 refers to linear Lagrange interpolation of centerline AND triad field
  defs["LIN2"]
      .AddIntVector("LIN2", 2)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 6)
      .AddOptionalTag("FAD");

  // note: LIN3 refers to quadratic Lagrange interpolation of centerline AND triad field
  defs["LIN3"]
      .AddIntVector("LIN3", 3)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 9)
      .AddOptionalTag("FAD");

  // note: LIN4 refers to cubic Lagrange interpolation of centerline AND triad field
  defs["LIN4"]
      .AddIntVector("LIN4", 4)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 12)
      .AddOptionalTag("FAD");

  // note: LIN5 refers to quartic Lagrange interpolation of centerline AND triad field
  defs["LIN5"]
      .AddIntVector("LIN5", 5)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 15)
      .AddOptionalTag("FAD");

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LIN2 refers to linear Lagrange interpolation of the triad field*/
  defs["HERM2LIN2"]
      .AddIntVector("HERM2LIN2", 2)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 6)
      .AddOptionalTag("FAD");

  /* note: HERM2 refers to cubic order Hermite interpolation of centerline (2 nodes)
   *       LIN3 refers to quadratic Lagrange interpolation of the triad field*/
  defs["HERM2LIN3"]
      .AddIntVector("HERM2LIN3", 3)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 9)
      .AddOptionalTag("FAD");

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LIN4 refers to cubic Lagrange interpolation of the triad field*/
  defs["HERM2LIN4"]
      .AddIntVector("HERM2LIN4", 4)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 12)
      .AddOptionalTag("FAD");

  /* note: HERM2 refers to cubic Hermite interpolation of centerline (2 nodes)
   *       LIN5 refers to quartic Lagrange interpolation of the triad field*/
  defs["HERM2LIN5"]
      .AddIntVector("HERM2LIN5", 5)
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 15)
      .AddOptionalTag("FAD");
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
    if (!currele) dserror("cast to Beam3r* failed");

    // reference node position
    std::vector<double> xrefe;
    std::vector<double> rotrefe;

    /* the triad field is discretized with Lagrange polynomials of order NumNode()-1;
     * the centerline is either discretized in the same way or with 3rd order Hermite polynomials;
     * in case of Hermite interpolation of the centerline, always the two boundary nodes are used
     * for centerline interpolation*/
    const bool centerline_hermite = currele->HermiteCenterlineInterpolation();

    // nnodetriad: number of nodes used for interpolation of triad field
    // nnodecl: number of nodes used for interpolation of centerline
    // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
    const int nnodetriad = currele->NumNode();
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
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> periodic_boundingbox =
        Teuchos::rcp(new GEO::MESHFREE::BoundingBox());
    periodic_boundingbox->Init();  // no Setup() call needed here

    std::vector<double> disp_shift;
    int numdof = currele->NumDofPerNode(*(currele->Nodes()[0]));
    disp_shift.resize(numdof * nnodecl);
    for (unsigned int i = 0; i < disp_shift.size(); ++i) disp_shift[i] = 0.0;
    if (periodic_boundingbox->HavePBC())
      currele->UnShiftNodePosition(disp_shift, *periodic_boundingbox);

    for (int node = 0; node < nnodecl; ++node)
    {
      if (currele->Nodes()[node] == NULL)
        dserror("beam3r: Cannot get nodes in order to compute reference configuration");

      for (unsigned int dim = 0; dim < 3; ++dim)
        xrefe[node * 3 + dim] = currele->Nodes()[node]->X()[dim] + disp_shift[node * numdof + dim];
    }

    // SetUpReferenceGeometry is a templated function
    switch (nnodetriad)
    {
      case 2:
      {
        if (!centerline_hermite)
          currele->SetUpReferenceGeometry<2, 2, 1>(xrefe, rotrefe);
        else
          currele->SetUpReferenceGeometry<2, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 3:
      {
        if (!centerline_hermite)
          currele->SetUpReferenceGeometry<3, 3, 1>(xrefe, rotrefe);
        else
          currele->SetUpReferenceGeometry<3, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 4:
      {
        if (!centerline_hermite)
          currele->SetUpReferenceGeometry<4, 4, 1>(xrefe, rotrefe);
        else
          currele->SetUpReferenceGeometry<4, 2, 2>(xrefe, rotrefe);
        break;
      }
      case 5:
      {
        if (!centerline_hermite)
          currele->SetUpReferenceGeometry<5, 5, 1>(xrefe, rotrefe);
        else
          currele->SetUpReferenceGeometry<5, 2, 2>(xrefe, rotrefe);
        break;
      }
      default:
        dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
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
      stiff_ptc_(true),
      useFAD_(false),
      centerline_hermite_(false),
      isinit_(false),
      jacobiGPelastf_(0),
      jacobiGPelastm_(0),
      jacobiGPmass_(0),
      jacobiGPdampstoch_(0),
      jacobiGPneumannline_(0),
      Eint_(0.0),
      Ekin_(0.0),
      Ekintorsion_(0.0),
      Ekinbending_(0.0),
      Ekintrans_(0.0),
      L_(true),
      P_(true)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3r::Beam3r(const DRT::ELEMENTS::Beam3r& old)
    : DRT::ELEMENTS::Beam3Base(old),
      useFAD_(old.useFAD_),
      centerline_hermite_(old.centerline_hermite_),
      isinit_(old.isinit_),
      reflength_(old.reflength_),
      theta0node_(old.theta0node_),
      Trefnode_(old.Trefnode_),
      Tcurrnode_(old.Tcurrnode_),
      KrefGP_(old.KrefGP_),
      GammarefGP_(old.GammarefGP_),
      jacobiGPelastf_(old.jacobiGPelastf_),
      jacobiGPelastm_(old.jacobiGPelastm_),
      jacobiGPmass_(old.jacobiGPmass_),
      jacobiGPdampstoch_(old.jacobiGPdampstoch_),
      jacobiGPneumannline_(old.jacobiGPneumannline_),
      Qconvnode_(old.Qconvnode_),
      Qnewnode_(old.Qnewnode_),
      QconvGPmass_(old.QconvGPmass_),
      QnewGPmass_(old.QnewGPmass_),
      wconvGPmass_(old.wconvGPmass_),
      wnewGPmass_(old.wnewGPmass_),
      aconvGPmass_(old.aconvGPmass_),
      anewGPmass_(old.anewGPmass_),
      amodconvGPmass_(old.amodconvGPmass_),
      amodnewGPmass_(old.amodnewGPmass_),
      rttconvGPmass_(old.rttconvGPmass_),
      rttnewGPmass_(old.rttnewGPmass_),
      rttmodconvGPmass_(old.rttmodconvGPmass_),
      rttmodnewGPmass_(old.rttmodnewGPmass_),
      rtconvGPmass_(old.rtconvGPmass_),
      rtnewGPmass_(old.rtnewGPmass_),
      rconvGPmass_(old.rconvGPmass_),
      rnewGPmass_(old.rnewGPmass_),
      QconvGPdampstoch_(old.QconvGPdampstoch_),
      QnewGPdampstoch_(old.QnewGPdampstoch_),
      Eint_(old.Eint_),
      Ekin_(old.Ekin_),
      Ekintorsion_(old.Ekintorsion_),
      Ekinbending_(old.Ekinbending_),
      Ekintrans_(old.Ekintrans_),
      L_(old.L_),
      P_(old.P_),
      Kmax_(old.Kmax_),
      axial_strain_GP_elastf_(old.axial_strain_GP_elastf_),
      shear_strain_2_GP_elastf_(old.shear_strain_2_GP_elastf_),
      shear_strain_3_GP_elastf_(old.shear_strain_3_GP_elastf_),
      twist_GP_elastm_(old.twist_GP_elastm_),
      curvature_2_GP_elastm_(old.curvature_2_GP_elastm_),
      curvature_3_GP_elastm_(old.curvature_3_GP_elastm_),
      material_axial_force_GP_elastf_(old.material_axial_force_GP_elastf_),
      material_shear_force_2_GP_elastf_(old.material_shear_force_2_GP_elastf_),
      material_shear_force_3_GP_elastf_(old.material_shear_force_3_GP_elastf_),
      material_torque_GP_elastm_(old.material_torque_GP_elastm_),
      material_bending_moment_2_GP_elastm_(old.material_bending_moment_2_GP_elastm_),
      material_bending_moment_3_GP_elastm_(old.material_bending_moment_3_GP_elastm_)
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
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3r::~Beam3r() { return; }

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
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3r::Shape() const
{
  int numnodes = NumNode();
  switch (numnodes)
  {
    case 2:
      return line2;
      break;
    case 3:
      return line3;
      break;
    case 4:
      return line4;
      break;
    case 5:
      return line5;
      break;
    default:
      dserror("Only Line2, Line3, Line4 and Line5 elements are implemented.");
      break;
  }

  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Beam3Base::Pack(data);

  // add all class variables of beam3r element
  AddtoPack(data, useFAD_);
  AddtoPack(data, centerline_hermite_);
  AddtoPack(data, isinit_);
  AddtoPack(data, reflength_);
  AddtoPack<3, 1>(data, theta0node_);
  AddtoPack<3, 1>(data, Trefnode_);
  AddtoPack<3, 1>(data, Tcurrnode_);
  AddtoPack<3, 1>(data, KrefGP_);
  AddtoPack<3, 1>(data, GammarefGP_);
  AddtoPack(data, jacobiGPelastf_);
  AddtoPack(data, jacobiGPelastm_);
  AddtoPack(data, jacobiGPmass_);
  AddtoPack(data, jacobiGPdampstoch_);
  AddtoPack(data, jacobiGPneumannline_);
  AddtoPack<4, 1>(data, Qconvnode_);
  AddtoPack<4, 1>(data, Qnewnode_);
  AddtoPack<4, 1>(data, QconvGPmass_);
  AddtoPack<4, 1>(data, QnewGPmass_);
  AddtoPack<3, 1>(data, wconvGPmass_);
  AddtoPack<3, 1>(data, wnewGPmass_);
  AddtoPack<3, 1>(data, aconvGPmass_);
  AddtoPack<3, 1>(data, anewGPmass_);
  AddtoPack<3, 1>(data, amodnewGPmass_);
  AddtoPack<3, 1>(data, amodconvGPmass_);
  AddtoPack<3, 1>(data, rttconvGPmass_);
  AddtoPack<3, 1>(data, rttnewGPmass_);
  AddtoPack<3, 1>(data, rttmodconvGPmass_);
  AddtoPack<3, 1>(data, rttmodnewGPmass_);
  AddtoPack<3, 1>(data, rtconvGPmass_);
  AddtoPack<3, 1>(data, rtnewGPmass_);
  AddtoPack<3, 1>(data, rconvGPmass_);
  AddtoPack<3, 1>(data, rnewGPmass_);
  AddtoPack<4, 1>(data, QconvGPdampstoch_);
  AddtoPack<4, 1>(data, QnewGPdampstoch_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Beam3Base::Unpack(basedata);

  // extract all class variables of beam3r element
  useFAD_ = ExtractInt(position, data);
  centerline_hermite_ = ExtractInt(position, data);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack(position, data, reflength_);
  ExtractfromPack<3, 1>(position, data, theta0node_);
  ExtractfromPack<3, 1>(position, data, Trefnode_);
  ExtractfromPack<3, 1>(position, data, Tcurrnode_);
  ExtractfromPack<3, 1>(position, data, KrefGP_);
  ExtractfromPack<3, 1>(position, data, GammarefGP_);
  ExtractfromPack(position, data, jacobiGPelastf_);
  ExtractfromPack(position, data, jacobiGPelastm_);
  ExtractfromPack(position, data, jacobiGPmass_);
  ExtractfromPack(position, data, jacobiGPdampstoch_);
  ExtractfromPack(position, data, jacobiGPneumannline_);
  ExtractfromPack<4, 1>(position, data, Qconvnode_);
  ExtractfromPack<4, 1>(position, data, Qnewnode_);
  ExtractfromPack<4, 1>(position, data, QconvGPmass_);
  ExtractfromPack<4, 1>(position, data, QnewGPmass_);
  ExtractfromPack<3, 1>(position, data, wconvGPmass_);
  ExtractfromPack<3, 1>(position, data, wnewGPmass_);
  ExtractfromPack<3, 1>(position, data, aconvGPmass_);
  ExtractfromPack<3, 1>(position, data, anewGPmass_);
  ExtractfromPack<3, 1>(position, data, amodconvGPmass_);
  ExtractfromPack<3, 1>(position, data, amodnewGPmass_);
  ExtractfromPack<3, 1>(position, data, rttconvGPmass_);
  ExtractfromPack<3, 1>(position, data, rttnewGPmass_);
  ExtractfromPack<3, 1>(position, data, rttmodconvGPmass_);
  ExtractfromPack<3, 1>(position, data, rttmodnewGPmass_);
  ExtractfromPack<3, 1>(position, data, rtconvGPmass_);
  ExtractfromPack<3, 1>(position, data, rtnewGPmass_);
  ExtractfromPack<3, 1>(position, data, rconvGPmass_);
  ExtractfromPack<3, 1>(position, data, rnewGPmass_);
  ExtractfromPack<4, 1>(position, data, QconvGPdampstoch_);
  ExtractfromPack<4, 1>(position, data, QnewGPdampstoch_);

  // NOT communicated
  Eint_ = 0.0;
  Ekin_ = 0.0;
  Ekintorsion_ = 0.0;
  Ekinbending_ = 0.0;
  Ekintrans_ = 0.0;
  L_.Clear();
  P_.Clear();
  Kmax_ = 0.0;
  axial_strain_GP_elastf_.clear();
  shear_strain_2_GP_elastf_.clear();
  shear_strain_3_GP_elastf_.clear();
  twist_GP_elastm_.clear();
  curvature_2_GP_elastm_.clear();
  curvature_3_GP_elastm_.clear();
  material_axial_force_GP_elastf_.clear();
  material_shear_force_2_GP_elastf_.clear();
  material_shear_force_3_GP_elastf_.clear();
  material_torque_GP_elastm_.clear();
  material_bending_moment_2_GP_elastm_.clear();
  material_bending_moment_3_GP_elastm_.clear();
  spatial_x_force_GP_elastf_.clear();
  spatial_y_force_2_GP_elastf_.clear();
  spatial_z_force_3_GP_elastf_.clear();
  spatial_x_moment_GP_elastm_.clear();
  spatial_y_moment_2_GP_elastm_.clear();
  spatial_z_moment_3_GP_elastm_.clear();

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Beam3r::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 | determine Gauss rule from purpose and interpolation scheme grill 03/16|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam3r::MyGaussRule(
    const IntegrationPurpose intpurpose) const
{
  const DRT::Element::DiscretizationType distype = this->Shape();

  switch (intpurpose)
  {
    // anti-locking: reduced integration of elastic residual contributions from forces (-> 'Gamma
    // terms')
    case res_elastic_force:
    {
      switch (distype)
      {
        case line2:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_1point;
          else
            return DRT::UTILS::intrule_line_lobatto3point;
        }
        case line3:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_2point;
          else
            return DRT::UTILS::intrule_line_lobatto3point;
        }
        case line4:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_3point;
          else
            return DRT::UTILS::intrule_line_lobatto3point;
        }
        case line5:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_4point;
          else
            return DRT::UTILS::intrule_line_lobatto3point;
        }
        default:
        {
          dserror("unknown discretization type!");
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
        case line2:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_1point;
          else
            return DRT::UTILS::intrule_line_2point;
        }
        case line3:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_2point;
          else
            return DRT::UTILS::intrule_line_3point;
        }
        case line4:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_3point;
          else
            return DRT::UTILS::intrule_line_4point;
        }
        case line5:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_4point;
          else
            return DRT::UTILS::intrule_line_5point;
        }
        default:
        {
          dserror("unknown discretization type!");
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
        case line2:
        {
          return DRT::UTILS::intrule_line_2point;
        }
        case line3:
        {
          return DRT::UTILS::intrule_line_3point;
        }
        case line4:
        {
          return DRT::UTILS::intrule_line_4point;
        }
        case line5:
        {
          return DRT::UTILS::intrule_line_5point;
        }
        default:
        {
          dserror("unknown discretization type!");
          break;
        }
      }
      break;
    }

    // 'full' integration of damping and stochastic contributions
    case res_damp_stoch:
    {
      return DRT::UTILS::intrule_line_4point;
    }

    /* 'full' integration of Neumann line loads
     * higher order Gauss quadrature scheme may prove useful in case of abnormal convergence
     * behaviour due to 'complex' line loads*/
    case neumann_lineload:
    {
      switch (distype)
      {
        case line2:
        {
          if (!centerline_hermite_)
            return DRT::UTILS::intrule_line_1point;
          else
            return DRT::UTILS::intrule_line_2point;
        }
        case line3:
        {
          return DRT::UTILS::intrule_line_2point;
        }
        case line4:
        {
          return DRT::UTILS::intrule_line_3point;
        }
        case line5:
        {
          return DRT::UTILS::intrule_line_4point;
        }
        default:
        {
          dserror("unknown discretization type!");
          break;
        }
      }
      break;
    }

    default:
    {
      dserror("beam3r: unknown purpose for numerical quadrature!");
      break;
    }
  }

  return DRT::UTILS::intrule1D_undefined;
  ;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry(
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
    dserror("Only 3rd order Hermite interpolation of beam centerline implemented!");

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
      dserror(
          "size mismatch in given position vector for stress-free reference geometry of beam3r:"
          " expected %d and got %d entries!",
          3 * nnodecl, xrefe.size());

    if (rotrefe.size() != 3 * nnodetriad)
      dserror(
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
    DRT::Element::DiscretizationType distype = Shape();

    /* Note: index i refers to the i-th shape function (i = 0 ... nnode*vpernode-1)
     * the vectors store individual shape functions, NOT an assembled matrix of shape functions) */
    /* vector whose numgp-th element is a 1xnnode-matrix with all Lagrange shape functions evaluated
     * at the numgp-th GP these shape functions are used for the interpolation of the triad field*/
    std::vector<LINALG::TMatrix<double, 1, nnodetriad>> I_i;
    // same for derivatives
    std::vector<LINALG::TMatrix<double, 1, nnodetriad>> I_i_xi;

    /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite)
     * shape functions evaluated at the numgp-th GP these shape functions are used for the
     * interpolation of the beam centerline*/
    std::vector<LINALG::TMatrix<double, 1, vpernode * nnodecl>> H_i;
    // same for the derivatives
    std::vector<LINALG::TMatrix<double, 1, vpernode * nnodecl>> H_i_xi;

    // beside the nodal reference positions from xrefe, this vector also holds the reference
    // tangents in case of Hermite interpolation of the beam centerline
    LINALG::Matrix<3 * vpernode * nnodecl, 1> disp_refe_centerline;

    // initial curve in physical space and derivative with respect to curve parameter xi \in [-1;1]
    // on element level
    LINALG::Matrix<3, 1> r0;
    LINALG::Matrix<3, 1> dr0dxi;

    // dummy 3x1 vector
    LINALG::Matrix<3, 1> dummy(true);


    /********************************** Compute nodal quantities
     *******************************************
     *****************************************************************************************************/

    /********************* store given nodal triads as quaternions in class variable
     * *********************/
    Qnewnode_.resize(nnodetriad);
    Qconvnode_.resize(nnodetriad);

    // nodal triads in stress-free configuration
    for (unsigned int node = 0; node < nnodetriad; node++)
    {
      LINALG::Matrix<3, 1> rotvec(&rotrefe[3 * node]);
      LARGEROTATIONS::angletoquaternion(rotvec, Qnewnode_[node]);
    }

    Qconvnode_ = Qnewnode_;

    std::vector<LINALG::TMatrix<double, 4, 1>> Qnewnode;

    for (unsigned int inode = 0; inode < nnodetriad; ++inode)
      Qnewnode.push_back(LINALG::TMatrix<double, 4, 1>(Qnewnode_[inode], true));

    // reset triad interpolation with nodal quaternions
    triad_interpolation_scheme_ptr->Reset(Qnewnode);

    LINALG::Matrix<3, 3> Gref;
    Trefnode_.resize(nnodecl);

    for (unsigned int node = 0; node < nnodecl; node++)
    {
      /* Calculate the (initial reference triads) = (initial material triads) at the nodes out of
       * the angles theta0node_. So far the initial value for the relative angle is set to zero,
       * i.e. material coordinate system and reference system in the reference configuration
       * coincidence (only at the nodes)*/
      Gref.Clear();
      LARGEROTATIONS::quaterniontotriad(Qnewnode_[node], Gref);
      // store initial nodal tangents in class variable
      for (int i = 0; i < 3; i++) (Trefnode_[node])(i) = (Gref)(i, 0);

      // fill disp_refe_centerline with reference nodal centerline positions and tangents
      for (int dim = 0; dim < 3; ++dim)
      {
        disp_refe_centerline(3 * vpernode * node + dim) = xrefe[3 * node + dim];
        if (centerline_hermite_)
          disp_refe_centerline(3 * vpernode * node + 3 + dim) = (Trefnode_[node])(dim);
      }
    }

    /***************************** Compute the initial length of the element
     * ******************************/

    // note: in case of Hermite centerline interpolation: iteratively via Newton's method
    Calculate_reflength<nnodecl, vpernode>(disp_refe_centerline, BEAM3RLENGTHCALCNEWTONTOL);


    /************************ Compute quantities required for elasticity
     ***********************************
     *****************************************************************************************************/

    // interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11),
    // Jelenic 1999
    LINALG::Matrix<3, 1> Psi_l;
    /* derivative of interpolated local relative rotation \Psi^l with respect to arc-length
     * parameter at a certain Gauss point according to (3.11), Jelenic 1999*/
    LINALG::Matrix<3, 1> Psi_l_s;
    // triad at GP
    LINALG::Matrix<3, 3> Lambda;

    //*********************** preparation for residual contributions from forces
    //***************************

    // Get the applied integration scheme
    DRT::UTILS::IntegrationPoints1D gausspoints_elast_force(MyGaussRule(res_elastic_force));

    jacobiGPelastf_.resize(gausspoints_elast_force.nquad);
    GammarefGP_.resize(gausspoints_elast_force.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_elast_force.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_elast_force, H_i_xi, distype, this->RefLength());


    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(axial_strain_GP_elastf_.begin(), axial_strain_GP_elastf_.end(), 0.0);
    shear_strain_2_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(shear_strain_2_GP_elastf_.begin(), shear_strain_2_GP_elastf_.end(), 0.0);
    shear_strain_3_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(shear_strain_3_GP_elastf_.begin(), shear_strain_3_GP_elastf_.end(), 0.0);

    material_axial_force_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(material_axial_force_GP_elastf_.begin(), material_axial_force_GP_elastf_.end(), 0.0);
    material_shear_force_2_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(
        material_shear_force_2_GP_elastf_.begin(), material_shear_force_2_GP_elastf_.end(), 0.0);
    material_shear_force_3_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(
        material_shear_force_3_GP_elastf_.begin(), material_shear_force_3_GP_elastf_.end(), 0.0);

    spatial_x_force_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_x_force_GP_elastf_.begin(), spatial_x_force_GP_elastf_.end(), 0.0);
    spatial_y_force_2_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_y_force_2_GP_elastf_.begin(), spatial_y_force_2_GP_elastf_.end(), 0.0);
    spatial_z_force_3_GP_elastf_.resize(gausspoints_elast_force.nquad);
    std::fill(spatial_z_force_3_GP_elastf_.begin(), spatial_z_force_3_GP_elastf_.end(), 0.0);

    dummy.Clear();

    // Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
    for (int numgp = 0; numgp < gausspoints_elast_force.nquad; ++numgp)
    {
      Calc_r_xi<nnodecl, vpernode, double>(disp_refe_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point for under-integration
      jacobiGPelastf_[numgp] = dr0dxi.Norm2();

      // we need dr0ds for computestrain, just reuse dr0dxi from above for simplicity
      dr0dxi.Scale(1.0 / jacobiGPelastf_[numgp]);

      triad_interpolation_scheme_ptr->GetInterpolatedTriadAtXi(
          Lambda, gausspoints_elast_force.qxg[numgp][0]);

      /* compute material strain Gamma according to Jelenic 1999, eq. (2.12) for reference
       * configuration, i.e. call this function with gammaref=zerovector*/
      computeGamma<double>(dr0dxi, Lambda, dummy, GammarefGP_[numgp]);
    }

    //*********************** preparation for residual contributions from moments
    //***************************

    // Get the applied integration scheme
    DRT::UTILS::IntegrationPoints1D gausspoints_elast_moment(MyGaussRule(res_elastic_moment));

    jacobiGPelastm_.resize(gausspoints_elast_moment.nquad);
    KrefGP_.resize(gausspoints_elast_moment.nquad);

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
    twist_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(twist_GP_elastm_.begin(), twist_GP_elastm_.end(), 0.0);
    curvature_2_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(curvature_2_GP_elastm_.begin(), curvature_2_GP_elastm_.end(), 0.0);
    curvature_3_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(curvature_3_GP_elastm_.begin(), curvature_3_GP_elastm_.end(), 0.0);

    material_torque_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_torque_GP_elastm_.begin(), material_torque_GP_elastm_.end(), 0.0);
    material_bending_moment_2_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_bending_moment_2_GP_elastm_.begin(),
        material_bending_moment_2_GP_elastm_.end(), 0.0);
    material_bending_moment_3_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(material_bending_moment_3_GP_elastm_.begin(),
        material_bending_moment_3_GP_elastm_.end(), 0.0);

    spatial_x_moment_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_x_moment_GP_elastm_.begin(), spatial_x_moment_GP_elastm_.end(), 0.0);
    spatial_y_moment_2_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_y_moment_2_GP_elastm_.begin(), spatial_y_moment_2_GP_elastm_.end(), 0.0);
    spatial_z_moment_3_GP_elastm_.resize(gausspoints_elast_moment.nquad);
    std::fill(spatial_z_moment_3_GP_elastm_.begin(), spatial_z_moment_3_GP_elastm_.end(), 0.0);


    dummy.Clear();

    // Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
    for (int numgp = 0; numgp < gausspoints_elast_moment.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(disp_refe_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobiGPelastm_[numgp] = dr0dxi.Norm2();

      // we need dr0ds for computestrain, just reuse dr0dxi from above for simplicity
      dr0dxi.Scale(1.0 / jacobiGPelastm_[numgp]);

      triad_interpolation_scheme_ptr->GetInterpolatedLocalRotationVector(Psi_l, I_i[numgp]);

      triad_interpolation_scheme_ptr->GetInterpolatedLocalRotationVectorDerivative(
          Psi_l_s, I_i_xi[numgp], jacobiGPelastm_[numgp]);

      /* compute material curvature K according to Jelenic 1999, eq. (2.12) for reference
       * configuration, i.e. call this function with kapparef=zerovector*/
      computeK<double>(Psi_l, Psi_l_s, dummy, KrefGP_[numgp]);
    }


    /******************************* Compute quantities required for inertia
     *******************************
     *****************************************************************************************************/

    // Get the applied integration scheme
    DRT::UTILS::GaussRule1D gaussrule_inertia = MyGaussRule(res_inertia);
    DRT::UTILS::IntegrationPoints1D gausspoints_inertia(gaussrule_inertia);

    // these quantities will later be used mainly for calculation of inertia terms -> named 'mass'
    jacobiGPmass_.resize(gausspoints_inertia.nquad);
    QconvGPmass_.resize(gausspoints_inertia.nquad);
    QnewGPmass_.resize(gausspoints_inertia.nquad);
    wconvGPmass_.resize(gausspoints_inertia.nquad);
    wnewGPmass_.resize(gausspoints_inertia.nquad);
    aconvGPmass_.resize(gausspoints_inertia.nquad);
    anewGPmass_.resize(gausspoints_inertia.nquad);
    rttconvGPmass_.resize(gausspoints_inertia.nquad);
    rttnewGPmass_.resize(gausspoints_inertia.nquad);
    rttmodconvGPmass_.resize(gausspoints_inertia.nquad);
    rttmodnewGPmass_.resize(gausspoints_inertia.nquad);
    rtconvGPmass_.resize(gausspoints_inertia.nquad);
    rtnewGPmass_.resize(gausspoints_inertia.nquad);
    rconvGPmass_.resize(gausspoints_inertia.nquad);
    rnewGPmass_.resize(gausspoints_inertia.nquad);
    amodconvGPmass_.resize(gausspoints_inertia.nquad);
    amodnewGPmass_.resize(gausspoints_inertia.nquad);

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
      Calc_r_xi<nnodecl, vpernode, double>(disp_refe_centerline, H_i_xi[numgp], dr0dxi);
      Calc_r<nnodecl, vpernode, double>(disp_refe_centerline, H_i[numgp], r0);

      // Store Jacobi determinant at this Gauss point
      jacobiGPmass_[numgp] = dr0dxi.Norm2();

      triad_interpolation_scheme_ptr->GetInterpolatedQuaternionAtXi(
          QnewGPmass_[numgp], gausspoints_inertia.qxg[numgp][0]);

      // copy QnewGPmass_ to QconvGPmass_
      QconvGPmass_[numgp] = QnewGPmass_[numgp];

      wconvGPmass_[numgp].Clear();
      wnewGPmass_[numgp].Clear();
      aconvGPmass_[numgp].Clear();
      anewGPmass_[numgp].Clear();
      amodconvGPmass_[numgp].Clear();
      amodnewGPmass_[numgp].Clear();
      rttconvGPmass_[numgp].Clear();
      rttnewGPmass_[numgp].Clear();
      rttmodconvGPmass_[numgp].Clear();
      rttmodnewGPmass_[numgp].Clear();
      rtconvGPmass_[numgp].Clear();
      rtnewGPmass_[numgp].Clear();
      rconvGPmass_[numgp] = r0;
      rnewGPmass_[numgp] = r0;
    }


    /********************* Compute quantities required for damping/stochastic forces
     **********************
     *****************************************************************************************************/

    // compute Jacobi determinant at GPs for integration of damping/stochastic forces

    // Get the applied integration scheme
    DRT::UTILS::GaussRule1D gaussrule_damp_stoch =
        MyGaussRule(res_damp_stoch);  // TODO reuse/copy quantities if same integration scheme has
                                      // been applied above
    DRT::UTILS::IntegrationPoints1D gausspoints_damp_stoch(gaussrule_damp_stoch);

    // these quantities will later be used mainly for calculation of damping/stochastic terms ->
    // named 'dampstoch'
    QconvGPdampstoch_.resize(gausspoints_damp_stoch.nquad);
    QnewGPdampstoch_.resize(gausspoints_damp_stoch.nquad);
    jacobiGPdampstoch_.resize(gausspoints_damp_stoch.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_damp_stoch.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_damp_stoch, H_i_xi, distype, this->RefLength());

    // Loop through all GPs
    for (int numgp = 0; numgp < gausspoints_damp_stoch.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(disp_refe_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobiGPdampstoch_[numgp] = dr0dxi.Norm2();

      triad_interpolation_scheme_ptr->GetInterpolatedQuaternionAtXi(
          QnewGPdampstoch_[numgp], gausspoints_damp_stoch.qxg[numgp][0]);

      // copy QnewGPdampstoch_ to QconvGPdampstoch_
      QconvGPdampstoch_[numgp] = QnewGPdampstoch_[numgp];
    }


    /********************* Compute quantities required for integration of Neumann lineloads
     ***************
     *****************************************************************************************************/

    // Get the applied integration scheme
    DRT::UTILS::GaussRule1D gaussrule_neumann = MyGaussRule(neumann_lineload);
    DRT::UTILS::IntegrationPoints1D gausspoints_neumann(gaussrule_neumann);

    // these quantities will later be used for calculation of Neumann lineloads
    jacobiGPneumannline_.resize(gausspoints_neumann.nquad);

    // reuse variables for individual shape functions and resize to new numgp
    H_i_xi.resize(gausspoints_neumann.nquad);

    // evaluate all shape functions and derivatives with respect to element parameter xi at all
    // specified Gauss points
    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnodecl, vpernode>(
        gausspoints_neumann, H_i_xi, distype, this->RefLength());

    // Loop through all GPs
    for (int numgp = 0; numgp < gausspoints_neumann.nquad; numgp++)
    {
      Calc_r_xi<nnodecl, vpernode, double>(disp_refe_centerline, H_i_xi[numgp], dr0dxi);

      // Store Jacobi determinant at this Gauss point
      jacobiGPneumannline_[numgp] = dr0dxi.Norm2();
    }
  }

  return;
}

/*----------------------------------------------------------------------------------*
 |  return current tangent at node                                   mukherjee 10/14|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> DRT::ELEMENTS::Beam3r::Tcurr(const int NodeID)
{
  // Attention: strictly speaking, this returns the first basis vector of the nodal triad
  // which is NOT identical to the tangent of the centerline in case of Reissner theory

  // TODO maybe handle calculation and update of class variable Tcurrnode_ elsewhere
  //      and make this a pure access function, i.e. set const (just like Tref() )

  // TODO
  //  if (NumNode()>2) dserror("Beam3r::Tcurr() not implemented yet for nnode>2");
  //  if (centerline_hermite_) dserror("Beam3r::Tcurr() not implemented yet for Hermite
  //  interpolation of centerline");
  //
  //  for (int node=0; node< 2; node++)    // TODO what happens here?
  //  {
  //    const int* nodeids=this->NodeIds();
  //    if (nodeids[this->nodeI_]==NodeID)
  //    {
  //      LINALG::Matrix<3,3>DummyLambda(true);
  //      LARGEROTATIONS::quaterniontotriad(Qnewnode_[this->nodeI_],DummyLambda);
  //      Tcurrnode_[0].Clear();
  //      for (int i=0; i<3; i++)
  //        Tcurrnode_[0](i)= DummyLambda(i,0);
  //    }
  //    else if (nodeids[this->nodeJ_]==NodeID)
  //    {
  //      LINALG::Matrix<3,3>DummyLambda(true);
  //      LARGEROTATIONS::quaterniontotriad(Qnewnode_[this->nodeJ_],DummyLambda);
  //      Tcurrnode_[0].Clear();
  //
  //      for (int i=0; i<3; i++)
  //        Tcurrnode_[0](i)= DummyLambda(i,0);
  //    }
  //    else
  //      for (int i=0; i<3; i++)
  //        Tcurrnode_[0](i)= 0;
  //  }

  // Fixme @grill
  dserror(
      "Beam3r::Tcurr(): the implementation of this method is highly questionable "
      "and needs rework!");

  return Tcurrnode_[0];
}

/*----------------------------------------------------------------------------------*
 |  return reference tangent at first node                            mukherjee 04/15|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> DRT::ELEMENTS::Beam3r::Treffirst() const
{
  // TODO @grill
  if (NumNode() > 2)
    dserror(
        "Beam3r::Treffirst() is not intended for nnode>2 since tangent vector varies along "
        "centerline!");
  if (centerline_hermite_)
    dserror(
        "Beam3r::Treffirst() is not intended for Hermite interpolation of centerline since tangent "
        "vector varies along centerline!");

  LINALG::Matrix<3, 1> Tref;
  double norm = Trefnode_[0].Norm2();

  if (norm <= 1e-14)
    dserror("beam3r: cannot normalize tangent vector because its norm is close to zero!");

  Tref.Update(1.0 / norm, Trefnode_[0]);

  return Tref;
}

/*--------------------------------------------------------------------------------------------*
 | Calculates the element length                                                   meier 01/16|
 *--------------------------------------------------------------------------------------------*/
template <unsigned int nnode, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::Calculate_reflength(
    const LINALG::TMatrix<double, 3 * vpernode * nnode, 1>& disp_totlag_centerline,
    const double tolerance)
{
  // nnode: number of nodes
  // vpernode: interpolated values per node (2: Hermite, i.e. value + derivative of value)

  /* in case of Hermite centerline interpolation,
   * the length is computed iteratively via Newton's method: f(l)=l-int(|N'd|)dxi=0*/

  // safety check
  if (vpernode == 2 and nnode != 2)
    dserror(
        "the function Calculate_length is implemented for 3rd order Hermite interpolation of the "
        "centerline only!");

  // initial value for iteration: difference vector of positions of boundary nodes (always ID 0 and
  // 1) also trivial solution in case of linear Lagrange interpolation: (nnode==2 && vpernode==1)
  if (vpernode == 2 || (nnode == 2 && vpernode == 1))
  {
    LINALG::Matrix<3, 1> tempvec(true);
    for (int dim = 0; dim < 3; dim++)
    {
      tempvec(dim) = disp_totlag_centerline(3 * vpernode * 1 + dim) - disp_totlag_centerline(dim);
    }
    reflength_ = tempvec.Norm2();
  }
  else
    reflength_ = 0.0;

  // non-trivial solution
  if (!(nnode == 2 && vpernode == 1))
  {
    // Get 'more than enough' integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints =
        DRT::UTILS::IntegrationPoints1D(DRT::UTILS::intrule_line_10point);

    // Newton Iteration - Tolerance and residual
    double res = 1.0;

    // Integral-value for Gauss Integration
    double int_length = 0.0;
    // Derivative value of the length integral for Newton Iteration (=weighted sum over deriv_int,
    // gauss quadrature of: int(d/dl(|N'd|))dxi)
    double deriv_length = 0.0;
    // value needed to store the derivative of the integral at the GP: d/dl(|N'd|)
    double deriv_int = 0.0;

    // Matrices to store the function values of the shape functions
    std::vector<LINALG::TMatrix<double, 1, nnode * vpernode>> H_i_xi(gausspoints.nquad);

    DRT::UTILS::BEAM::EvaluateShapeFunctionDerivsAllGPs<nnode, vpernode>(
        gausspoints, H_i_xi, this->Shape(), this->RefLength());

    // current value of the derivative at the GP
    LINALG::Matrix<3, 1> r_xi;

    int numiter = 0;

    // in case of Lagrange interpolation, one integration loop is sufficient
    do
    {
      numiter++;
      int_length = 0.0;
      deriv_length = 0.0;

      // Loop through all GPs and compute the length and the derivative of the length
      for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
      {
        deriv_int = 0;

        // integral of the length
        deriv_int = 0;

        Calc_r_xi<nnode, vpernode, double>(disp_totlag_centerline, H_i_xi[numgp], r_xi);

        int_length += gausspoints.qwgt[numgp] * r_xi.Norm2();

        // derivative only needed for Newton's method in case of Hermite interpolation
        if (vpernode == 2)
        {
          // derivative of the integral of the length at GP
          for (int dim = 0; dim < 3; dim++)
          {
            deriv_int += (disp_totlag_centerline(3 + dim) * H_i_xi[numgp](1) / reflength_ +
                             disp_totlag_centerline(3 * vpernode * 1 + 3 + dim) * H_i_xi[numgp](3) /
                                 reflength_) *
                         r_xi(dim);
          }
          deriv_length += gausspoints.qwgt[numgp] * deriv_int / r_xi.Norm2();
        }
      }

      res = reflength_ - int_length;
      // Update
      reflength_ =
          reflength_ - res / (1 - deriv_length);  // the derivative of f(l)=l-int(|N'd|)dxi=0 is
                                                  // f'(l)=1-int(d/dl(|N'd|))dxi
    } while (vpernode != 1 && std::fabs(res) > tolerance && numiter < 100);

    if (numiter > 100)
      dserror(
          "failed to compute length of element in reference configuration iteratively: "
          "Newton unconverged!");
  }

  return;
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::GetPosAtXi(
    LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int numnodalvalues = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->NumNode();

  std::vector<double> disp_centerline(3 * numnodalvalues * nnodecl, 0.0);

  /* we assume that either the full disp vector of this element or
   * disp_centerline (without rotational DoFs) is passed in this function call */
  if (disp.size() == 3 * numnodalvalues * nnodecl)
    disp_centerline = disp;
  else if (disp.size() == 3 * numnodalvalues * nnodecl + 3 * nnodetriad)
    ExtractCenterlineDofValuesFromElementStateVector(disp, disp_centerline);
  else
    dserror(
        "size mismatch: expected either %d values for disp_centerline or "
        "%d values for full disp state vector of this element and got %d",
        3 * numnodalvalues * nnodecl, 3 * numnodalvalues * nnodecl + 3 * nnodetriad, disp.size());

  switch (nnodecl)
  {
    case 2:
    {
      if (this->HermiteCenterlineInterpolation())
      {
        LINALG::Matrix<12, 1> disp_totlag_centerline_fixedsize(&disp_centerline[0]);
        AddRefValuesDispCenterline<2, 2, double>(disp_totlag_centerline_fixedsize);
        pos = this->GetPosAtXi<2, 2>(xi, disp_totlag_centerline_fixedsize);
      }
      else
      {
        LINALG::Matrix<6, 1> disp_totlag_centerline_fixedsize(&disp_centerline[0]);
        AddRefValuesDispCenterline<2, 1, double>(disp_totlag_centerline_fixedsize);
        pos = this->GetPosAtXi<2, 1>(xi, disp_totlag_centerline_fixedsize);
      }
      break;
    }
    case 3:
    {
      LINALG::Matrix<9, 1> disp_totlag_centerline_fixedsize(&disp_centerline[0]);
      AddRefValuesDispCenterline<3, 1, double>(disp_totlag_centerline_fixedsize);
      pos = this->GetPosAtXi<3, 1>(xi, disp_totlag_centerline_fixedsize);
      break;
    }
    case 4:
    {
      LINALG::Matrix<12, 1> disp_totlag_centerline_fixedsize(&disp_centerline[0]);
      AddRefValuesDispCenterline<4, 1, double>(disp_totlag_centerline_fixedsize);
      pos = this->GetPosAtXi<4, 1>(xi, disp_totlag_centerline_fixedsize);
      break;
    }
    case 5:
    {
      LINALG::Matrix<15, 1> disp_totlag_centerline_fixedsize(&disp_centerline[0]);
      AddRefValuesDispCenterline<5, 1, double>(disp_totlag_centerline_fixedsize);
      pos = this->GetPosAtXi<5, 1>(xi, disp_totlag_centerline_fixedsize);
      break;
    }
    default:
      dserror("no valid number for number of centerline nodes");
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
      if (this->HermiteCenterlineInterpolation())
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
      dserror("no valid number for number of centerline nodes");
  }

  return jacfac;
}

void DRT::ELEMENTS::Beam3r::GetTriadAtXi(
    LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int numnodalvalues = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->NumNode();

  std::vector<LINALG::TMatrix<double, 3, 1>> nodal_rotvecs(nnodetriad);

  /* we assume that either the full disp vector of this element or only
   * values for nodal rotation vectors are passed in this function call */
  if (disp.size() == 3 * nnodetriad)
  {
    for (unsigned int node = 0; node < nnodetriad; ++node)
      for (unsigned int i = 0; i < 3; ++i) nodal_rotvecs[node](i) = disp[node * 3 + i];
  }
  else if (disp.size() == 3 * numnodalvalues * nnodecl + 3 * nnodetriad)
  {
    ExtractRotVecDofValues(disp, nodal_rotvecs);
  }
  else
  {
    dserror(
        "size mismatch: expected either %d values for psi (rotation vecs) or "
        "%d values for for full disp state vector of this element and got %d",
        3 * nnodetriad, 3 * numnodalvalues * nnodecl + 3 * nnodetriad, disp.size());
  }

  // nodal triads
  std::vector<LINALG::TMatrix<double, 4, 1>> Qnode(nnodetriad);

  switch (nnodetriad)
  {
    case 2:
    {
      GetNodalTriadsFromDispTheta<2, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<2, double>(triad, xi, Qnode);
      break;
    }
    case 3:
    {
      GetNodalTriadsFromDispTheta<3, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<3, double>(triad, xi, Qnode);
      break;
    }
    case 4:
    {
      GetNodalTriadsFromDispTheta<4, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<4, double>(triad, xi, Qnode);
      break;
    }
    case 5:
    {
      GetNodalTriadsFromDispTheta<5, double>(nodal_rotvecs, Qnode);
      this->GetTriadAtXi<5, double>(triad, xi, Qnode);
      break;
    }
    default:
      dserror("%d is no valid number of nodes for beam3r triad interpolation", nnodetriad);
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::GetGeneralizedInterpolationMatrixVariationsAtXi(
    LINALG::SerialDenseMatrix& Ivar, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int vpernode = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->NumNode();

  // safety check
  if ((unsigned int)Ivar.M() != 6 or (unsigned int)
                                             Ivar.N() != 3 * vpernode * nnodecl + 3 * nnodetriad)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6,
        3 * vpernode * nnodecl + 3 * nnodetriad, Ivar.M(), Ivar.N());

  switch (nnodetriad)
  {
    case 2:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 12> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<2, 2, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        LINALG::TMatrix<double, 6, 18> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<2, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 3:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 18> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<3, 3, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        LINALG::TMatrix<double, 6, 21> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<3, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 4:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 24> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<4, 4, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        LINALG::TMatrix<double, 6, 24> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<4, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    case 5:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 30> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<5, 5, 1>(Ivar_fixedsize, xi);
      }
      else
      {
        LINALG::TMatrix<double, 6, 27> Ivar_fixedsize(&Ivar(0, 0), true);
        GetGeneralizedInterpolationMatrixVariationsAtXi<5, 2, 2>(Ivar_fixedsize, xi);
      }
      break;
    }
    default:
      dserror("Beam3r: no valid number of nodes specified");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::GetGeneralizedInterpolationMatrixVariationsAtXi(
    LINALG::TMatrix<double, 6, 3 * vpernode * nnodecl + 3 * nnodetriad>& Ivar,
    const double& xi) const
{
  const unsigned int dofperclnode = 3 * vpernode;
  const unsigned int dofpertriadnode = 3;
  const unsigned int dofpercombinode = dofperclnode + dofpertriadnode;

  // these shape functions are used for the interpolation of the triad field
  // (so far always Lagrange polynomials of order 1...5)
  LINALG::TMatrix<double, 1, nnodetriad> I_i;
  // these shape functions are used for the interpolation of the beam centerline
  // (either cubic Hermite or Lagrange polynomials of order 1...5)
  LINALG::TMatrix<double, 1, vpernode * nnodecl> H_i;

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
void DRT::ELEMENTS::Beam3r::GetGeneralizedInterpolationMatrixIncrementsAtXi(
    LINALG::SerialDenseMatrix& Iinc, const double& xi, const std::vector<double>& disp) const
{
  const unsigned int vpernode = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();
  const unsigned int nnodetriad = this->NumNode();

  // safety check
  if ((unsigned int)Iinc.M() != 6 or (unsigned int)
                                             Iinc.N() != 3 * vpernode * nnodecl + 3 * nnodetriad)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6,
        3 * vpernode * nnodecl + 3 * nnodetriad, Iinc.M(), Iinc.N());

  switch (nnodetriad)
  {
    case 2:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 12> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<2, 2, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        LINALG::TMatrix<double, 6, 18> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<2, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 3:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 18> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<3, 3, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        LINALG::TMatrix<double, 6, 21> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<3, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 4:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 24> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<4, 4, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        LINALG::TMatrix<double, 6, 24> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<4, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    case 5:
    {
      if (vpernode == 1)
      {
        LINALG::TMatrix<double, 6, 30> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<5, 5, 1>(Iinc_fixedsize, xi, disp);
      }
      else
      {
        LINALG::TMatrix<double, 6, 27> Iinc_fixedsize(&Iinc(0, 0), true);
        GetGeneralizedInterpolationMatrixIncrementsAtXi<5, 2, 2>(Iinc_fixedsize, xi, disp);
      }
      break;
    }
    default:
      dserror("Beam3r: no valid number of nodes specified");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::GetGeneralizedInterpolationMatrixIncrementsAtXi(
    LINALG::TMatrix<double, 6, 3 * vpernode * nnodecl + 3 * nnodetriad>& Iinc, const double& xi,
    const std::vector<double>& disp) const
{
  const unsigned int dofperclnode = 3 * vpernode;
  const unsigned int dofpertriadnode = 3;
  const unsigned int dofpercombinode = dofperclnode + dofpertriadnode;

  // these shape functions are used for the interpolation of the beam centerline
  // (either cubic Hermite or Lagrange polynomials of order 1...5)
  LINALG::TMatrix<double, 1, vpernode * nnodecl> H_i;

  DRT::UTILS::BEAM::EvaluateShapeFunctionsAtXi<nnodecl, vpernode>(
      xi, H_i, this->Shape(), this->RefLength());

  // nodal triads in form of quaternions
  std::vector<LINALG::TMatrix<double, 4, 1>> Qnode(nnodetriad);

  GetNodalTriadsFromFullDispVecOrFromDispTheta<nnodetriad, double>(disp, Qnode);

  // vector with nnodetriad elements, who represent the 3x3-matrix-shaped interpolation
  // function \tilde{I}^nnode at a certain point xi according to (3.18), Jelenic 1999
  std::vector<LINALG::TMatrix<double, 3, 3>> Itilde(nnodetriad);
  ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads<nnodetriad, double>(
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
void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads(const std::vector<double>& disp,
    LINALG::TMatrix<T, 3 * vpernode * nnodecl, 1>& disp_totlag_centerline,
    std::vector<LINALG::TMatrix<T, 4, 1>>& Q_i)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  // get current values of translational nodal DOFs in total Lagrangean manner (initial value +
  // disp) rotational DOFs need different handling, depending on whether FAD is used or not (see
  // comment below)
  ExtractCenterlineDofValuesFromElementStateVector<nnodecl, vpernode, T>(
      disp, disp_totlag_centerline);
  AddRefValuesDispCenterline<nnodecl, vpernode, T>(disp_totlag_centerline);

  // get current displacement values of rotational DOFs (i.e. relative rotation with respect to
  // reference config)
  std::vector<LINALG::TMatrix<double, 3, 1>> disptheta;
  disptheta.resize(nnodetriad);
  ExtractRotVecDofValues<nnodetriad, nnodecl, vpernode, double>(disp, disptheta);

  // Compute current nodal triads
  GetNodalTriadsFromDispTheta<nnodetriad, T>(disptheta, Q_i);

  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // copy quaternions of nodal triads to class variable
    for (unsigned int i = 0; i < 4; ++i) Qnewnode_[node](i) = FADUTILS::CastToDouble(Q_i[node](i));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables(
    LINALG::TMatrix<FAD, 3 * vpernode * nnodecl, 1>& disp_totlag_centerline,
    std::vector<LINALG::TMatrix<FAD, 4, 1>>& Q_i) const
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
  std::vector<LINALG::TMatrix<FAD, 3, 1>> theta_totlag_i(nnodetriad);

  // compute nodal quaternions based on multiplicative increments of rotational DOFs
  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // compute physical total angle theta_totlag
    LARGEROTATIONS::quaterniontoangle(Q_i[node], theta_totlag_i[node]);
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
    LARGEROTATIONS::angletoquaternion(theta_totlag_i[node], Q_i[node]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector(
    const std::vector<double>& dofvec,
    LINALG::TMatrix<T, 3 * vpernode * nnodecl, 1>& dofvec_centerline,
    bool add_reference_values) const
{
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperclnode = 3 * vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode + dofpertriadnode;

  if (dofvec.size() != dofperclnode * nnodecl + dofpertriadnode * this->NumNode())
    dserror("size mismatch: expected %d values for element state vector and got %d",
        dofperclnode * nnodecl + dofpertriadnode * this->NumNode(), dofvec.size());

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

  if (add_reference_values) AddRefValuesDispCenterline<nnodecl, vpernode, T>(dofvec_centerline);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector(
    const std::vector<double>& dofvec, std::vector<double>& dofvec_centerline,
    bool add_reference_values) const
{
  const unsigned int vpernode = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();

  dofvec_centerline.resize(3 * vpernode * nnodecl, 0.0);

  switch (nnodecl)
  {
    case 2:
    {
      if (vpernode == 2)
      {
        // we use the method for LINALG fixed size matrix and create it as a view on the STL vector
        LINALG::Matrix<12, 1> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);
        this->ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
            dofvec, dofvec_centerline_fixedsize, add_reference_values);
      }
      else
      {
        LINALG::Matrix<6, 1> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);
        this->ExtractCenterlineDofValuesFromElementStateVector<2, 1, double>(
            dofvec, dofvec_centerline_fixedsize, add_reference_values);
      }
      break;
    }
    case 3:
    {
      LINALG::Matrix<9, 1> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);
      this->ExtractCenterlineDofValuesFromElementStateVector<3, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    case 4:
    {
      LINALG::Matrix<12, 1> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);
      this->ExtractCenterlineDofValuesFromElementStateVector<4, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    case 5:
    {
      LINALG::Matrix<15, 1> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);
      this->ExtractCenterlineDofValuesFromElementStateVector<5, 1, double>(
          dofvec, dofvec_centerline_fixedsize, add_reference_values);
      break;
    }
    default:
      dserror("no valid number for number of centerline nodes");
  }
}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues(
    const std::vector<double>& dofvec, std::vector<LINALG::TMatrix<T, 3, 1>>& dofvec_rotvec) const
{
  // nnodetriad: number of nodes used for triad interpolation
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperclnode = 3 * vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode + dofpertriadnode;

  if (dofvec.size() != dofperclnode * nnodecl + dofpertriadnode * nnodetriad)
    dserror("size mismatch: expected %d values for element state vector and got %d",
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
void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues(const std::vector<double>& dofvec,
    std::vector<LINALG::TMatrix<double, 3, 1>>& dofvec_rotvec) const
{
  switch (this->NumNode())
  {
    case 2:
    {
      if (this->HermiteCenterlineInterpolation())
      {
        this->ExtractRotVecDofValues<2, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->ExtractRotVecDofValues<2, 2, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 3:
    {
      if (this->HermiteCenterlineInterpolation())
      {
        this->ExtractRotVecDofValues<3, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->ExtractRotVecDofValues<3, 3, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 4:
    {
      if (this->HermiteCenterlineInterpolation())
      {
        this->ExtractRotVecDofValues<4, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->ExtractRotVecDofValues<4, 4, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    case 5:
    {
      if (this->HermiteCenterlineInterpolation())
      {
        this->ExtractRotVecDofValues<5, 2, 2, double>(dofvec, dofvec_rotvec);
      }
      else
      {
        this->ExtractRotVecDofValues<5, 5, 1, double>(dofvec, dofvec_rotvec);
      }
      break;
    }
    default:
      dserror("no valid number for number of centerline nodes");
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline(
    LINALG::TMatrix<T, 3 * vpernode * nnodecl, 1>& dofvec_centerline) const
{
  for (unsigned int dim = 0; dim < 3; ++dim)
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_centerline(3 * vpernode * node + dim) += Nodes()[node]->X()[dim];

      // have Hermite interpolation? then update tangent DOFs as well
      if (vpernode == 2) dofvec_centerline(3 * vpernode * node + 3 + dim) += Trefnode_[node](dim);
    }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromDispTheta(
    const std::vector<LINALG::TMatrix<double, 3, 1>>& disptheta,
    std::vector<LINALG::TMatrix<T, 4, 1>>& Qnode) const
{
  // initial nodal rotation vector in quaternion form
  LINALG::Matrix<4, 1> Q0;
  // rotational displacement at a certain node in quaternion form
  LINALG::Matrix<4, 1> deltaQ;

  // Compute nodal triads in quaternion form
  for (unsigned int node = 0; node < nnodetriad; ++node)
  {
    // get initial nodal rotation vectors and transform to quaternions
    LARGEROTATIONS::angletoquaternion(theta0node_[node], Q0);

    // rotate initial triads by relative rotation vector from displacement vector (via quaternion
    // product)
    LARGEROTATIONS::angletoquaternion(disptheta[node], deltaQ);
    LARGEROTATIONS::quaternionproduct(Q0, deltaQ, Qnode[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and
    // intricate calculations
    Qnode[node].Scale(1.0 / FADUTILS::VectorNorm(Qnode[node]));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromFullDispVecOrFromDispTheta(
    const std::vector<T>& dispvec, std::vector<LINALG::TMatrix<T, 4, 1>>& Qnode) const
{
  const unsigned int vpernode = this->HermiteCenterlineInterpolation() ? 2 : 1;
  const unsigned int nnodecl = this->NumCenterlineNodes();

  std::vector<LINALG::TMatrix<double, 3, 1>> nodal_rotvecs(nnodetriad);

  /* we assume that either the full disp vector of this element or only
   * values for nodal rotation vectors are passed in this function call */
  if (dispvec.size() == 3 * nnodetriad)
  {
    for (unsigned int node = 0; node < nnodetriad; ++node)
      for (unsigned int i = 0; i < 3; ++i) nodal_rotvecs[node](i) = dispvec[node * 3 + i];
  }
  else if (dispvec.size() == 3 * vpernode * nnodecl + 3 * nnodetriad)
  {
    ExtractRotVecDofValues(dispvec, nodal_rotvecs);
  }
  else
  {
    dserror(
        "size mismatch: expected either %d values for psi (rotation vecs) or "
        "%d values for for full disp state vector of this element and got %d",
        3 * nnodetriad, 3 * vpernode * nnodecl + 3 * nnodetriad, dispvec.size());
  }

  GetNodalTriadsFromDispTheta<nnodetriad, double>(nodal_rotvecs, Qnode);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodetriad, typename T>
void DRT::ELEMENTS::Beam3r::ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads(
    const std::vector<LINALG::TMatrix<T, 4, 1>>& Qnode, const double xi,
    std::vector<LINALG::TMatrix<T, 3, 3>>& Itilde) const
{
  // create object of triad interpolation scheme
  Teuchos::RCP<LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, T>>
      triad_interpolation_scheme_ptr =
          Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<nnodetriad, T>());

  // reset triad interpolation scheme with nodal quaternions
  triad_interpolation_scheme_ptr->Reset(Qnode);

  triad_interpolation_scheme_ptr->GetNodalGeneralizedRotationInterpolationMatricesAtXi(Itilde, xi);
}

// explicit template instantations (some compilers do not export symboles defined above)
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<2, 2, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<2, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<3, 3, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<3, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<4, 4, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<4, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<5, 5, 1>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry<5, 2, 2>(
    const std::vector<double>&, const std::vector<double>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<2, 2, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 6, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<3, 3, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 9, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<4, 4, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<5, 5, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 15, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<2, 2, 2, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<3, 2, 2, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<4, 2, 2, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<5, 2, 2, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<2, 2, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 6, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<3, 3, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 9, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<4, 4, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 12, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<5, 5, 1, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 15, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<2, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 12, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<3, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 12, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<4, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 12, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void
DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads<5, 2, 2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&, LINALG::TMatrix<Sacado::Fad::DFad<double>, 12, 1>&,
    std::vector<LINALG::TMatrix<Sacado::Fad::DFad<double>, 4, 1>>&);
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<2, 2, 1>(
    LINALG::TMatrix<FAD, 6, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<3, 3, 1>(
    LINALG::TMatrix<FAD, 9, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<4, 4, 1>(
    LINALG::TMatrix<FAD, 12, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<5, 5, 1>(
    LINALG::TMatrix<FAD, 15, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<2, 2, 2>(
    LINALG::TMatrix<FAD, 12, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<3, 2, 2>(
    LINALG::TMatrix<FAD, 12, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<4, 2, 2>(
    LINALG::TMatrix<FAD, 12, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::SetAutomaticDifferentiationVariables<5, 2, 2>(
    LINALG::TMatrix<FAD, 12, 1>&, std::vector<LINALG::TMatrix<FAD, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector<2, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 6, 1>&, bool) const;
template void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector<3, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 9, 1>&, bool) const;
template void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector<4, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&, bool) const;
template void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector<5, 1, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 15, 1>&, bool) const;
template void DRT::ELEMENTS::Beam3r::ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
    const std::vector<double>&, LINALG::TMatrix<double, 12, 1>&, bool) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<2, 2, 1, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<2, 2, 2, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<3, 3, 1, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<3, 2, 2, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<4, 4, 1, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<4, 2, 2, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<5, 5, 1, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::ExtractRotVecDofValues<5, 2, 2, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 3, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline<2, 1, double>(
    LINALG::TMatrix<double, 6, 1>&) const;
template void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline<3, 1, double>(
    LINALG::TMatrix<double, 9, 1>&) const;
template void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline<4, 1, double>(
    LINALG::TMatrix<double, 12, 1>&) const;
template void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline<5, 1, double>(
    LINALG::TMatrix<double, 15, 1>&) const;
template void DRT::ELEMENTS::Beam3r::AddRefValuesDispCenterline<2, 2, double>(
    LINALG::TMatrix<double, 12, 1>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromDispTheta<2, double>(
    const std::vector<LINALG::TMatrix<double, 3, 1>>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromDispTheta<3, double>(
    const std::vector<LINALG::TMatrix<double, 3, 1>>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromDispTheta<4, double>(
    const std::vector<LINALG::TMatrix<double, 3, 1>>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromDispTheta<5, double>(
    const std::vector<LINALG::TMatrix<double, 3, 1>>&,
    std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromFullDispVecOrFromDispTheta<2, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromFullDispVecOrFromDispTheta<3, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromFullDispVecOrFromDispTheta<4, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void DRT::ELEMENTS::Beam3r::GetNodalTriadsFromFullDispVecOrFromDispTheta<5, double>(
    const std::vector<double>&, std::vector<LINALG::TMatrix<double, 4, 1>>&) const;
template void
DRT::ELEMENTS::Beam3r::ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads<2, double>(
    const std::vector<LINALG::TMatrix<double, 4, 1>>&, const double,
    std::vector<LINALG::TMatrix<double, 3, 3>>&) const;
template void
DRT::ELEMENTS::Beam3r::ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads<3, double>(
    const std::vector<LINALG::TMatrix<double, 4, 1>>&, const double,
    std::vector<LINALG::TMatrix<double, 3, 3>>&) const;
template void
DRT::ELEMENTS::Beam3r::ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads<4, double>(
    const std::vector<LINALG::TMatrix<double, 4, 1>>&, const double,
    std::vector<LINALG::TMatrix<double, 3, 3>>&) const;
template void
DRT::ELEMENTS::Beam3r::ComputeGeneralizedNodalRotationInterpolationMatrixFromNodalTriads<5, double>(
    const std::vector<LINALG::TMatrix<double, 4, 1>>&, const double,
    std::vector<LINALG::TMatrix<double, 3, 3>>&) const;
