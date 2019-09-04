/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3k.H"

#include "triad_interpolation_local_rotation_vectors.H"

#include "../drt_beaminteraction/periodic_boundingbox.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

// Todo @grill: check for obsolete header inclusions
#include "../headers/FAD_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"

#include <Teuchos_TimeMonitor.hpp>

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3kType DRT::ELEMENTS::Beam3kType::instance_;

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3kType& DRT::ELEMENTS::Beam3kType::Instance() { return instance_; }

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Beam3kType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Beam3k* object = new DRT::ELEMENTS::Beam3k(-1, -1);
  object->Unpack(data);
  return object;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "BEAM3K")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3k(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create(const int id, const int owner)

{
  return Teuchos::rcp(new Beam3k(id, owner));
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3kType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  dserror("method 'NodalBlockInformation' not implemented for element type beam3k!");
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3kType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  dserror("Function not implemented yet.");
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3kType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["BEAM3K"];

  defs["LIN2"]
      .AddIntVector("LIN2", 2)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 6)
      .AddOptionalTag("FAD");

  defs["LIN3"]
      .AddIntVector("LIN3", 3)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 9)
      .AddOptionalTag("FAD");

  defs["LIN4"]
      .AddIntVector("LIN4", 4)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDoubleVector("TRIADS", 12)
      .AddOptionalTag("FAD");
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 01/16|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3kType::Initialize(DRT::Discretization& dis)
{
  // setting up geometric variables for Beam3k elements
  for (int num = 0; num < dis.NumMyColElements(); ++num)
  {
    // in case that current element is not a Beam3k element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    // if we get so far current element is a Beam3k element and  we get a pointer at it
    DRT::ELEMENTS::Beam3k* currele = dynamic_cast<DRT::ELEMENTS::Beam3k*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3k* failed");

    // reference node position
    std::vector<LINALG::Matrix<3, 1>> xrefe;

    const int nnode = currele->NumNode();

    // resize xrefe for the number of nodes to store
    xrefe.resize(nnode);

    // the next section is needed in case of periodic boundary conditions and a shifted
    // configuration (i.e. elements cut by the periodic boundary) in the input file
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> periodic_boundingbox =
        Teuchos::rcp(new GEO::MESHFREE::BoundingBox());
    periodic_boundingbox->Init();  // no Setup() call needed here

    std::vector<double> disp_shift;
    int numdof = currele->NumDofPerNode(*(currele->Nodes()[0]));
    disp_shift.resize(numdof * nnode);
    for (unsigned int i = 0; i < disp_shift.size(); ++i) disp_shift[i] = 0.0;
    if (periodic_boundingbox->HavePBC())
      currele->UnShiftNodePosition(disp_shift, *periodic_boundingbox);

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int node = 0; node < nnode; ++node)  // element has k nodes
        for (int dof = 0; dof < 3; ++dof)       // element node has three coordinates x1, x2 and x3
        {
          xrefe[node](dof) = currele->Nodes()[node]->X()[dof] + disp_shift[node * numdof + dof];
        }
    }

    // Set up all geometrical (triads, curvatures, jacobians etc.) quantities describing the
    // (initial) reference geometry
    currele->SetUpReferenceGeometry(xrefe);
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(int id, int owner)
    : DRT::ELEMENTS::Beam3Base(id, owner),
      useFAD_(false),
      isinit_(false),
      T0_(0),
      T_(0),
      theta0_(0),
      Qrefconv_(0),
      Qrefnew_(0),
      K0_(0),
      length_(0.0),
      jacobi_(0.0),
      jacobi2_(0.0),
      jacobi_cp_(0.0),
      jacobi2_cp_(0.0),
      rotvec_(false),
      weakkirchhoff_(false),
      Eint_(0.0),
      Ekin_(0.0),
      Qconvmass_(0),
      Qnewmass_(0),
      wconvmass_(0),
      wnewmass_(0),
      aconvmass_(0),
      anewmass_(0),
      amodconvmass_(0),
      amodnewmass_(0),
      rttconvmass_(0),
      rttnewmass_(0),
      rttmodconvmass_(0),
      rttmodnewmass_(0),
      rtconvmass_(0),
      rtnewmass_(0),
      rconvmass_(0),
      rnewmass_(0),
      axial_strain_GP_(0),
      twist_GP_(0),
      curvature_2_GP_(0),
      curvature_3_GP_(0),
      axial_force_GP_(0),
      torque_GP_(0),
      bending_moment_2_GP_(0),
      bending_moment_3_GP_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(const DRT::ELEMENTS::Beam3k& old)
    : DRT::ELEMENTS::Beam3Base(old),
      useFAD_(old.useFAD_),
      isinit_(old.isinit_),
      T0_(old.T0_),
      T_(old.T_),
      theta0_(old.theta0_),
      Qrefconv_(old.Qrefconv_),
      Qrefnew_(old.Qrefnew_),
      K0_(old.K0_),
      length_(old.length_),
      jacobi_(old.jacobi_),
      jacobi2_(old.jacobi2_),
      jacobi_cp_(old.jacobi_cp_),
      jacobi2_cp_(old.jacobi_cp_),
      rotvec_(old.rotvec_),
      weakkirchhoff_(old.weakkirchhoff_),
      Eint_(old.Eint_),
      Ekin_(old.Ekin_),
      Qconvmass_(old.Qconvmass_),
      Qnewmass_(old.Qnewmass_),
      wconvmass_(old.wconvmass_),
      wnewmass_(old.wnewmass_),
      aconvmass_(old.aconvmass_),
      anewmass_(old.anewmass_),
      amodconvmass_(old.amodconvmass_),
      amodnewmass_(old.amodnewmass_),
      rttconvmass_(old.rttconvmass_),
      rttnewmass_(old.rttnewmass_),
      rttmodconvmass_(old.rttmodconvmass_),
      rttmodnewmass_(old.rttmodnewmass_),
      rtconvmass_(old.rtconvmass_),
      rtnewmass_(old.rtnewmass_),
      rconvmass_(old.rconvmass_),
      rnewmass_(old.rnewmass_),
      axial_strain_GP_(old.axial_strain_GP_),
      twist_GP_(old.twist_GP_),
      curvature_2_GP_(old.curvature_2_GP_),
      curvature_3_GP_(old.curvature_3_GP_),
      axial_force_GP_(old.axial_force_GP_),
      torque_GP_(old.torque_GP_),
      bending_moment_2_GP_(old.bending_moment_2_GP_),
      bending_moment_3_GP_(old.bending_moment_3_GP_)
{
  return;
}
/*--------------------------------------------------------------------------------*
 |  Deep copy this instance of Beam3k and return pointer to it (public) |
 |                                                                    meier 05/12 |
 *--------------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3k::Clone() const
{
  DRT::ELEMENTS::Beam3k* newelement = new DRT::ELEMENTS::Beam3k(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::~Beam3k() { return; }

/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Print(std::ostream& os) const { return; }

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3k::Shape() const
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
    default:
      dserror("Only Line2, Line3 and Line4 elements are implemented.");
      break;
  }

  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Beam3Base::Pack(data);

  // add all class variables
  AddtoPack(data, useFAD_);
  AddtoPack(data, isinit_);
  AddtoPack<3, 1>(data, T0_);
  AddtoPack<3, 1>(data, T_);
  AddtoPack<3, 1>(data, theta0_);
  AddtoPack<4, 1>(data, Qrefconv_);
  AddtoPack<4, 1>(data, Qrefnew_);
  AddtoPack<3, 1>(data, K0_);
  AddtoPack(data, length_);
  AddtoPack(data, jacobi_);
  AddtoPack(data, jacobi2_);
  AddtoPack(data, jacobi_cp_);
  AddtoPack(data, jacobi2_cp_);
  AddtoPack(data, rotvec_);
  AddtoPack(data, weakkirchhoff_);
  AddtoPack(data, Eint_);
  AddtoPack(data, Ekin_);
  AddtoPack<4, 1>(data, Qconvmass_);
  AddtoPack<4, 1>(data, Qnewmass_);
  AddtoPack<3, 1>(data, wconvmass_);
  AddtoPack<3, 1>(data, wnewmass_);
  AddtoPack<3, 1>(data, aconvmass_);
  AddtoPack<3, 1>(data, anewmass_);
  AddtoPack<3, 1>(data, amodconvmass_);
  AddtoPack<3, 1>(data, amodnewmass_);
  AddtoPack<3, 1>(data, rttconvmass_);
  AddtoPack<3, 1>(data, rttnewmass_);
  AddtoPack<3, 1>(data, rttmodconvmass_);
  AddtoPack<3, 1>(data, rttmodnewmass_);
  AddtoPack<3, 1>(data, rtconvmass_);
  AddtoPack<3, 1>(data, rtnewmass_);
  AddtoPack<3, 1>(data, rconvmass_);
  AddtoPack<3, 1>(data, rnewmass_);
  AddtoPack(data, axial_strain_GP_);
  AddtoPack(data, twist_GP_);
  AddtoPack(data, curvature_2_GP_);
  AddtoPack(data, curvature_3_GP_);
  AddtoPack(data, axial_force_GP_);
  AddtoPack(data, torque_GP_);
  AddtoPack(data, bending_moment_2_GP_);
  AddtoPack(data, bending_moment_3_GP_);

  return;
}
/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Unpack(const std::vector<char>& data)
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

  // extract all class variables of beam3k element
  useFAD_ = ExtractInt(position, data);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<3, 1>(position, data, T0_);
  ExtractfromPack<3, 1>(position, data, T_);
  ExtractfromPack<3, 1>(position, data, theta0_);
  ExtractfromPack<4, 1>(position, data, Qrefconv_);
  ExtractfromPack<4, 1>(position, data, Qrefnew_);
  ExtractfromPack<3, 1>(position, data, K0_);
  ExtractfromPack(position, data, length_);
  ExtractfromPack(position, data, jacobi_);
  ExtractfromPack(position, data, jacobi2_);
  ExtractfromPack(position, data, jacobi_cp_);
  ExtractfromPack(position, data, jacobi2_cp_);
  rotvec_ = ExtractInt(position, data);
  weakkirchhoff_ = ExtractInt(position, data);
  ExtractfromPack(position, data, Eint_);
  ExtractfromPack(position, data, Ekin_);
  ExtractfromPack<4, 1>(position, data, Qconvmass_);
  ExtractfromPack<4, 1>(position, data, Qnewmass_);
  ExtractfromPack<3, 1>(position, data, wconvmass_);
  ExtractfromPack<3, 1>(position, data, wnewmass_);
  ExtractfromPack<3, 1>(position, data, aconvmass_);
  ExtractfromPack<3, 1>(position, data, anewmass_);
  ExtractfromPack<3, 1>(position, data, amodconvmass_);
  ExtractfromPack<3, 1>(position, data, amodnewmass_);
  ExtractfromPack<3, 1>(position, data, rttconvmass_);
  ExtractfromPack<3, 1>(position, data, rttnewmass_);
  ExtractfromPack<3, 1>(position, data, rttmodconvmass_);
  ExtractfromPack<3, 1>(position, data, rttmodnewmass_);
  ExtractfromPack<3, 1>(position, data, rtconvmass_);
  ExtractfromPack<3, 1>(position, data, rtnewmass_);
  ExtractfromPack<3, 1>(position, data, rconvmass_);
  ExtractfromPack<3, 1>(position, data, rnewmass_);
  ExtractfromPack(position, data, axial_strain_GP_);
  ExtractfromPack(position, data, twist_GP_);
  ExtractfromPack(position, data, curvature_2_GP_);
  ExtractfromPack(position, data, curvature_3_GP_);
  ExtractfromPack(position, data, axial_force_GP_);
  ExtractfromPack(position, data, torque_GP_);
  ExtractfromPack(position, data, bending_moment_2_GP_);
  ExtractfromPack(position, data, bending_moment_3_GP_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Beam3k::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   meier 01/16|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometry(
    const std::vector<LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  if (weakkirchhoff_)
    SetUpReferenceGeometryWK(xrefe, secondinit);
  else
    SetUpReferenceGeometrySK(xrefe, secondinit);
}

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a weak Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometryWK(
    const std::vector<LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe and rotrefe; note: the
   *isinit_ flag is important for avoiding reinitialization upon restart. However, it should be
   *possible to conduct a second initilization in principle (e.g. for periodic boundary conditions*/

  // TODO: It turned out that the approximation of initially curved centerlines via third order
  // Hermite polynomials by simply setting the nodal positions and tangents at the boundary nodes to
  // the corresponding values of the analytical geometry at these nodes yields a worse approximation
  // than that following from a third order Lagrange polynomial interpolation based on nodal
  // positions that coincide with the analytical values at the four element nodes. This worse
  // approximation appears in terms of a larger error in the jacobian, the derivative of the
  // relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial
  // configurations where the error in the initial geometry representation might dominate the
  // discretization error it can be useful to consider strategies which yield in a better Hermite
  // approximation of the initial geometry e.g. by determining the initial nodal values (or the
  // shape function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if (!isinit_ || secondinit)
  {
    const int nnode = 2;  // number of nodes

    // Calculate the (initial reference triads) = (initial material triads) at the CPs out of the
    // angles theta0_. So far the initial value for the relative angle is set to zero, i.e. material
    // coordinate system and reference system in the reference configuration coincidence.
    std::vector<LINALG::Matrix<3, 3>> Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node], Gref[node]);
    }

    T0_.resize(2);
    T_.resize(2);
    // write initial nodal tangents in extra vector
    for (int i = 0; i < 3; i++)
    {
      (T0_[0])(i) = (Gref[0])(i, 0);
      (T0_[1])(i) = (Gref[1])(i, 0);
      (T_[0])(i) = (Gref[0])(i, 0);
      (T_[1])(i) = (Gref[1])(i, 0);
    }

    // Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints =
        DRT::UTILS::IntegrationPoints1D(DRT::UTILS::MYGAUSSRULEBEAM3K);

    // Vector holding angle theta of triads
    std::vector<LINALG::Matrix<3, 1>> theta_cp;
    theta_cp.resize(BEAM3K_COLLOCATION_POINTS);
    LINALG::Matrix<3, 1> theta(true);
    LINALG::Matrix<3, 1> theta_s(true);
    LINALG::Matrix<3, 3> triad_mat(true);  // material triad at gp

    // Resize vectors for storage of time integration quantities
    ResizeClassVariables(gausspoints.nquad);

    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_GP_.resize(gausspoints.nquad);
    std::fill(axial_strain_GP_.begin(), axial_strain_GP_.end(), 0.0);
    twist_GP_.resize(gausspoints.nquad);
    std::fill(twist_GP_.begin(), twist_GP_.end(), 0.0);
    curvature_2_GP_.resize(gausspoints.nquad);
    std::fill(curvature_2_GP_.begin(), curvature_2_GP_.end(), 0.0);
    curvature_3_GP_.resize(gausspoints.nquad);
    std::fill(curvature_3_GP_.begin(), curvature_3_GP_.end(), 0.0);

    axial_force_GP_.resize(gausspoints.nquad);
    std::fill(axial_force_GP_.begin(), axial_force_GP_.end(), 0.0);
    torque_GP_.resize(gausspoints.nquad);
    std::fill(torque_GP_.begin(), torque_GP_.end(), 0.0);
    bending_moment_2_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_2_GP_.begin(), bending_moment_2_GP_.end(), 0.0);
    bending_moment_3_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_3_GP_.begin(), bending_moment_3_GP_.end(), 0.0);


    // calculate the length of the element via Newton iteration
    Calculate_length(xrefe, T0_, LENGTHCALCNEWTONTOL);

    // Matrices to store the function values of the Lagrange shape functions used to interpolate
    // theta
    LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i;
    LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i_xi;

    // Matrices to store the (derivative of) the Hermite shape functions
    LINALG::Matrix<1, 2 * nnode> N_i;
    LINALG::Matrix<1, 2 * nnode> N_i_xi;

    // Matrices to store r and dr/dxi
    LINALG::Matrix<3, 1> r;
    LINALG::Matrix<3, 1> r_xi;

    // storage index for collocation point
    unsigned int ind = 0;

    // Calculate initial material triads at the collocation points
    for (unsigned int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      // colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi = (double)node / (BEAM3K_COLLOCATION_POINTS - 1) * 2 - 1.0;

      // Get values of shape functions
      N_i_xi.Clear();
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);

      // Determine storage position for the node colpt
      ind = LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

      // current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();

      for (int i = 0; i < 3; i++)
      {
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + T0_[0](i) * N_i_xi(1) +
                   T0_[1](i) * N_i_xi(3);
      }

      jacobi_cp_[ind] = r_xi.Norm2();

      // rotate (initial reference triad) = (initial material triad) at the interior CPs on
      // tangential line resulting from the Hermite interpolation. This is necessary for initially
      // curved geometries for which the Hermite interpolation does not deliver the exact tangent
      // values of the analytical representation of the initial curve at the CPs.
      if (ind > 1)  // only for internal CPs
      {
        LINALG::Matrix<3, 3> G_aux(true);
        LARGEROTATIONS::CalculateSRTriads<double>(r_xi, Gref[ind], G_aux);
        // rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind] = G_aux;
        LARGEROTATIONS::triadtoquaternion(G_aux, Qrefconv_[ind]);
        LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind], theta0_[ind]);
      }
      else
      {
        LARGEROTATIONS::triadtoquaternion(Gref[ind], Qrefconv_[ind]);
      }
      Qrefnew_[ind] = Qrefconv_[ind];
    }

    // SETUP INTERPOLATION via calculation of difference angle
    for (int colpt = 0; colpt < BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      theta_cp[colpt].Clear();
      LARGEROTATIONS::triadtoangleright(theta_cp[colpt], Gref[REFERENCE_NODE], Gref[colpt]);
    }

    // Loop through all GPs and computation of all relevant values at each gp
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      // Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i_xi.Clear();
      N_i.Clear();
      DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      DRT::UTILS::shape_function_1D_deriv1(L_i_xi, xi, Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);

      // current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      theta.Clear();
      theta_s.Clear();

      for (int i = 0; i < 3; i++)
      {
        r(i) +=
            xrefe[0](i) * N_i(0) + xrefe[1](i) * N_i(2) + T0_[0](i) * N_i(1) + T0_[1](i) * N_i(3);
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + T0_[0](i) * N_i_xi(1) +
                   T0_[1](i) * N_i_xi(3);
      }

      // calculate jacobi jacobi_=|r'_0|
      jacobi_[numgp] = r_xi.Norm2();

      // calculate interpolated angle
      for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          theta(j) += L_i(i) * theta_cp[i](j);
          theta_s(j) += L_i_xi(i) * theta_cp[i](j) / jacobi_[numgp];
        }
      }
      computestrain(theta, theta_s, K0_[numgp]);

      triad_mat.Clear();
      LARGEROTATIONS::angletotriad(theta, Gref[REFERENCE_NODE], triad_mat);
      SetInitialDynamicClassVariables(numgp, triad_mat, r);
    }

    isinit_ = true;
  }
}

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a strong Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK(
    const std::vector<LINALG::Matrix<3, 1>>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe and rotrefe; note: the
   *isinit_ flag is important for avoiding reinitialization upon restart. However, it should be
   *possible to conduct a second initilization in principle (e.g. for periodic boundary conditions*/

  // TODO: It turned out that the approximation of initially curved centerlines via third order
  // Hermite polynomials by simply setting the nodal positions and tangents at the boundary nodes to
  // the corresponding values of the analytical geometry at these nodes yields a worse approximation
  // than that following from a third order Lagrange polynomial interpolation based on nodal
  // positions that coincide with the analytical values at the four element nodes. This worse
  // approximation appears in terms of a larger error in the jacobian, the derivative of the
  // relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial
  // configurations where the error in the initial geometry representation might dominate the
  // discretization error it can be useful to consider strategies which yield in a better Hermite
  // approximation of the initial geometry e.g. by determining the initial nodal values (or the
  // shape function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if (!isinit_ || secondinit)
  {
    const int nnode = 2;  // number of nodes

    // Calculate the (initial reference triads) = (initial material triads) at the CPs out of the
    // angles theta0_. So far the initial value for the relative angle is set to zero, i.e. material
    // coordinate system and reference system in the reference configuration coincidence.
    std::vector<LINALG::Matrix<3, 3>> Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node], Gref[node]);
    }

    T0_.resize(2);
    T_.resize(2);
    // write initial nodal tangents in extra vector
    for (int i = 0; i < 3; i++)
    {
      (T0_[0])(i) = (Gref[0])(i, 0);
      (T0_[1])(i) = (Gref[1])(i, 0);
      (T_[0])(i) = (Gref[0])(i, 0);
      (T_[1])(i) = (Gref[1])(i, 0);
    }

    // Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints =
        DRT::UTILS::IntegrationPoints1D(DRT::UTILS::MYGAUSSRULEBEAM3K);

    // Vector holding angle theta of triads
    std::vector<double> phi_cp;
    phi_cp.resize(BEAM3K_COLLOCATION_POINTS);
    double phi(true);
    double phi_s(true);
    LINALG::Matrix<3, 3> triad_mat(true);  // material triad at gp

    // Resize vectors for storage of time integration quantities
    ResizeClassVariables(gausspoints.nquad);

    // assure correct size of strain and stress resultant class variables and fill them
    // with zeros (by definition, the reference configuration is undeformed and stress-free)
    axial_strain_GP_.resize(gausspoints.nquad);
    std::fill(axial_strain_GP_.begin(), axial_strain_GP_.end(), 0.0);
    twist_GP_.resize(gausspoints.nquad);
    std::fill(twist_GP_.begin(), twist_GP_.end(), 0.0);
    curvature_2_GP_.resize(gausspoints.nquad);
    std::fill(curvature_2_GP_.begin(), curvature_2_GP_.end(), 0.0);
    curvature_3_GP_.resize(gausspoints.nquad);
    std::fill(curvature_3_GP_.begin(), curvature_3_GP_.end(), 0.0);

    axial_force_GP_.resize(gausspoints.nquad);
    std::fill(axial_force_GP_.begin(), axial_force_GP_.end(), 0.0);
    torque_GP_.resize(gausspoints.nquad);
    std::fill(torque_GP_.begin(), torque_GP_.end(), 0.0);
    bending_moment_2_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_2_GP_.begin(), bending_moment_2_GP_.end(), 0.0);
    bending_moment_3_GP_.resize(gausspoints.nquad);
    std::fill(bending_moment_3_GP_.begin(), bending_moment_3_GP_.end(), 0.0);


    // calculate the length of the element via Newton iteration
    Calculate_length(xrefe, T0_, LENGTHCALCNEWTONTOL);

    // Matrices to store the function values of the Lagrange shape functions used to interpolate
    // theta
    LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i;
    LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS> L_i_xi;

    // Matrices to store the (derivative of) the Hermite shape functions
    LINALG::Matrix<1, 2 * nnode> N_i;
    LINALG::Matrix<1, 2 * nnode> N_i_xi;
    LINALG::Matrix<1, 2 * nnode> N_i_xixi;

    // Matrices to store r, dr/dxi and d^2r/dxi^2
    LINALG::Matrix<3, 1> r;
    LINALG::Matrix<3, 1> r_xi;
    LINALG::Matrix<3, 1> r_xixi;
    LINALG::Matrix<3, 1> r_s;
    LINALG::Matrix<3, 1> r_ss;
    // centerline curvature
    LINALG::Matrix<3, 1> kappacl;

    // storage index for collocation point
    unsigned int ind = 0;

    // Calculate initial material triads at the collocation points
    for (unsigned int node = 0; node < BEAM3K_COLLOCATION_POINTS; node++)
    {
      // colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi = (double)node / (BEAM3K_COLLOCATION_POINTS - 1) * 2 - 1.0;

      // Get values of shape functions
      L_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();
      DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi, xi, length_, line2);

      // Determine storage position for the node colpt
      ind = LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

      // current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();
      r_xixi.Clear();

      for (int i = 0; i < 3; i++)
      {
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + T0_[0](i) * N_i_xi(1) +
                   T0_[1](i) * N_i_xi(3);
        r_xixi(i) += xrefe[0](i) * N_i_xixi(0) + xrefe[1](i) * N_i_xixi(2) +
                     T0_[0](i) * N_i_xixi(1) + T0_[1](i) * N_i_xixi(3);
      }

      // calculate jacobi_=||r'_0|| and jacobi2_=r'_0^T r''_0
      jacobi_cp_[ind] = r_xi.Norm2();
      jacobi2_cp_[ind] = r_xi.Dot(r_xixi);

      // rotate (initial reference triad) = (initial material triad) at the interior CPs on
      // tangential line resulting from the Hermite interpolation. This is necessary for initially
      // curved geometries for which the Hermite interpolation does not deliver the exact tangent
      // values of the analytical representation of the initial curve at the CPs.
      if (ind > 1)  // only for internal CPs
      {
        LINALG::Matrix<3, 3> G_aux(true);
        LARGEROTATIONS::CalculateSRTriads<double>(r_xi, Gref[ind], G_aux);
        // rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind] = G_aux;
        LARGEROTATIONS::triadtoquaternion(G_aux, Qrefconv_[ind]);
        LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind], theta0_[ind]);
      }
      else
      {
        LARGEROTATIONS::triadtoquaternion(Gref[ind], Qrefconv_[ind]);
      }
      Qrefnew_[ind] = Qrefconv_[ind];
    }  //(int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)

    // SETUP INTERPOLATION via calculation of difference angle
    for (int colpt = 0; colpt < BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      LINALG::Matrix<3, 3> Lambdabarref(true);
      LINALG::Matrix<3, 1> tangentref(true);
      LINALG::Matrix<3, 1> phivec(true);
      for (int i = 0; i < 3; i++)
      {
        tangentref(i) = Gref[colpt](i, 0);
      }
      LARGEROTATIONS::CalculateSRTriads<double>(tangentref, Gref[REFERENCE_NODE], Lambdabarref);
      LARGEROTATIONS::triadtoangleleft(phivec, Lambdabarref, Gref[colpt]);
      phi_cp[colpt] = 0.0;
      for (int i = 0; i < 3; i++)
      {
        phi_cp[colpt] += tangentref(i) * phivec(i);
      }
    }

    // Loop through all GPs and computation of all relevant values at each gp
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      // Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();

      DRT::UTILS::shape_function_1D(L_i, xi, Shape());
      DRT::UTILS::shape_function_1D_deriv1(L_i_xi, xi, Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
      DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi, xi, length_, line2);
      DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);

      // current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      r_xixi.Clear();
      r_s.Clear();
      r_ss.Clear();
      kappacl.Clear();
      phi = 0.0;
      phi_s = 0.0;

      for (int i = 0; i < 3; i++)
      {
        r(i) +=
            xrefe[0](i) * N_i(0) + xrefe[1](i) * N_i(2) + T0_[0](i) * N_i(1) + T0_[1](i) * N_i(3);
        r_xi(i) += xrefe[0](i) * N_i_xi(0) + xrefe[1](i) * N_i_xi(2) + T0_[0](i) * N_i_xi(1) +
                   T0_[1](i) * N_i_xi(3);
        r_xixi(i) += xrefe[0](i) * N_i_xixi(0) + xrefe[1](i) * N_i_xixi(2) +
                     T0_[0](i) * N_i_xixi(1) + T0_[1](i) * N_i_xixi(3);
      }

      // calculate jacobi_=||r'_0|| and jacobi2_=r'_0^T r''_0
      jacobi_[numgp] = r_xi.Norm2();
      jacobi2_[ind] = r_xi.Dot(r_xixi);

      // calculate interpolated angle
      for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
      {
        phi += L_i(i) * phi_cp[i];
        phi_s += L_i_xi(i) * phi_cp[i] / jacobi_[numgp];
      }

      // calculate derivatives in s
      r_s = r_xi;
      r_s.Scale(1 / jacobi_[numgp]);
      for (int i = 0; i < 3; i++)
      {
        r_ss(i) = r_xixi(i) / pow(jacobi_[numgp], 2.0) -
                  r_xi(i) * jacobi2_[numgp] / pow(jacobi_[numgp], 4.0);
      }

      triad_mat.Clear();
      ComputeTriadSK(phi, r_s, Gref[REFERENCE_NODE], triad_mat);
      Calculate_clcurvature(r_s, r_ss, kappacl);

      computestrainSK(phi_s, kappacl, Gref[REFERENCE_NODE], triad_mat, K0_[numgp]);

      SetInitialDynamicClassVariables(numgp, triad_mat, r);
    }

    isinit_ = true;
  }  // if(!isinit_)

}  // DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK()

/*--------------------------------------------------------------------------------------------*
 |  Calculates the element length via a Newton Iteration: f(l)=l-int(|N'd|)dxi=0   meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Calculate_length(const std::vector<LINALG::Matrix<3, 1>>& xrefe,
    const std::vector<LINALG::Matrix<3, 1>>& trefe, double tolerance)
{
  const int nnode = 2;  // number of nodes
  const int vnode = 2;  // interpolated values per node (2: value + derivative of value)

  // Get integration points for exact integration
  // DRT::UTILS::IntegrationPoints1D gausspoints =
  //    DRT::UTILS::IntegrationPoints1D(DRT::UTILS::MYGAUSSRULEBEAM3K);
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

  // inital value for iteration
  {
    LINALG::Matrix<3, 1> tempvec;
    tempvec.Clear();
    for (int i = 0; i < 3; i++)
    {
      tempvec(i) = xrefe[1](i) - xrefe[0](i);
    }
    length_ = tempvec.Norm2();
  }

  // Matrices to store the function values of the shape functions
  LINALG::Matrix<1, nnode * vnode> shapefuncderiv;

  shapefuncderiv.Clear();

  // current value of the derivative at the GP
  LINALG::Matrix<3, 1> r_xi;

  while (std::fabs(res) > tolerance)
  {
    int_length = 0;
    deriv_length = 0;
    // Loop through all GPs and computation of the length and the derivative of the length
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      deriv_int = 0;
      // Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      // Get derivatives of the shape functions
      DRT::UTILS::shape_function_hermite_1D_deriv1(shapefuncderiv, xi, length_, line2);

      // integral of the length
      r_xi.Clear();
      deriv_int = 0;
      for (int i = 0; i < 3; i++)
      {
        r_xi(i) += xrefe[0](i) * shapefuncderiv(0) + xrefe[1](i) * shapefuncderiv(2) +
                   trefe[0](i) * shapefuncderiv(1) + trefe[1](i) * shapefuncderiv(3);
      }
      int_length += gausspoints.qwgt[numgp] * r_xi.Norm2();

      // derivative of the integral of the length at GP
      for (int i = 0; i < 3; i++)
      {
        deriv_int += (trefe[0](i) * shapefuncderiv(1) / length_ +
                         trefe[1](i) * shapefuncderiv(3) / length_) *
                     r_xi(i);
      }
      deriv_length += gausspoints.qwgt[numgp] * deriv_int / r_xi.Norm2();
    }
    // cout << endl << "Länge:" << length_ << "\tIntegral:" << int_length << "\tAbleitung:" <<
    // deriv_length << endl;
    res = length_ - int_length;
    // Update
    length_ = length_ - res / (1 - deriv_length);  // the derivative of f(l)=l-int(|N'd|)dxi=0 is
                                                   // f'(l)=1-int(d/dl(|N'd|))dxi
  }

  //  //*************************************begin: Determine optimal Hermite constant for circle
  //  segment****************************** if(fabs(trefe[0].Norm2()-1.0)>1.0e-12 or
  //  fabs(trefe[1].Norm2()-1.0)>1.0e-12)
  //    dserror("Tangents have to be unit vectors!");
  //
  //  double copt=0.0;
  //  double int1=0.0;
  //  double int2=0.0;
  //
  //  LINALG::Matrix<3,1> v0(true);
  //  LINALG::Matrix<3,1> n0(true);
  //  LINALG::Matrix<3,3> Sv0(true);
  //  LINALG::Matrix<3,3> Sn0(true);
  //  LINALG::Matrix<3,1> r0(true);
  //  LINALG::Matrix<3,1> R0vec(true);
  //  LINALG::Matrix<3,3> triad0(true);
  //  double R0 = 0.0;
  //  double alpha0 = 0.0;
  //  double d0 = 0.0;
  //  LINALG::Matrix<1,1> scalarproduct10(true);
  //  LINALG::Matrix<1,1> scalarproduct20(true);
  //
  //  v0.Update(1.0,xrefe[1],1.0);
  //  v0.Update(-1.0,xrefe[0],1.0);
  //  d0 = v0.Norm2();
  //  v0.Scale(1.0/d0);
  //
  //  LARGEROTATIONS::computespin(Sv0,v0);
  //  n0.Multiply(Sv0,trefe[0]);
  //  n0.Scale(1.0/n0.Norm2());
  //  LARGEROTATIONS::computespin(Sn0,n0);
  //  scalarproduct10.MultiplyTN(n0,trefe[1]);
  //  if(scalarproduct10.Norm2()>1.0e-12)
  //    dserror("Only plane geometries possible at this point!");
  //
  //  scalarproduct20.MultiplyTN(v0,trefe[0]);
  //  alpha0 = acos(scalarproduct20(0,0));
  //  R0=d0/(2*sin(alpha0));
  //  R0vec.Multiply(Sn0,trefe[0]);
  //  R0vec.Scale(1/R0vec.Norm2());
  //  R0vec.Scale(R0);
  //  r0.Update(1.0,xrefe[0],0.0);
  //  r0.Update(-1.0,R0vec,1.0);
  //
  //  for(int i=0;i<3;i++)
  //  {
  //    triad0(i,0)=R0vec(i);
  //    triad0(i,1)=n0(i);
  //    triad0(i,2)=trefe[0](i);
  //  }
  //
  ////  std::cout << "xrefe[0]: " << xrefe[0] << std::endl;
  ////  std::cout << "xrefe[1]: " << xrefe[1] << std::endl;
  ////  std::cout << "trefe[0]: " << trefe[0] << std::endl;
  ////  std::cout << "trefe[1]: " << trefe[1] << std::endl;
  ////  std::cout << "alpha0: " << alpha0 << std::endl;
  ////  std::cout << "r0: " << r0 << std::endl;
  ////  std::cout << "R0vec: " << R0vec << std::endl;
  ////  std::cout << "d0: " << d0 << std::endl;
  ////  std::cout << "R0: " << R0 << std::endl;
  //
  //  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  //  {
  //    //Get position xi of GP
  //    const double xi = gausspoints.qxg[numgp][0];
  //    const double wgt = gausspoints.qwgt[numgp];
  //    LINALG::Matrix<3,1> a(true);
  //    LINALG::Matrix<3,1> b(true);
  //    LINALG::Matrix<3,1> r(true);
  //    LINALG::Matrix<3,1> r_hermite(true);
  //    LINALG::Matrix<3,1> auxvec(true);
  //    LINALG::Matrix<1,1> auxscal(true);
  //
  //    //Get derivatives of the shape functions
  //    LINALG::Matrix<1,nnode*vnode> shapefunc(true);
  //    DRT::UTILS::shape_function_hermite_1D(shapefunc,xi,1.0,line2);
  //
  //    for (int i=0; i<3; i++)
  //    {
  //      b(i)+=trefe[0](i)*shapefunc(1)+trefe[1](i)*shapefunc(3);
  //      a(i)+=xrefe[0](i)*shapefunc(0)+xrefe[1](i)*shapefunc(2);
  //      r_hermite(i)+=xrefe[0](i)*shapefunc(0)+xrefe[1](i)*shapefunc(2)+trefe[0](i)*length_*shapefunc(1)+trefe[1](i)*length_*shapefunc(3);
  //    }
  //
  //    LINALG::Matrix<3,3> triad(true);
  //    LINALG::Matrix<3,3> triadalpha(true);
  //    LINALG::Matrix<3,1> Rvec(true);
  //    LINALG::Matrix<3,1> alphavec(true);
  //
  //    alphavec.Update(-(1.0+xi)*alpha0,n0,0.0);
  //    LARGEROTATIONS::angletotriad(alphavec,triadalpha);
  //    triad.Multiply(triadalpha,triad0);
  //    for(int i=0;i<3;i++)
  //    {
  //      Rvec(i)=triad(i,0);
  //    }
  //    r.Update(1.0,r0,0.0);
  //    r.Update(1.0,Rvec,1.0);
  //
  ////    std::cout << "r: " << r << std::endl;
  ////    std::cout << "r_hermite: " << r_hermite << std::endl;
  ////
  ////    std::cout << "alphavec: " << alphavec << std::endl;
  ////    std::cout << "Rvec: " << Rvec << std::endl;
  //
  //
  //    auxvec.Update(1.0,r,0.0);
  //    auxvec.Update(-1.0,a,1.0);
  //    auxscal.MultiplyTN(b,auxvec);
  //
  //    int1+= auxscal(0,0)*wgt;
  //
  //    auxscal=0.0;
  //    auxscal.MultiplyTN(b,b);
  //    int2+= auxscal(0,0)*wgt;
  //  }
  //
  //  copt=int1/int2;
  //  double length_approx=d0;
  //  double length_analyt=2*alpha0*R0;
  //  std::cout << "copt length_ length_approx length_analyt" << std::endl;
  //
  //  std::cout << std::setprecision(16) << copt << " " << length_ << " " << length_approx << " " <<
  //  length_analyt << std::endl;
  //  //*************************************end: Determine optimal Hermite constant for circle
  //  segment******************************

  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3k::GetJacobiFacAtXi(const double& xi) const
{
  const int nnode = 2;

  // Matrices to store the the Hermite shape function derivative values
  LINALG::Matrix<1, 2 * nnode> N_i_xi;
  DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);

  // jacobi = ds/dxi = ||r'_0||
  LINALG::Matrix<3, 1> r_xi;

  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    r_xi(dim) += Nodes()[0]->X()[dim] * N_i_xi(0) + Nodes()[1]->X()[dim] * N_i_xi(2) +
                 T0_[0](dim) * N_i_xi(1) + T0_[1](dim) * N_i_xi(3);
  }

  return r_xi.Norm2();
}

/*----------------------------------------------------------------------------------------------------------*
 | Get position vector at xi for given nodal displacements popp 02/16|
 *----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> DRT::ELEMENTS::Beam3k::GetPos(
    const double& xi, const LINALG::Matrix<12, 1>& disp_totlag_centerline) const
{
  // note: this method expects the absolute ("total Lagrangean") values for positions and tangents
  // of both centerline nodes (local numbering 0 and 1)
  LINALG::Matrix<3, 1> r(true);
  LINALG::Matrix<1, 4> N_i(true);

  DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);

  this->Calc_r<2, 2, double>(disp_totlag_centerline, N_i, r);

  return (r);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetPosAtXi(
    LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const
{
  /* we expect the (centerline) displacement state vector (positions and tangents at boundary nodes)
   * here; for flexibility, we also accept complete DoF vector of this element (2*7+1 DoFs) and do
   * the rest automatically here NOTE: the latter option is favorable in case of rotvec_==true: in
   * case of ROTVEC=true, we need 3 positions, 3 absolute rotvec DOFs AND the length of tangent
   * vector (7th DoF) */
  LINALG::Matrix<12, 1> disp_totlag_centerline(true);

  if (disp.size() == 15)
  {
    // in this case, we need to "add" reference values first, because if rotvec_==true,
    // we can extract tangent vectors only from total rotation vectors
    LINALG::Matrix<15, 1, double> disp_totlag(&disp[0]);
    AddRefValuesDisp<2, double>(disp_totlag);
    this->ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
        disp_totlag, disp_totlag_centerline);
  }
  else if (disp.size() == 12)
  {
    /* in this case, we expect the position and tangent DOF values for both boundary nodes;
     * for rotvec_==true, the tangents are NOT nodal DOFs, so they need to be pre-calculated
     * before calling this method */
    disp_totlag_centerline = LINALG::Matrix<12, 1>(&disp[0]);
    AddRefValuesDispCenterline<2, 2, double>(disp_totlag_centerline);
  }
  else
  {
    dserror(
        "size mismatch: expected either 12 values for disp_centerline or 15 for "
        "full element disp vector and got %d",
        disp.size());
  }

  pos = this->GetPos(xi, disp_totlag_centerline);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetTriadAtXi(
    LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_)
    dserror("method GetTriadAtXi is limited to WK so far! extend to SK if needed");

  if (disp.size() != 15)
    dserror("size mismatch: expected 15 values for element disp vector and got %d", disp.size());


  // Dof vector in total Lagrangian style, i.e. "displacement + reference values"
  LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double> disp_totlag(true);

  UpdateDispTotlag<2, double>(disp, disp_totlag);

  // material triads at collocation points
  std::vector<LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<3, 3, double>(true));
  LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double> dummy(true);
  std::vector<LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<4, 1>(true));

  // Todo @grill:
  //    this method uses Qrefconv_[node] as reference triads; so we can only call this method if
  //    those were calculated correctly in a preceding time step. (be careful in post-processing!)
  UpdateNodalVariables<2, double>(disp_totlag, dummy, triad_mat_cp,
      Qref_dummy);  // Todo @grill split/adapt method and avoid dummy variables

  // create object of triad interpolation scheme
  Teuchos::RCP<
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS, double>>
      triad_interpolation_scheme_ptr = Teuchos::rcp(
          new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>());

  // reset scheme with nodal triads
  triad_interpolation_scheme_ptr->Reset(triad_mat_cp);

  triad_interpolation_scheme_ptr->GetInterpolatedTriadAtXi(triad, xi);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetScaledSecondAndThirdBaseVectorAtXi(const double& xi,
    const std::vector<double>& disp, LINALG::Matrix<3, 2>& scaledbasevectors) const
{
  // Todo @grill: delete or update this method. kept it for now as we might need it soon
  dserror(
      "Beam3k: method GetScaledSecondAndThirdBaseVectorAtXi is deprecated for now because it "
      "was only valid for the implicit assumption of a rectangular cross-section! If need be, "
      "generalize the definition of cross-section shape and dimensions and adapt this method. "
      "For now, the cross-section shape and dimensions are only required explicitly if beam "
      "interactions (contact, potentials, drag in background fluid) are to be evaluated. In this "
      "case, we only support and hence assume a circular cross-section with a radius which is "
      "either "
      "explicitly specified in the material definition line of the input file or per default "
      "computed "
      "from the area moment of inertia Iyy.");

  //  LINALG::Matrix<3,3> triad(true);
  //
  //  GetTriadAtXi(triad,xi,disp);
  //
  //  // ToDo careful, this is a hack for rectangular cross-sections ?!?
  //  double Ly=std::pow(pow(12.0*Izz_,3)/(12.0*Iyy_),1.0/8.0);
  //  double Lz=12.0*Izz_/pow(Ly,3);
  //
  //  for(int i=0;i<3;i++)
  //  {
  //    scaledbasevectors(i,0) = 0.5 * Ly * triad(i,1);
  //    scaledbasevectors(i,1) = 0.5 * Lz * triad(i,2);
  //  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetGeneralizedInterpolationMatrixVariationsAtXi(
    LINALG::SerialDenseMatrix& Ivar, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");

  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)Ivar.M() != 2 * ndim or (unsigned int) Ivar.N() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6, numdof, Ivar.M(), Ivar.N());

  Ivar.Zero();

  // *******************************************************************************
  // concerning interpolation of variation of CENTERLINE POSITION \vardelta r:
  // *******************************************************************************

  LINALG::Matrix<ndim, numdof, double> N(true);
  LINALG::Matrix<1, vpernode * nnodecl, double> N_i(true);


  DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);
  AssembleShapefunctionsN(N_i, N);

  // this part is associated with the variation of the centerline position
  // (first three rows of Ivar)
  for (unsigned int irow = 0; irow < N.M(); ++irow)
    for (unsigned int icol = 0; icol < N.N(); ++icol) Ivar(irow, icol) += N(irow, icol);

  // *******************************************************************************
  // concerning interpolation of variation of CENTERLINE ORIENTATION \vardelta \theta:
  // *******************************************************************************

  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  LINALG::Matrix<ndim, numdof, double> N_s(true);
  LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);

  LINALG::Matrix<1, numdof, double> L(true);
  LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<LINALG::Matrix<3, 3, double>> triad_dummy(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<3, 3, double>(true));
  std::vector<LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<4, 1>(true));


  LINALG::Matrix<3, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                    // ||r'||


  std::vector<LINALG::Matrix<numdof, ndim, double>> v_thetaperp_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<numdof, ndim, double>(true));
  std::vector<LINALG::Matrix<numdof, ndim, double>> v_thetapar_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<numdof, ndim, double>(true));


  // re-interpolated spin vector variation: v_theta_bar
  LINALG::Matrix<numdof, ndim, double> v_theta_bar(true);


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_dummy,
      Qref_dummy);  // Todo @grill: make this nicer and avoid dummies !


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);


    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);


    // Calculation of r' at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = FADUTILS::Norm(r_s);


    Calc_v_thetaperp<nnodecl>(v_thetaperp_cp[ind], N_s, r_s, abs_r_s);

    Calc_v_thetapartheta<nnodecl>(v_thetapar_cp[ind], L, r_s, abs_r_s);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************

  L_i.Clear();
  DRT::UTILS::shape_function_1D(L_i, xi, Shape());

  v_theta_bar.Clear();
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    v_theta_bar.Update(L_i(icp), v_thetaperp_cp[icp], 1.0);
    v_theta_bar.Update(L_i(icp), v_thetapar_cp[icp], 1.0);
  }

  // finally assemble the generalized interpolation matrix for the variations Ivar
  // *******************************************************************************

  // this part is associated with the increment of the centerline orientation
  // (expressed as rotation vector theta)
  // (rows 4-6 of Iinc)
  // note: we need the transposed of v_theta_bar (rows <-> cols)
  for (unsigned int irow = 0; irow < v_theta_bar.N(); ++irow)
    for (unsigned int icol = 0; icol < v_theta_bar.M(); ++icol)
      Ivar(ndim + irow, icol) += v_theta_bar(icol, irow);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
    LINALG::SerialDenseMatrix& stiffmat, const double& xi, const std::vector<double>& disp,
    const LINALG::SerialDenseVector& force) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");


  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)stiffmat.M() != numdof or (unsigned int) stiffmat.N() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", numdof, numdof, stiffmat.M(),
        stiffmat.N());

  stiffmat.Zero();

  // create an auxiliary fixed size matrix and set as a view to original data in stiffmat
  LINALG::Matrix<numdof, numdof, double> stiffmat_fixedsize(stiffmat, true);


  LINALG::Matrix<ndim, 1, double> moment(&(force(3)));

  LINALG::Matrix<ndim, ndim, double> S_of_moment(true);
  LARGEROTATIONS::computespin<double>(S_of_moment, moment);


  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  LINALG::Matrix<ndim, numdof, double> N_s(true);
  LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);


  LINALG::Matrix<1, numdof, double> L(true);
  LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<3, 3, double>(true));
  std::vector<LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<4, 1>(true));


  LINALG::Matrix<ndim, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                       // ||r'||

  // first base vector at CP
  LINALG::Matrix<ndim, 1, double> g_1_cp(true);


  std::vector<LINALG::Matrix<numdof, numdof, double>> lin_v_thetaperp_moment_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<numdof, numdof, double>(true));

  std::vector<LINALG::Matrix<numdof, numdof, double>> lin_v_thetapar_moment_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<numdof, numdof, double>(true));


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_mat_cp,
      Qref_dummy);  // Todo make this nicer and avoid dummies!


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);


    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);


    // Calculation of first base vector at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = FADUTILS::Norm(r_s);

    g_1_cp.Clear();
    g_1_cp.Update(std::pow(abs_r_s, -1.0), r_s);


    Calc_lin_v_thetaperp_moment<nnodecl>(
        lin_v_thetaperp_moment_cp[ind], N_s, g_1_cp, abs_r_s, S_of_moment);

    Calc_lin_v_thetapar_moment<nnodecl>(
        lin_v_thetapar_moment_cp[ind], L, N_s, g_1_cp, abs_r_s, moment);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************
  L_i.Clear();
  DRT::UTILS::shape_function_1D(L_i, xi, Shape());

  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    stiffmat_fixedsize.Update(L_i(icp), lin_v_thetaperp_moment_cp[ind], 1.0);
    stiffmat_fixedsize.Update(L_i(icp), lin_v_thetapar_moment_cp[ind], 1.0);
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetGeneralizedInterpolationMatrixIncrementsAtXi(
    LINALG::SerialDenseMatrix& Iinc, const double& xi, const std::vector<double>& disp) const
{
  if (not weakkirchhoff_) dserror("method is limited to WK so far! extend to SK if needed");

  if (rotvec_)
    dserror("method is limited to tangent-based parametrization! extend to rotvec if needed");


  const unsigned int ndim = 3;
  const unsigned int nnodecl = 2;
  const unsigned int vpernode = 2;
  const unsigned int numdof = ndim * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS;

  // safety check
  if ((unsigned int)Iinc.M() != 2 * ndim or (unsigned int) Iinc.N() != numdof)
    dserror("size mismatch! expected %dx%d matrix and got %dx%d", 6, numdof, Iinc.M(), Iinc.N());

  Iinc.Zero();


  // *******************************************************************************
  // concerning interpolation of increment of CENTERLINE POSITION \Delta r:
  // *******************************************************************************

  LINALG::Matrix<ndim, numdof, double> N(true);
  LINALG::Matrix<1, vpernode * nnodecl, double> N_i(true);


  DRT::UTILS::shape_function_hermite_1D(N_i, xi, length_, line2);
  AssembleShapefunctionsN(N_i, N);

  // this part is associated with the increment of the centerline position
  // (first three rows of Iinc)
  for (unsigned int irow = 0; irow < N.M(); ++irow)
    for (unsigned int icol = 0; icol < N.N(); ++icol) Iinc(irow, icol) += N(irow, icol);


  // *******************************************************************************
  // concerning interpolation of increment of CENTERLINE ORIENTATION \Delta \theta:
  // *******************************************************************************

  // define and initialize variables
  // *******************************************************************************
  // position index where CP quantities have to be stored (according to numbering convention)
  unsigned int ind = 0;
  double xi_cp = 0.0;


  LINALG::Matrix<ndim, numdof, double> N_s(true);
  LINALG::Matrix<1, vpernode * nnodecl, double> N_i_xi(true);

  LINALG::Matrix<1, numdof, double> L(true);
  LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, double> L_i(true);


  LINALG::Matrix<numdof, 1, double> disp_totlag(true);
  LINALG::Matrix<numdof, 1, double> disp_totlag_centerline(true);

  std::vector<LINALG::Matrix<3, 3, double>> triad_mat_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<3, 3, double>(true));
  std::vector<LINALG::Matrix<4, 1>> Qref_dummy(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<4, 1>(true));


  LINALG::Matrix<3, 1, double> r_s(true);  // r' vector
  double abs_r_s = 0.0;                    // ||r'||

  // first base vector at CP
  LINALG::Matrix<3, 1, double> g_1_cp(true);
  // first base vector of the triad, from which the new intermediate triad is obtained via
  // smallest rotation (SR); this triad is arbitrary, but we choose the intermediate triad
  // of the last time step; see Dissertation Meier, p.25
  LINALG::Matrix<3, 1, double> g_1_cp_bar(true);


  LINALG::Matrix<ndim, numdof, double> lin_theta_perp_cp(true), lin_theta_par_cp(true);

  // lin_theta_cp = lin_theta_perp_cp + lin_theta_par_cp
  std::vector<LINALG::Matrix<ndim, numdof, double>> lin_theta_cp(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<ndim, numdof, double>(true));

  // re-interpolated lin_theta:
  LINALG::Matrix<ndim, numdof, double> lin_theta_bar(true);


  // set nodal / cp quantities: positions, tangents, triads
  // *******************************************************************************

  // Set current positions and orientations at all nodes:
  UpdateDispTotlag<nnodecl, double>(disp, disp_totlag);


  UpdateNodalVariables<nnodecl, double>(disp_totlag, disp_totlag_centerline, triad_mat_cp,
      Qref_dummy);  // Todo @grill: make this nicer and avoid dummies!


  // compute quantities at collocation points
  // *******************************************************************************
  for (unsigned int icp = 0; icp < BEAM3K_COLLOCATION_POINTS; ++icp)
  {
    // Determine storage position for this cp
    ind = LARGEROTATIONS::NumberingTrafo(icp + 1, BEAM3K_COLLOCATION_POINTS);

    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi_cp = (double)icp / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;

    // get value of interpolating function for theta (Lagrange polynomials) at xi_cp
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i, xi_cp, Shape());

    L.Clear();
    AssembleShapefunctionsL(L_i, L);

    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi_cp, length_, line2);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);

    // Calculation of r' at xi
    r_s.Clear();
    r_s.Multiply(N_s, disp_totlag_centerline);

    abs_r_s = FADUTILS::Norm(r_s);


    // lin_thetaperp
    Calc_lin_thetaperp<nnodecl>(lin_theta_perp_cp, N_s, r_s, abs_r_s);


    // v_lin_thetapar
    g_1_cp.Clear();
    g_1_cp.Update(std::pow(abs_r_s, -1.0), r_s);

    LINALG::Matrix<3, 3, double> triad_ref_conv_cp(true);
    LARGEROTATIONS::quaterniontotriad(Qrefconv_[ind], triad_ref_conv_cp);

    g_1_cp_bar.Clear();
    for (unsigned int idim = 0; idim < ndim; ++idim) g_1_cp_bar(idim) = triad_ref_conv_cp(idim, 0);

    Calc_lin_thetapar<nnodecl>(lin_theta_par_cp, L, N_s, g_1_cp, g_1_cp_bar, abs_r_s);

    // lin_theta
    lin_theta_cp[ind].Clear();
    lin_theta_cp[ind].Update(1.0, lin_theta_par_cp, 1.0, lin_theta_perp_cp);
  }


  // re-interpolation of quantities at xi based on CP values
  // *******************************************************************************
  std::vector<LINALG::Matrix<3, 3, double>> Itilde(
      BEAM3K_COLLOCATION_POINTS, LINALG::Matrix<3, 3, double>(true));

  // create object of triad interpolation scheme
  Teuchos::RCP<
      LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS, double>>
      triad_interpolation_scheme_ptr = Teuchos::rcp(
          new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<BEAM3K_COLLOCATION_POINTS,
              double>());

  // reset scheme with nodal triads
  triad_interpolation_scheme_ptr->Reset(triad_mat_cp);

  // compute Itilde matrices required for re-interpolation of CP values of lin_theta
  triad_interpolation_scheme_ptr->GetNodalGeneralizedRotationInterpolationMatricesAtXi(Itilde, xi);


  LINALG::Matrix<3, numdof, double> auxmatrix(true);

  lin_theta_bar.Clear();
  for (unsigned int inode = 0; inode < BEAM3K_COLLOCATION_POINTS; ++inode)
  {
    auxmatrix.Clear();

    auxmatrix.Multiply(Itilde[inode], lin_theta_cp[inode]);

    lin_theta_bar.Update(1.0, auxmatrix, 1.0);
  }


  // finally assemble the generalized interpolation matrix for the increments Iinc
  // *******************************************************************************

  // this part is associated with the increment of the centerline orientation
  // (expressed as rotation vector theta)
  // (rows 4-6 of Iinc)
  for (unsigned int irow = 0; irow < lin_theta_bar.M(); ++irow)
    for (unsigned int icol = 0; icol < lin_theta_bar.N(); ++icol)
      Iinc(ndim + irow, icol) += lin_theta_bar(irow, icol);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector(
    const std::vector<double>& dofvec, std::vector<double>& dofvec_centerline,
    bool add_reference_values) const
{
  if (dofvec.size() != 15)
    dserror("size mismatch: expected 15 values for element state vector and got %d", dofvec.size());

  dofvec_centerline.resize(12, 0.0);

  // we use the method for LINALG fixed size matrix and create it as a view on the STL vector
  LINALG::Matrix<15, 1, double> dofvec_fixedsize(&dofvec[0]);
  LINALG::Matrix<12, 1, double> dofvec_centerline_fixedsize(&dofvec_centerline[0], true);

  this->ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
      dofvec_fixedsize, dofvec_centerline_fixedsize, add_reference_values);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector(
    const LINALG::Matrix<3 * vpernode * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& dofvec,
    LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& dofvec_centerline,
    bool add_reference_values) const
{
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value
  // + derivative of value (i.e. Hermite))

  const int dofperboundarynode = 3 * vpernode + 1;

  // get current values for position DOFs directly
  for (unsigned int dim = 0; dim < 3; ++dim)
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_centerline(3 * vpernode * node + dim) = dofvec(dofperboundarynode * node + dim);
    }

  // Hermite interpolation: get values of tangent DOFs as well
  if (rotvec_ == false)
  {
    for (unsigned int dim = 0; dim < 3; ++dim)
      for (unsigned int node = 0; node < nnodecl; ++node)
      {
        // in case of rotvec_==false, tangents are equivalent to nodal DOFs
        dofvec_centerline(3 * vpernode * node + 3 + dim) =
            dofvec(dofperboundarynode * node + 3 + dim);
      }
  }
  else
  {
    // values for tangent DOFs must be transformed in case of rotvec_==true
    LINALG::Matrix<3, 1, T> theta(true);
    LINALG::Matrix<3, 3, T> unity(true);
    LINALG::Matrix<3, 3, T> triad(true);

    for (unsigned int dim = 0; dim < 3; ++dim) unity(dim, dim) = 1.0;

    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      for (unsigned int dim = 0; dim < 3; ++dim)
      {
        // get values for rotation vector DOFs
        theta(dim) = dofvec(dofperboundarynode * node + 3 + dim);
      }
      // transform to triad
      LARGEROTATIONS::angletotriad(theta, unity, triad);

      // direction of tangent is equivalent to first base vector of triad; length of tangent is 7th
      // nodal DOF
      for (unsigned int dim = 0; dim < 3; ++dim)
      {
        dofvec_centerline(3 * vpernode * node + 3 + dim) =
            triad(dim, 0) * dofvec(dofperboundarynode * node + 6);
      }
    }
  }

  if (add_reference_values) AddRefValuesDispCenterline<nnodecl, vpernode, T>(dofvec_centerline);
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline(
    LINALG::Matrix<3 * vpernode * nnodecl, 1, T>& dofvec_centerline) const
{
  for (unsigned int dim = 0; dim < 3; ++dim)
    for (unsigned int node = 0; node < nnodecl; ++node)
    {
      dofvec_centerline(3 * vpernode * node + dim) += Nodes()[node]->X()[dim];

      // Hermite interpolation: update tangent DOFs as well
      if (vpernode == 2) dofvec_centerline(3 * vpernode * node + 3 + dim) += T0_[node](dim);
    }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::AddRefValuesDisp(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& dofvec) const
{
  // For the boundary nodes we have to add the initial values. For the interior nodes we are already
  // done since the initial relative angles are zero: alpha_i(t=0)=0;
  if (rotvec_ == false)
  {
    // Calculate total displacements = positions vectors, tangents and relative angles
    for (unsigned int node = 0; node < 2; ++node)  // loop over boundary nodes
    {
      for (unsigned int ndof = 0; ndof < 7; ++ndof)  // loop over dofs per node
      {
        if (ndof < 3)
        {
          dofvec(7 * node + ndof) += (Nodes()[node])->X()[ndof];
        }
        else if (ndof < 6)
        {
          dofvec(7 * node + ndof) += (Tref()[node])(ndof - 3);
        }
        else
        {
          // nothing to do here: alpha_i(t=0)=0;
        }
      }
    }
  }
  else
  {
    // Calculate total displacements = positions vectors, tangents and relative angles
    for (unsigned int node = 0; node < 2; ++node)  // loop over boundary nodes
    {
      for (unsigned int ndof = 0; ndof < 7; ++ndof)  // loop over dofs per node
      {
        if (ndof < 3)
        {
          dofvec(7 * node + ndof) += (Nodes()[node])->X()[ndof];
        }
        else if (ndof < 6)
        {
          // Nothing to do here, rotations are treated below
        }
        else
        {
          // here we have to add the initial length of the tangents at the boundary nodes,
          // i.e. ||r'_i(t=0)||=1:
          dofvec(7 * node + ndof) += 1.0;
        }
      }  // for (int ndof=0;ndof<7;ndof++)//loop over dofs per node

      LINALG::Matrix<3, 1> disptheta(true);
      LINALG::Matrix<3, 1> thetanew(true);
      LINALG::Matrix<4, 1> deltaQ(true);
      LINALG::Matrix<4, 1> Qnew(true);
      for (unsigned int i = 0; i < 3; ++i)
      {
        disptheta(i) = FADUTILS::CastToDouble(dofvec(7 * node + 3 + i));
      }

      LARGEROTATIONS::angletoquaternion(disptheta, deltaQ);

      LINALG::Matrix<4, 1> Q0;
      LARGEROTATIONS::angletoquaternion(Theta0()[node], Q0);
      LARGEROTATIONS::quaternionproduct(Q0, deltaQ, Qnew);

      // renormalize quaternion to keep its absolute value one even in case of long simulations
      // and intricate calculations
      Qnew.Scale(1.0 / Qnew.Norm2());

      // Calculate the new nodal angle thetanew \in ]-PI,PI] -> Here, thetanew \in ]-PI,PI] by
      // quaterniontoangle()
      LARGEROTATIONS::quaterniontoangle(Qnew, thetanew);

      // Finally set rotation values in disp_totlag
      for (unsigned int i = 0; i < 3; ++i)
      {
        dofvec(7 * node + 3 + i) = thetanew(i);
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::UpdateDispTotlag(const std::vector<double>& disp,
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag) const
{
  disp_totlag.Clear();
  for (unsigned int dof = 0; dof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++dof)
    disp_totlag(dof) = disp[dof];

  AddRefValuesDisp<nnodecl, T>(disp_totlag);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::UpdateNodalVariables(
    const LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline,
    std::vector<LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<LINALG::Matrix<4, 1>>& Qref_new) const
{
  // Set positions vectors and tangents and triads at boundary nodes
  SetPositionsAtBoundaryNodes<nnodecl, T>(disp_totlag, disp_totlag_centerline);

  // next, set triads and tangents at boundary nodes
  SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<nnodecl, T>(
      disp_totlag, disp_totlag_centerline, triad_mat_cp, Qref_new);

  // finally, set triads at remaining CPs (all except boundary nodes)
  SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<nnodecl, T>(
      disp_totlag, disp_totlag_centerline, triad_mat_cp, Qref_new);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes(
    const LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline) const
{
  for (unsigned int i = 0; i < 3; ++i)
  {
    disp_totlag_centerline(i) = disp_totlag(i);
    disp_totlag_centerline(7 + i) = disp_totlag(7 + i);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes(
    const LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline,
    std::vector<LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<LINALG::Matrix<4, 1>>& Qref_new) const
{
  if (rotvec_ == false)
  {
    LINALG::Matrix<3, 1, T> tangent(true);
    LINALG::Matrix<3, 3, T> triad_ref(true);
    LINALG::Matrix<3, 3> triad_aux(true);
    // Todo @grill: get rid of auxiliary matrix of different type
    LINALG::Matrix<3, 3, T> triad_aux2(true);
    T alpha = 0.0;

    for (unsigned int node = 0; node < 2; ++node)
    {
      // tangent
      for (unsigned int i = 0; i < 3; ++i)
      {
        disp_totlag_centerline(7 * node + 3 + i) = disp_totlag(7 * node + 3 + i);
        tangent(i) = disp_totlag(7 * node + 3 + i);
      }

      alpha = disp_totlag(7 * node + 6);

      // calculate new sr triad
      triad_ref.Clear();
      triad_aux.Clear();
      triad_aux2.Clear();
      LARGEROTATIONS::quaterniontotriad(Qrefconv_[node], triad_aux);
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j) triad_aux2(i, j) = triad_aux(i, j);

      LARGEROTATIONS::CalculateSRTriads<T>(tangent, triad_aux2, triad_ref);

      // Store nodal reference triad
      LINALG::Matrix<4, 1, T> Qref(true);
      LARGEROTATIONS::triadtoquaternion(triad_ref, Qref);
      for (unsigned int i = 0; i < 4; ++i) Qref_new[node](i) = FADUTILS::CastToDouble(Qref(i));

      // calculate material triad
      triad_mat_cp[node].Clear();
      LARGEROTATIONS::RotateTriad(triad_ref, alpha, triad_mat_cp[node]);
    }
  }
  else
  {
    LINALG::Matrix<3, 1, T> theta(true);
    LINALG::Matrix<3, 3, T> unity(true);
    for (unsigned int node = 0; node < 2; ++node)
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        theta(i) = disp_totlag(7 * node + 3 + i);
        unity(i, i) = 1.0;
      }
      triad_mat_cp[node].Clear();
      LARGEROTATIONS::angletotriad(theta, unity, triad_mat_cp[node]);

      // tangent
      for (unsigned int i = 0; i < 3; ++i)
      {
        disp_totlag_centerline(7 * node + 3 + i) =
            (triad_mat_cp[node])(i, 0) * disp_totlag(7 * node + 6);
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints(
    const LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag,
    const LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, T>& disp_totlag_centerline,
    std::vector<LINALG::Matrix<3, 3, T>>& triad_mat_cp,
    std::vector<LINALG::Matrix<4, 1>>& Qref_new) const
{
  LINALG::Matrix<3, 1, T> tangent(true);
  LINALG::Matrix<3, 3, T> triad_ref(true);
  LINALG::Matrix<3, 3> triad_aux(true);
  // Todo @grill: get rid of auxiliary matrix of different type
  LINALG::Matrix<3, 3, T> triad_aux2(true);
  T alpha = 0.0;

  LINALG::Matrix<1, 4, T> N_i_xi(true);
  LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, T> N_s(true);
  LINALG::Matrix<1, BEAM3K_COLLOCATION_POINTS, T> L_i;
  double xi = 0.0;
  unsigned int ind = 0;

  //********begin: evaluate quantities at collocation points********************************
  // intermediate nodes: all but first (node=0) and last node (node=BEAM3K_COLLOCATION_POINTS)
  for (unsigned int node = 1; node < BEAM3K_COLLOCATION_POINTS - 1; ++node)
  {
    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0 node=2->xi=1
    xi = (double)node / (double)(BEAM3K_COLLOCATION_POINTS - 1) * 2.0 - 1.0;
    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi, xi, length_, line2);
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i, xi, Shape());

    // Determine storage position for the node node
    ind = LARGEROTATIONS::NumberingTrafo(node + 1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi, jacobi_cp_[ind], N_s);

    tangent.Clear();
    // Calculation of r' at xi
    tangent.Multiply(N_s, disp_totlag_centerline);

    alpha = disp_totlag(7 * 2 + ind - 2);

    // calculate new sr triads
    triad_ref.Clear();
    triad_aux.Clear();
    triad_aux2.Clear();
    LARGEROTATIONS::quaterniontotriad(Qrefconv_[ind], triad_aux);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j) triad_aux2(i, j) = triad_aux(i, j);

    LARGEROTATIONS::CalculateSRTriads<T>(tangent, triad_aux2, triad_ref);

    // Store nodal reference triad
    LINALG::Matrix<4, 1, T> Qref(true);
    LARGEROTATIONS::triadtoquaternion(triad_ref, Qref);
    for (unsigned int i = 0; i < 4; ++i) Qref_new[ind](i) = FADUTILS::CastToDouble(Qref(i));

    // calculate material triad
    triad_mat_cp[ind].Clear();
    LARGEROTATIONS::RotateTriad(triad_ref, alpha, triad_mat_cp[ind]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::SetAutomaticDifferentiationVariables(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 1, FAD>& disp_totlag) const
{
  for (unsigned int dof = 0; dof < nnodecl * 6 + BEAM3K_COLLOCATION_POINTS; dof++)
  {
    disp_totlag(dof).diff(dof, nnodecl * 6 + BEAM3K_COLLOCATION_POINTS);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3k::Calc_velocity(
    const LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
    const LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
    LINALG::Matrix<ndim, 1, double>& velocity, const LINALG::Matrix<ndim, 1, double>& position,
    int gausspoint_index) const
{
  DRT::UTILS::BEAM::CalcInterpolation<nnodecl, vpernode, ndim, double>(
      velocity_dofvec, N_i, velocity);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3k::Calc_velocity(
    const LINALG::Matrix<ndim * vpernode * nnodecl, 1, double>& velocity_dofvec,
    const LINALG::Matrix<1, vpernode * nnodecl, double>& N_i,
    LINALG::Matrix<ndim, 1, FAD>& velocity, const LINALG::Matrix<ndim, 1, FAD>& position,
    int gausspoint_index)
{
  /* if we use FAD, we need to track the dependency of velocity vector on primary variables, i.e.
   * we must calculate it like this using the class variable rconvmass_ for position vector at GP
   * of last converged state this again is tricky in case of periodic boundary conditions, because
   * the position of this element (nodes) might have been shifted outside; we therefore manually
   * adapt the calculated velocity and compare it with the velocity calculated in time integrator
   * and handed in from outside for safety reasons */
  LINALG::Matrix<ndim, 1, double> velocity_test;
  DRT::UTILS::BEAM::CalcInterpolation<nnodecl, vpernode, ndim, double>(
      velocity_dofvec, N_i, velocity_test);

  // get time step size
  const double dt = ParamsInterface().GetDeltaTime();

  LINALG::Matrix<3, 1> diff(true);

  LINALG::Matrix<ndim, 1, FAD> delta_r_ost(true);
  Teuchos::RCP<GEO::MESHFREE::BoundingBox> pbb =
      BrownianDynParamsInterface().GetPeriodicBoundingBox();

  LINALG::Matrix<3, 1> unshiftedrconvmass_i(true), position_i_double(true);

  for (unsigned int idim = 0; idim < ndim; ++idim)
  {
    unshiftedrconvmass_i(idim) = rconvmass_[gausspoint_index](idim);
    position_i_double(idim) = FADUTILS::CastToDouble(position(idim));
  }

  // difference in position of this GP as compared to last time step
  pbb->UnShift3D(unshiftedrconvmass_i, position_i_double);

  for (unsigned int idim = 0; idim < ndim; ++idim)
  {
    // difference in position of this GP as compared to last time step
    delta_r_ost(idim) = position(idim) - unshiftedrconvmass_i(idim);

    // velocity according to Backward Euler scheme
    velocity(idim) = delta_r_ost(idim) / dt;

    diff(idim) = velocity(idim).val() - velocity_test(idim);

    // set class variable, such that rconvmass_ is available in next time step
    rnewmass_[gausspoint_index](idim) = position(idim).val();
  }

  // safety check
  if (diff.NormInf() > 1e-11)
  {
    std::cout << "\nrnewmass = " << position;
    std::cout << "\nrconvmass_ = " << rconvmass_[gausspoint_index];
    std::cout << "\nvel = " << velocity;
    std::cout << "\nvel_test = " << velocity_test;
    std::cout << "\nabs(diff) = " << diff.Norm2();
    std::cout << "\n*** SERIOUS WARNING: ***";
    std::cout << "\nvelocity vector at GP computed locally in beam3k differs from OneStepTheta "
                 "velocity vector (see values above)!\n\n";
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetaperp,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& N_s,
    const LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const
{
  v_thetaperp.Clear();

  LINALG::Matrix<3, 3, T> S_of_r_s(true);
  LARGEROTATIONS::computespin<T>(S_of_r_s, r_s);

  v_thetaperp.MultiplyTN(N_s, S_of_r_s);
  v_thetaperp.Scale(-1.0 * std::pow(abs_r_s, -2.0));
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl, typename T>
void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, T>& v_thetapartheta,
    const LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, T>& L,
    const LINALG::Matrix<3, 1, T>& r_s, T abs_r_s) const
{
  v_thetapartheta.Clear();

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      v_thetapartheta(idof, jdim) = L(idof) * r_s(jdim) / abs_r_s;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_thetaperp(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetaperp,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& r_s, double abs_r_s) const
{
  // Todo @grill: maybe re-use method for v_thetaperp here, which is simply the transpose of this
  // term

  lin_thetaperp.Clear();

  LINALG::Matrix<3, 3, double> S_of_r_s(true);
  LARGEROTATIONS::computespin<double>(S_of_r_s, r_s);

  lin_thetaperp.Multiply(S_of_r_s, N_s);
  lin_thetaperp.Scale(std::pow(abs_r_s, -2.0));
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_thetapar(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_thetapar,
    const LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_bar,
    double abs_r_s) const
{
  const unsigned int numdof = 6 * nnodecl + BEAM3K_COLLOCATION_POINTS;

  lin_thetapar.Clear();

  LINALG::Matrix<3, 3, double> S_of_g_1(true);
  LARGEROTATIONS::computespin<double>(S_of_g_1, g_1);


  // Todo @grill: decide about alternatives, apparently no change in results

  // *********** alternative 1 ************
  LINALG::Matrix<3, 3, double> g_1_dyadicproduct_g_1_bar_T(true);

  for (unsigned int irow = 0; irow < 3; ++irow)
    for (unsigned int icol = 0; icol < 3; ++icol)
      g_1_dyadicproduct_g_1_bar_T(irow, icol) = g_1(irow) * g_1_bar(icol);

  LINALG::Matrix<3, 3, double> g_1_dyadicproduct_g_1_bar_T_S_of_g_1(true);
  g_1_dyadicproduct_g_1_bar_T_S_of_g_1.Multiply(g_1_dyadicproduct_g_1_bar_T, S_of_g_1);

  lin_thetapar.Multiply(g_1_dyadicproduct_g_1_bar_T_S_of_g_1, N_s);

  // *********** alternative 2 ************

  //  LINALG::Matrix<1,3,double> g_1_bar_T_S_of_g_1(true);
  //  g_1_bar_T_S_of_g_1.MultiplyTN(g_1_bar, S_of_g_1);
  //
  //
  //  LINALG::Matrix<1,numdof,double> g_1_bar_T_S_of_g_1_N_s(true);
  //  g_1_bar_T_S_of_g_1_N_s.Multiply(g_1_bar_T_S_of_g_1, N_s);
  //
  //
  //  // dyadic product
  //  for (unsigned int idim=0; idim<3; ++idim)
  //    for (unsigned int jdof=0; jdof<numdof; ++jdof)
  //      lin_thetapar(idim,jdof) += g_1(idim) * g_1_bar_T_S_of_g_1_N_s(jdof) ;

  // *************************************

  lin_thetapar.Scale(-1.0 / (1.0 + g_1.Dot(g_1_bar)));
  lin_thetapar.Scale(1.0 / abs_r_s);


  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdof = 0; jdof < numdof; ++jdof)
      lin_thetapar(idim, jdof) += g_1(idim) * L(jdof);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_tangent_tilde,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_tangent_tilde.Clear();

  LINALG::Matrix<3, 3, double> auxmatrix;

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -2.0));

  lin_tangent_tilde.Multiply(auxmatrix, N_s);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde_s(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_tangent_tilde_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_s,
    const LINALG::Matrix<3, 1, double>& r_s, const LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s) const
{
  lin_tangent_tilde_s.Clear();

  LINALG::Matrix<3, 3, double> auxmatrix;
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> summand(true);

  // first summand
  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(2.0 * r_s.Dot(r_ss) * std::pow(abs_r_s, -4.0));

  lin_tangent_tilde_s.Multiply(auxmatrix, N_s);

  // second summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = g_1_s(idim) * g_1(jdim) + g_1(idim) * g_1_s(jdim);

  auxmatrix.Scale(-2.0 * std::pow(abs_r_s, -2.0));

  summand.Multiply(auxmatrix, N_s);

  lin_tangent_tilde_s.Update(1.0, summand, 1.0);

  // third summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - 2.0 * g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -2.0));

  summand.Clear();

  summand.Multiply(auxmatrix, N_ss);

  lin_tangent_tilde_s.Update(1.0, summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_g_1(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_g_1.Clear();

  LINALG::Matrix<3, 3, double> auxmatrix;

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -1.0));

  lin_g_1.Multiply(auxmatrix, N_s);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_g_1_s(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_g_1_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_s,
    const LINALG::Matrix<3, 1, double>& r_s, const LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s) const
{
  lin_g_1_s.Clear();

  LINALG::Matrix<3, 3, double> auxmatrix;
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> summand(true);

  // first summand
  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(-1.0 * r_s.Dot(r_ss) * std::pow(abs_r_s, -3.0));

  lin_g_1_s.Multiply(auxmatrix, N_s);

  // second summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = g_1_s(idim) * g_1(jdim) + g_1(idim) * g_1_s(jdim);

  auxmatrix.Scale(-1.0 * std::pow(abs_r_s, -1.0));

  summand.Multiply(auxmatrix, N_s);

  lin_g_1_s.Update(1.0, summand, 1.0);

  // third summand
  auxmatrix.Clear();

  for (unsigned int idim = 0; idim < 3; ++idim)
    for (unsigned int jdim = 0; jdim < 3; ++jdim)
      auxmatrix(idim, jdim) = 1.0 * (idim == jdim) - g_1(idim) * g_1(jdim);

  auxmatrix.Scale(std::pow(abs_r_s, -1.0));

  summand.Clear();

  summand.Multiply(auxmatrix, N_ss);

  lin_g_1_s.Update(1.0, summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_epsilon(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_epsilon,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s) const
{
  lin_v_epsilon.Clear();

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  lin_v_epsilon.MultiplyTN(N_s, lin_g_1);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_resultant(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_resultant,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta_s,
    const LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
    const LINALG::Matrix<3, 3, double>& cm) const
{
  lin_moment_resultant.Clear();

  // first summand
  lin_moment_resultant.Multiply(spinmatrix_of_moment, lin_theta);
  lin_moment_resultant.Scale(-1.0);

  // second summand
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix(true);

  auxmatrix.Multiply(cm, lin_theta_s);

  lin_moment_resultant.Update(1.0, auxmatrix, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_inertia(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_inertia,
    const LINALG::Matrix<3, 3, double>& triad_mat,
    const LINALG::Matrix<3, 3, double>& triad_mat_conv,
    const LINALG::Matrix<3, 1, double>& deltatheta,
    const LINALG::Matrix<3, 1, double>& angular_velocity_material,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const LINALG::Matrix<3, 3, double>& spinmatrix_of_moment,
    const LINALG::Matrix<3, 3, double>& C_rho, double lin_prefactor_acc,
    double lin_prefactor_vel) const
{
  lin_moment_inertia.Clear();

  // first summand
  lin_moment_inertia.Multiply(spinmatrix_of_moment, lin_theta);

  // second summand
  LINALG::Matrix<3, 3, double> auxmatrix(true);

  LINALG::Matrix<3, 3, double> spinmatrix(true);
  LARGEROTATIONS::computespin(spinmatrix, angular_velocity_material);

  auxmatrix.Multiply(spinmatrix, C_rho);

  LINALG::Matrix<3, 1, double> C_rho_W(true);
  C_rho_W.Multiply(C_rho, angular_velocity_material);

  spinmatrix.Clear();
  LARGEROTATIONS::computespin(spinmatrix, C_rho_W);
  auxmatrix.Update(-1.0, spinmatrix, 1.0);
  auxmatrix.Scale(lin_prefactor_vel);

  auxmatrix.Update(lin_prefactor_acc, C_rho, 1.0);


  LINALG::Matrix<3, 3, double> Tmat_of_deltatheta = LARGEROTATIONS::Tmatrix(deltatheta);
  LINALG::Matrix<3, 3, double> Lambda_conv_Tmat_of_deltatheta(true);
  Lambda_conv_Tmat_of_deltatheta.MultiplyTN(triad_mat_conv, Tmat_of_deltatheta);

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_THETA_tilde(true);
  lin_THETA_tilde.Multiply(Lambda_conv_Tmat_of_deltatheta, lin_theta);

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix2(true);
  auxmatrix2.Multiply(auxmatrix, lin_THETA_tilde);


  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> auxmatrix3(true);
  auxmatrix3.Multiply(triad_mat, auxmatrix2);

  lin_moment_inertia.Update(1.0, auxmatrix3, 1.0);

  lin_moment_inertia.Scale(-1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_moment_viscous(
    LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_moment_viscous,
    const LINALG::Matrix<3, 3, double>& triad_mat,
    const LINALG::Matrix<3, 3, double>& triad_mat_conv,
    const LINALG::Matrix<3, 1, double>& deltatheta,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& lin_theta,
    const LINALG::Matrix<3, 3, double>& spinmatrix_of_moment, double gamma_polar, double dt) const
{
  lin_moment_viscous.Clear();

  // first summand
  lin_moment_viscous.Multiply(spinmatrix_of_moment, lin_theta);
  lin_moment_viscous.Scale(-1.0);

  // second summand
  LINALG::Matrix<3, 3, double> auxmatrix(true);

  LINALG::Matrix<3, 3, double> gamma_g1_g1_conv;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      gamma_g1_g1_conv(i, j) = triad_mat(i, 0) * triad_mat_conv(j, 0) * gamma_polar / dt;

  LINALG::Matrix<3, 3, double> Tmat_of_deltatheta = LARGEROTATIONS::Tmatrix(deltatheta);

  auxmatrix.Multiply(gamma_g1_g1_conv, Tmat_of_deltatheta);

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> second_summand(true);
  second_summand.Multiply(auxmatrix, lin_theta);

  lin_moment_viscous.Update(1.0, second_summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_moment(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetaperp_moment,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const
{
  lin_v_thetaperp_moment.Clear();

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde(true);

  Calc_lin_tangent_tilde<nnodecl>(lin_tangent_tilde, N_s, g_1, abs_r_s);


  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  // Todo @grill: check: is the order of matrix products relevant?
  auxmatrix.MultiplyTN(N_s, spinmatrix_of_moment);

  lin_v_thetaperp_moment.Multiply(auxmatrix, lin_tangent_tilde);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_s_moment(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetaperp_s_moment,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_s,
    const LINALG::Matrix<3, 1, double>& r_s, const LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const LINALG::Matrix<3, 3, double>& spinmatrix_of_moment) const
{
  lin_v_thetaperp_s_moment.Clear();

  // first summand
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde(true);

  Calc_lin_tangent_tilde<nnodecl>(lin_tangent_tilde, N_s, g_1, abs_r_s);

  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  auxmatrix.MultiplyTN(N_ss, spinmatrix_of_moment);

  lin_v_thetaperp_s_moment.Multiply(auxmatrix, lin_tangent_tilde);

  // second summand
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_tangent_tilde_s(true);

  Calc_lin_tangent_tilde_s<nnodecl>(lin_tangent_tilde_s, N_s, N_ss, g_1, g_1_s, r_s, r_ss, abs_r_s);

  auxmatrix.Clear();
  auxmatrix.MultiplyTN(N_s, spinmatrix_of_moment);

  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
      double>
      second_summand(true);

  second_summand.Multiply(auxmatrix, lin_tangent_tilde_s);

  lin_v_thetaperp_s_moment.Update(1.0, second_summand, 1.0);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_moment(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetapar_moment,
    LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const LINALG::Matrix<3, 1, double>& moment) const
{
  lin_v_thetapar_moment.Clear();

  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  // Todo @grill: check: is the order of matrix products relevant?
  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L(idof) * moment(jdim);

  lin_v_thetapar_moment.Multiply(auxmatrix, lin_g_1);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int nnodecl>
void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_s_moment(
    LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
        double>& lin_v_thetapar_s_moment,
    LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L,
    LINALG::Matrix<1, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& L_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_s,
    const LINALG::Matrix<3, 1, double>& r_s, const LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const LINALG::Matrix<3, 1, double>& moment) const
{
  lin_v_thetapar_s_moment.Clear();

  // first summand
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1(true);

  Calc_lin_g_1<nnodecl>(lin_g_1, N_s, g_1, abs_r_s);

  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 3, double> auxmatrix(true);

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L_s(idof) * moment(jdim);

  lin_v_thetapar_s_moment.Multiply(auxmatrix, lin_g_1);

  // second summand
  LINALG::Matrix<3, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS, double> lin_g_1_s(true);

  Calc_lin_g_1_s<nnodecl>(lin_g_1_s, N_s, N_ss, g_1, g_1_s, r_s, r_ss, abs_r_s);

  auxmatrix.Clear();

  for (unsigned int idof = 0; idof < 6 * nnodecl + BEAM3K_COLLOCATION_POINTS; ++idof)
    for (unsigned int jdim = 0; jdim < 3; ++jdim) auxmatrix(idof, jdim) = L(idof) * moment(jdim);

  LINALG::Matrix<6 * nnodecl + BEAM3K_COLLOCATION_POINTS, 6 * nnodecl + BEAM3K_COLLOCATION_POINTS,
      double>
      second_summand(true);

  second_summand.Multiply(auxmatrix, lin_g_1_s);

  lin_v_thetapar_s_moment.Update(1.0, second_summand, 1.0);
}

// explicit template instantiations
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector<2, 2,
    Sacado::Fad::DFad<double>>(
    const LINALG::Matrix<12 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&, bool) const;
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValuesFromElementStateVector<2, 2, double>(
    const LINALG::Matrix<12 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    LINALG::Matrix<12, 1, double>&, bool) const;

template void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline<2, 2, Sacado::Fad::DFad<double>>(
    LINALG::Matrix<12, 1, Sacado::Fad::DFad<double>>&) const;
template void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline<2, 2, double>(
    LINALG::Matrix<12, 1, double>&) const;

template void DRT::ELEMENTS::Beam3k::AddRefValuesDisp<2, double>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::AddRefValuesDisp<2, Sacado::Fad::DFad<double>>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::UpdateDispTotlag<2, double>(const std::vector<double>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::UpdateDispTotlag<2, Sacado::Fad::DFad<double>>(
    const std::vector<double>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::UpdateNodalVariables<2, double>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<LINALG::Matrix<3, 3, double>>&, std::vector<LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::UpdateNodalVariables<2, Sacado::Fad::DFad<double>>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes<2, double>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&) const;
template void DRT::ELEMENTS::Beam3k::SetPositionsAtBoundaryNodes<2, Sacado::Fad::DFad<double>>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<2,
    double>(const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<LINALG::Matrix<3, 3, double>>&, std::vector<LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::SetTangentsAndTriadsAndReferenceTriadsAtBoundaryNodes<2,
    Sacado::Fad::DFad<double>>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<2,
    double>(const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, double>&,
    std::vector<LINALG::Matrix<3, 3, double>>&, std::vector<LINALG::Matrix<4, 1>>&) const;
template void DRT::ELEMENTS::Beam3k::SetTriadsAndReferenceTriadsAtRemainingCollocationPoints<2,
    Sacado::Fad::DFad<double>>(
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    const LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, Sacado::Fad::DFad<double>>&,
    std::vector<LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>>&,
    std::vector<LINALG::Matrix<4, 1>>&) const;

template void DRT::ELEMENTS::Beam3k::SetAutomaticDifferentiationVariables<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 1, FAD>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_velocity<2, 2, 3>(
    const LINALG::Matrix<3 * 2 * 2, 1, double>&, const LINALG::Matrix<1, 2 * 2, double>&,
    LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&, int) const;

template void DRT::ELEMENTS::Beam3k::Calc_velocity<2, 2, 3>(
    const LINALG::Matrix<3 * 2 * 2, 1, double>&, const LINALG::Matrix<1, 2 * 2, double>&,
    LINALG::Matrix<3, 1, FAD>&, const LINALG::Matrix<3, 1, FAD>&, int);

template void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp<2, double>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;
template void DRT::ELEMENTS::Beam3k::Calc_v_thetaperp<2, Sacado::Fad::DFad<double>>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, Sacado::Fad::DFad<double>>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, Sacado::Fad::DFad<double>>&,
    const LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&, Sacado::Fad::DFad<double>) const;

template void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta<2, double>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, double>&,
    const LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;
template void DRT::ELEMENTS::Beam3k::Calc_v_thetapartheta<2, Sacado::Fad::DFad<double>>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 3, Sacado::Fad::DFad<double>>&,
    const LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, Sacado::Fad::DFad<double>>&,
    const LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&, Sacado::Fad::DFad<double>) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_thetaperp<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_thetapar<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_tangent_tilde_s<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_g_1<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_g_1_s<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_epsilon<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_resultant<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 3, double>&, const LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_inertia<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 3, double>&, const LINALG::Matrix<3, 3, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 3, double>&, const LINALG::Matrix<3, 3, double>&, double, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_moment_viscous<2>(
    LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 3, double>&, const LINALG::Matrix<3, 3, double>&,
    const LINALG::Matrix<3, 1, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 3, double>&, double, double) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_moment<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, double, const LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_moment<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&
        lin_v_thetapar_moment,
    LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 1, double>& g_1, double abs_r_s,
    const LINALG::Matrix<3, 1, double>& moment) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetaperp_s_moment<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&,
    const LINALG::Matrix<3, 1, double>&, const LINALG::Matrix<3, 1, double>&, double,
    const LINALG::Matrix<3, 3, double>&) const;

template void DRT::ELEMENTS::Beam3k::Calc_lin_v_thetapar_s_moment<2>(
    LINALG::Matrix<6 * 2 + BEAM3K_COLLOCATION_POINTS, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>&
        lin_v_thetapar_s_moment,
    LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L,
    LINALG::Matrix<1, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& L_s,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_s,
    const LINALG::Matrix<3, 6 * 2 + BEAM3K_COLLOCATION_POINTS, double>& N_ss,
    const LINALG::Matrix<3, 1, double>& g_1, const LINALG::Matrix<3, 1, double>& g_1_s,
    const LINALG::Matrix<3, 1, double>& r_s, const LINALG::Matrix<3, 1, double>& r_ss,
    double abs_r_s, const LINALG::Matrix<3, 1, double>& moment) const;
