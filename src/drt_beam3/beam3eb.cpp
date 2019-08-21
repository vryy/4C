/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear torsionless rod based on a C1 curve

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3eb.H"

// Todo check for obsolete header inclusions
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_browniandyn.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"


DRT::ELEMENTS::Beam3ebType DRT::ELEMENTS::Beam3ebType::instance_;

DRT::ELEMENTS::Beam3ebType& DRT::ELEMENTS::Beam3ebType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Beam3ebType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Beam3eb* object = new DRT::ELEMENTS::Beam3eb(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "BEAM3EB")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3eb(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebType::Create(const int id, const int owner)

{
  return Teuchos::rcp(new Beam3eb(id, owner));
}


void DRT::ELEMENTS::Beam3ebType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 6;
  nv = 6;
  dimns = 3;
}

// TODO: the function ComputeNullSpace has still to be implemented
void DRT::ELEMENTS::Beam3ebType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  dserror("Function not implemented yet.");
  // DRT::UTILS::ComputeXFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam3ebType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["BEAM3EB"];

  defs["LIN2"].AddIntVector("LIN2", 2).AddNamedInt("MAT");

  defs["LIN3"].AddIntVector("LIN3", 3).AddNamedInt("MAT");
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebType::Initialize(DRT::Discretization& dis)
{
  // setting up geometric variables for beam3eb elements
  for (int num = 0; num < dis.NumMyColElements(); ++num)
  {
    // in case that current element is not a Truss3CL element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    // if we get so far current element is a beam3eb element and  we get a pointer at it
    DRT::ELEMENTS::Beam3eb* currele = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3eb* failed");

    // reference node position
    std::vector<double> xrefe;

    const int nnode = currele->NumNode();

    // resize xrefe for the number of coordinates we need to store
    xrefe.resize(3 * nnode);

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
          xrefe[node * 3 + dof] =
              currele->Nodes()[node]->X()[dof] + disp_shift[node * numdof + dof];
        }
    }

    currele->SetUpReferenceGeometry(xrefe);

  }  // for (int num=0; num<dis_.NumMyColElements(); ++num)
  return 0;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3eb::Beam3eb(int id, int owner)
    : DRT::ELEMENTS::Beam3Base(id, owner),
      isinit_(false),
      jacobi_(0.0),
      firstcall_(true),
      Ekin_(0.0),
      Eint_(0.0),
      L_(LINALG::Matrix<3, 1>(true)),
      P_(LINALG::Matrix<3, 1>(true)),
      t0_(LINALG::Matrix<3, 2>(true)),
      t_(LINALG::Matrix<3, 2>(true)),
      kappa_max_(0.0),
      epsilon_max_(0.0),
      axial_strain_GP_(0),
      curvature_GP_(0),
      axial_force_GP_(0),
      bending_moment_GP_(0)
{
#if defined(INEXTENSIBLE)
  if (ANSVALUES != 3 or NODALDOFS != 2)
    dserror("Flag INEXTENSIBLE only possible in combination with ANSVALUES=3 and NODALDOFS=2!");
#endif

  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3eb::Beam3eb(const DRT::ELEMENTS::Beam3eb& old)
    : DRT::ELEMENTS::Beam3Base(old),
      isinit_(old.isinit_),
      jacobi_(old.jacobi_),
      Ekin_(old.Ekin_),
      Eint_(old.Eint_),
      Tref_(old.Tref_),
      L_(old.L_),
      P_(old.P_),
      t0_(old.t0_),
      t_(old.t_),
      kappa_max_(old.kappa_max_),
      epsilon_max_(old.epsilon_max_),
      axial_strain_GP_(old.axial_strain_GP_),
      curvature_GP_(old.curvature_GP_),
      axial_force_GP_(old.axial_force_GP_),
      bending_moment_GP_(old.bending_moment_GP_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3eb and return pointer to it (public) |
 |                                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3eb::Clone() const
{
  DRT::ELEMENTS::Beam3eb* newelement = new DRT::ELEMENTS::Beam3eb(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3eb::~Beam3eb() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::Print(std::ostream& os) const
{
  os << "beam3eb ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3eb::Shape() const { return line2; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Beam3Base::Pack(data);

  // add all class variables
  AddtoPack(data, jacobi_);
  AddtoPack(data, isinit_);
  AddtoPack(data, Ekin_);
  AddtoPack(data, Eint_);
  AddtoPack(data, Tref_);
  AddtoPack<3, 1>(data, L_);
  AddtoPack<3, 1>(data, P_);
  AddtoPack<3, 2>(data, t0_);
  AddtoPack<3, 2>(data, t_);
  AddtoPack(data, kappa_max_);
  AddtoPack(data, epsilon_max_);
  AddtoPack(data, axial_strain_GP_);
  AddtoPack(data, curvature_GP_);
  AddtoPack(data, axial_force_GP_);
  AddtoPack(data, bending_moment_GP_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::Unpack(const std::vector<char>& data)
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

  // extract all class variables of beam3 element
  ExtractfromPack(position, data, jacobi_);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack(position, data, Ekin_);
  ExtractfromPack(position, data, Eint_);
  ExtractfromPack(position, data, Tref_);
  ExtractfromPack<3, 1>(position, data, L_);
  ExtractfromPack<3, 1>(position, data, P_);
  ExtractfromPack<3, 2>(position, data, t0_);
  ExtractfromPack<3, 2>(position, data, t_);
  ExtractfromPack(position, data, kappa_max_);
  ExtractfromPack(position, data, epsilon_max_);
  ExtractfromPack(position, data, axial_strain_GP_);
  ExtractfromPack(position, data, curvature_GP_);
  ExtractfromPack(position, data, axial_force_GP_);
  ExtractfromPack(position, data, bending_moment_GP_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Beam3eb::Lines()
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
 | element nodes are already known (public)                   meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::SetUpReferenceGeometry(
    const std::vector<double>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe and rotrefe; note: the
   *isinit_ flag is important for avoiding reinitialization upon restart. However, it should be
   *possible to conduct a second initilization in principle (e.g. for periodic boundary conditions*/

  const int nnode = 2;

  // low precission

  {
    if (!isinit_ || secondinit)
    {
      isinit_ = true;

      // Get DiscretizationType
      DRT::Element::DiscretizationType distype = Shape();

      // Get integrationpoints for exact integration
      DRT::UTILS::IntegrationPoints1D gausspoints =
          DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleeb);

      Tref_.resize(gausspoints.nquad);

      // assure correct size of strain and stress resultant class variables and fill them
      // with zeros (by definition, the reference configuration is undeformed and stress-free)
      axial_strain_GP_.resize(gausspoints.nquad);
      std::fill(axial_strain_GP_.begin(), axial_strain_GP_.end(), 0.0);

      curvature_GP_.resize(gausspoints.nquad);
      std::fill(curvature_GP_.begin(), curvature_GP_.end(), 0.0);

      axial_force_GP_.resize(gausspoints.nquad);
      std::fill(axial_force_GP_.begin(), axial_force_GP_.end(), 0.0);

      bending_moment_GP_.resize(gausspoints.nquad);
      std::fill(bending_moment_GP_.begin(), bending_moment_GP_.end(), 0.0);


      // create Matrix for the derivates of the shapefunctions at the GP
      LINALG::Matrix<1, nnode> shapefuncderiv;

      // Loop through all GPs and compute jacobi at the GPs
      for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
      {
        // Get position xi of GP
        const double xi = gausspoints.qxg[numgp][0];

        // Get derivatives of shapefunctions at GP --> for simplicity here are Lagrange polynomials
        // instead of Hermite polynomials used to calculate the reference geometry. Since the
        // reference geometry for this beam element must always be a straight line there is no
        // difference between theses to types of interpolation functions.
        DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv, xi, distype);

        Tref_[numgp].Clear();

        // calculate vector dxdxi
        for (int node = 0; node < nnode; node++)
        {
          for (int dof = 0; dof < 3; dof++)
          {
            Tref_[numgp](dof) += shapefuncderiv(node) * xrefe[3 * node + dof];
          }  // for(int dof=0; dof<3 ; dof++)
        }    // for(int node=0; node<nnode; node++)

        // Store length factor for every GP
        // note: the length factor jacobi replaces the determinant and refers to the reference
        // configuration by definition
        jacobi_ = Tref_[numgp].Norm2();

        Tref_[numgp].Scale(1 / jacobi_);
      }

      // compute tangent at each node
      double norm2 = 0.0;

      Tref_.resize(nnode);
#if NODALDOFS == 3
      Kref_.resize(gausspoints.nquad);
#endif

      for (int node = 0; node < nnode; node++)
      {
        Tref_[node].Clear();
#if NODALDOFS == 3
        Kref_[node].Clear();
#endif
        for (int dof = 0; dof < 3; dof++)
        {
          Tref_[node](dof) = xrefe[3 + dof] - xrefe[dof];
        }
        norm2 = Tref_[node].Norm2();
        Tref_[node].Scale(1 / norm2);

        for (int i = 0; i < 3; i++) t0_(i, node) = Tref_[node](i);
      }

    }  // if(!isinit_)
  }
  // end low precission
  /*//begin high precission
  {
    if(!isinit_ || secondinit)
    {
      isinit_ = true;

        for (int i=0; i<3;i++)
        {
          Trefprec_(i)="0.0e+0_40";
        }
          Trefprec_(0) = "1.0e+0_40";
          jacobiprec_ = "10.0_40";
          Izzprec_ = "0.785398_40";
          crosssecprec_ = "3.14159_40";

  //Xrefprec_.resize(6);


  for (int i=0;i<6;i++)
  {

    std::ostringstream strs;
    strs << std::scientific << std::setw( 12 ) << xrefe[i];
    std::string str = strs.str();
    //std::string str2 = str.substr( 0, 8 )+ "_40";
    std::string str2 = str + "_40";
    //cout << "xrefe" << i << ": "<< xrefe[i] << endl;
    //cout << "str2: " << str2 << endl;
    if(xrefe[i]==0.0)
      Xrefprec_[i]="0.0_40";
    else
    Xrefprec_[i]= str2.c_str();
    //cout << "Xrefprec_: " << Xrefprec_[i] << endl;

  }


        jacobiprec_ = cl_float(0.0,float_format(40));
        for (int i=0;i<3;i++)
        {
          jacobiprec_+= Trefprec_(i)*Trefprec_(i);
        }
        jacobiprec_= sqrt(jacobiprec_);

        Trefprec_.Scale(cl_float(1.0,float_format(40))/jacobiprec_);



    }//if(!isinit_)
  }//end high precission*/
  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::vector<LINALG::Matrix<3, 1>> DRT::ELEMENTS::Beam3eb::Tref() const { return Tref_; }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3eb::jacobi() const { return jacobi_; }

///*-----------------------------------------------------------------------------------------------*
// *-----------------------------------------------------------------------------------------------*/
// double DRT::ELEMENTS::Beam3eb::MomentOfInertiaY() const
//{
//  return Iyy_;
//}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::GetPosAtXi(
    LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const
{
  if (disp.size() != 12)
    dserror(
        "size mismatch: expected 12 values for element displacement vector "
        "and got %d",
        disp.size());

  // add reference positions and tangents => total Lagrangean state vector
  LINALG::Matrix<12, 1> disp_totlag(true);
  UpdateDispTotlag<2, 6>(disp, disp_totlag);

  // Todo @grill: rename this to GetPosAtXiFromDispTotlag and make it private; replace call from
  // outside
  //      by a call to this more general method (GetPosAtXi) based on displacement instead of
  //      absolute values
  pos = this->GetPos(xi, disp_totlag);

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::GetTriadAtXi(
    LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
{
  if (disp.size() != 12)
    dserror(
        "size mismatch: expected 12 values for element displacement vector "
        "and got %d",
        disp.size());

  // add reference positions and tangents => total Lagrangean state vector
  LINALG::Matrix<12, 1> disp_totlag(true);
  UpdateDispTotlag<2, 6>(disp, disp_totlag);

  triad.Clear();

  /* note: this beam formulation (Beam3eb = torsion-free, isotropic Kirchhoff beam)
   *       does not need to track material triads and therefore can not provide it here;
   *       instead, we return the unit tangent vector as first base vector; both are
   *       identical in the case of Kirchhoff beams (shear-free);
   *
   * Todo @grill: what to do with second and third base vector?
   *
   */

  dserror(
      "\nBeam3eb::GetTriadAtXi(): by definition, this element can not return "
      "a full triad; think about replacing it by GetTangentAtXi or another solution.");
}
