/*-----------------------------------------------------------------------*/
/*!
\file beam3k.cpp

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

\maintainer Maximilian Grill

*/
/*-----------------------------------------------------------------------*/

#include "beam3k.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include <Teuchos_TimeMonitor.hpp>

DRT::ELEMENTS::Beam3kType DRT::ELEMENTS::Beam3kType::instance_;

DRT::ELEMENTS::Beam3kType& DRT::ELEMENTS::Beam3kType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Beam3kType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3k* object = new DRT::ELEMENTS::Beam3k(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create(const std::string eletype,
                                   const std::string eledistype,
                                 const int id,
                                 const int owner )
{
  if ( eletype=="BEAM3K" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3k(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3kType::Create( const int id, const int owner )

{
  return Teuchos::rcp( new Beam3k( id, owner ) );
}

void DRT::ELEMENTS::Beam3kType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
      dserror("method 'NodalBlockInformation' not implemented for element type beam3k!");
}

void DRT::ELEMENTS::Beam3kType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  dserror("Function not implemented yet.");
}

void DRT::ELEMENTS::Beam3kType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3K"];
  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("WK")
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddOptionalNamedDouble("IT")
    .AddOptionalNamedDouble("IR1")
    .AddOptionalNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("WK")
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddOptionalNamedDouble("IT")
    .AddOptionalNamedDouble("IR1")
    .AddOptionalNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("WK")
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddOptionalNamedDouble("IT")
    .AddOptionalNamedDouble("IR1")
    .AddOptionalNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("WK")
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddOptionalNamedDouble("IT")
    .AddOptionalNamedDouble("IR1")
    .AddOptionalNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LINE4"]
      .AddIntVector("LINE4",4)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddOptionalNamedDouble("IT")
      .AddOptionalNamedDouble("IR1")
      .AddOptionalNamedDouble("IR2")
      .AddNamedDoubleVector("TRIADS",12)
      ;

    defs["LIN4"]
      .AddIntVector("LIN4",4)
      .AddNamedInt("WK")
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddOptionalNamedDouble("IT")
      .AddOptionalNamedDouble("IR1")
      .AddOptionalNamedDouble("IR2")
      .AddNamedDoubleVector("TRIADS",12)
      ;

}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(int id, int owner) :
 DRT::ELEMENTS::Beam3Base(id,owner),
isinit_(false),
statmechprob_(false),
crosssec_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
T0_(0),
T_(0),
theta0_(0),
dispthetanew_(0),
Qnew_(0),
Qrefconv_(0),
Qrefnew_(0),
K0_(0),
length_(0),
jacobi_(0),
jacobi2_(0),
jacobi_cp_(0),
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
inertscaletrans_(0),
inertscalerot1_(0),
inertscalerot2_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3k::Beam3k(const DRT::ELEMENTS::Beam3k& old) :
 DRT::ELEMENTS::Beam3Base(old),
 isinit_(old.isinit_),
 statmechprob_(old.statmechprob_),
 crosssec_(old.crosssec_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 T0_(old.T0_),
 T_(old.T_),
 theta0_(old.theta0_),
 dispthetanew_(old.dispthetanew_),
 Qnew_(old.Qnew_),
 Qrefconv_(old.Qrefconv_),
 Qrefnew_(old.Qrefnew_),
 K0_(old.K0_),
 length_(old.length_),
 jacobi_(old.jacobi_),
 jacobi2_(old.jacobi2_),
 jacobi_cp_(old.jacobi_cp_),
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
 inertscaletrans_(old.inertscaletrans_),
 inertscalerot1_(old.inertscalerot1_),
 inertscalerot2_(old.inertscalerot2_)
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
DRT::ELEMENTS::Beam3k::~Beam3k()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Print(std::ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3k::Shape() const
{
  int numnodes = NumNode();
  switch(numnodes)
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
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Beam3Base::Pack(data);

  //add all class variables
  AddtoPack(data,isinit_);
  AddtoPack(data,statmechprob_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,Irr_);
  AddtoPack<3,1>(data,T0_);
  AddtoPack<3,1>(data,T_);
  AddtoPack<3,1>(data,theta0_);
  AddtoPack<3,1>(data,dispthetanew_);
  AddtoPack<4,1>(data,Qnew_);
  AddtoPack<4,1>(data,Qrefconv_);
  AddtoPack<4,1>(data,Qrefnew_);
  AddtoPack<3,1>(data,K0_);
  AddtoPack(data,length_);
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobi2_);
  AddtoPack(data,jacobi_cp_);
  AddtoPack(data,rotvec_);
  AddtoPack(data,weakkirchhoff_);
  AddtoPack(data,Eint_);
  AddtoPack(data,Ekin_);
  AddtoPack<4,1>(data,Qconvmass_);
  AddtoPack<4,1>(data,Qnewmass_);
  AddtoPack<3,1>(data,wconvmass_);
  AddtoPack<3,1>(data,wnewmass_);
  AddtoPack<3,1>(data,aconvmass_);
  AddtoPack<3,1>(data,anewmass_);
  AddtoPack<3,1>(data,amodconvmass_);
  AddtoPack<3,1>(data,amodnewmass_);
  AddtoPack<3,1>(data,rttconvmass_);
  AddtoPack<3,1>(data,rttnewmass_);
  AddtoPack<3,1>(data,rttmodconvmass_);
  AddtoPack<3,1>(data,rttmodnewmass_);
  AddtoPack<3,1>(data,rtconvmass_);
  AddtoPack<3,1>(data,rtnewmass_);
  AddtoPack<3,1>(data,rconvmass_);
  AddtoPack<3,1>(data,rnewmass_);
  AddtoPack(data,inertscaletrans_);
  AddtoPack(data,inertscalerot1_);
  AddtoPack(data,inertscalerot2_);

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
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Beam3Base::Unpack(basedata);

  //extract all class variables of beam3 element
  isinit_ = ExtractInt(position,data);
  statmechprob_ = ExtractInt(position,data);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack<3,1>(position,data,T0_);
  ExtractfromPack<3,1>(position,data,T_);
  ExtractfromPack<3,1>(position,data,theta0_);
  ExtractfromPack<3,1>(position,data,dispthetanew_);
  ExtractfromPack<4,1>(position,data,Qnew_);
  ExtractfromPack<4,1>(position,data,Qrefconv_);
  ExtractfromPack<4,1>(position,data,Qrefnew_);
  ExtractfromPack<3,1>(position,data,K0_);
  ExtractfromPack(position,data,length_);
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobi2_);
  ExtractfromPack(position,data,jacobi_cp_);
  rotvec_ = ExtractInt(position,data);
  weakkirchhoff_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Eint_);
  ExtractfromPack(position,data,Ekin_);
  ExtractfromPack<4,1>(position,data,Qconvmass_);
  ExtractfromPack<4,1>(position,data,Qnewmass_);
  ExtractfromPack<3,1>(position,data,wconvmass_);
  ExtractfromPack<3,1>(position,data,wnewmass_);
  ExtractfromPack<3,1>(position,data,aconvmass_);
  ExtractfromPack<3,1>(position,data,anewmass_);
  ExtractfromPack<3,1>(position,data,amodconvmass_);
  ExtractfromPack<3,1>(position,data,amodnewmass_);
  ExtractfromPack<3,1>(position,data,rttconvmass_);
  ExtractfromPack<3,1>(position,data,rttnewmass_);
  ExtractfromPack<3,1>(position,data,rttmodconvmass_);
  ExtractfromPack<3,1>(position,data,rttmodnewmass_);
  ExtractfromPack<3,1>(position,data,rtconvmass_);
  ExtractfromPack<3,1>(position,data,rtnewmass_);
  ExtractfromPack<3,1>(position,data,rconvmass_);
  ExtractfromPack<3,1>(position,data,rnewmass_);
  ExtractfromPack(position,data,inertscaletrans_);
  ExtractfromPack(position,data,inertscalerot1_);
  ExtractfromPack(position,data,inertscalerot2_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3k::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   meier 01/16|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometry(const std::vector<LINALG::Matrix<3,1> >& xrefe, const bool secondinit)
{
  if(weakkirchhoff_)
    SetUpReferenceGeometryWK(xrefe,secondinit);
  else
    SetUpReferenceGeometrySK(xrefe,secondinit);

}

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a weak Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometryWK(const std::vector<LINALG::Matrix<3,1> >& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
   *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
   *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
   *the method and is true the element is initialized another time with respective xrefe and rotrefe;
   *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
   *second initilization in principle (e.g. for periodic boundary conditions*/

  //TODO: It turned out that the approximation of initially curved centerlines via third order Hermite polynomials by simply setting
  //the nodal positions and tangents at the boundary nodes to the corresponding values of the analytical geometry at these nodes yields
  //a worse approximation than that following from a third order Lagrange polynomial interpolation based on nodal positions that coincide
  //with the analytical values at the four element nodes. This worse approximation appears in terms of a larger error in the jacobian,
  //the derivative of the relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial configurations
  //where the error in the initial geometry representation might dominate the discretization error it can be useful to consider strategies
  //which yield in a better Hermite approximation of the initial geometry e.g. by determining the initial nodal values (or the shape
  //function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if(!isinit_ || secondinit)
  {
    const int nnode = 2; //number of nodes

    // set the flag statmechprob_
    statmechprob_ = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(), "STATMECHPROB");

    //Calculate the (initial reference triads) = (initial material triads) at the CPs out of the angles theta0_.
    //So far the initial value for the relative angle is set to zero, i.e.
    //material coordinate system and reference system in the reference configuration coincidence.
    std::vector< LINALG::Matrix<3,3> > Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node],Gref[node]);
    }

    T0_.resize(2);
    T_.resize(2);
    //write initial nodal tangents in extra vector
    for (int i=0;i<3;i++)
    {
      (T0_[0])(i)=(Gref[0])(i,0);
      (T0_[1])(i)=(Gref[1])(i,0);
      (T_[0])(i)=(Gref[0])(i,0);
      (T_[1])(i)=(Gref[1])(i,0);
    }

    //Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

    //Vector holding angle theta of triads
    std::vector<LINALG::Matrix<3,1> > theta_cp;
    theta_cp.resize(BEAM3K_COLLOCATION_POINTS);
    LINALG::Matrix<3,1> theta(true);
    LINALG::Matrix<3,1> theta_s(true);
    LINALG::Matrix<3,3> triad_mat(true); //material triad at gp

    //Resize vectors
    ResizeClassVariables(gausspoints.nquad);

    //calculate the length of the element via Newton iteration
    Calculate_length(xrefe,T0_,LENGTHCALCNEWTONTOL);

    //Matrices to store the function values of the Lagrange shape functions used to interpolate theta
    LINALG::Matrix<1,BEAM3K_COLLOCATION_POINTS> L_i;
    LINALG::Matrix<1,BEAM3K_COLLOCATION_POINTS> L_i_xi;

    //Matrices to store the (derivative of) the Hermite shape functions
    LINALG::Matrix<1,2*nnode> N_i;
    LINALG::Matrix<1,2*nnode> N_i_xi;

    //Matrices to store r and dr/dxi
    LINALG::Matrix<3,1> r;
    LINALG::Matrix<3,1> r_xi;

    //storage index for collocation point
    int ind=0;

    //Calculate initial material triads at the collocation points
    for (int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)
    {
      //colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;

      //Get values of shape functions
      L_i.Clear();
      N_i_xi.Clear();
      DRT::UTILS::shape_function_1D(L_i,xi,Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

      //Determine storage position for the node colpt
      ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

      //current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();

      for (int i=0; i<3; i++)
      {
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+T0_[0](i)*N_i_xi(1)+T0_[1](i)*N_i_xi(3);
      }

      jacobi_cp_[ind]=r_xi.Norm2();

      //rotate (initial reference triad) = (initial material triad) at the interior CPs on tangential line
      //resulting from the Hermite interpolation. This is necessary for initially curved geometries for which
      //the Hermite interpolation does not deliver the exact tangent values of the analytical representation of the
      //initial curve at the CPs.
      if(ind>1)//only for internal CPs
      {
        LINALG::Matrix<3,3> G_aux(true);
        LARGEROTATIONS::CalculateSRTriads<double>(r_xi,Gref[ind],G_aux);
        //rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind]=G_aux;
        LARGEROTATIONS::triadtoquaternion(G_aux,Qrefconv_[ind]);
        LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind],theta0_[ind]);
      }
      else
      {
        LARGEROTATIONS::triadtoquaternion(Gref[ind],Qrefconv_[ind]);
      }
      //Set initial values for the quaternions describing the orientation of the current Newton iteration
      Qnew_[ind]=Qrefconv_[ind];
      Qrefnew_[ind]=Qrefconv_[ind];
    }//(int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)

    //SETUP INTERPOLATION via calculation of difference angle
    for (int colpt=0; colpt<BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      theta_cp[colpt].Clear();
      LARGEROTATIONS::triadtoangleright(theta_cp[colpt],Gref[REFERENCE_NODE],Gref[colpt]);
    }

    //Loop through all GPs and computation of all relevant values at each gp
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i_xi.Clear();
      N_i.Clear();
      DRT::UTILS::shape_function_1D(L_i,xi,Shape());
      DRT::UTILS::shape_function_1D_deriv1(L_i_xi,xi,Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

      //current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      theta.Clear();
      theta_s.Clear();

      for (int i=0; i<3; i++)
      {
        r(i)+=xrefe[0](i)*N_i(0)+xrefe[1](i)*N_i(2)+T0_[0](i)*N_i(1)+T0_[1](i)*N_i(3);
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+T0_[0](i)*N_i_xi(1)+T0_[1](i)*N_i_xi(3);
      }

      //calculate jacobi jacobi_=|r'_0|
      jacobi_[numgp]=r_xi.Norm2();

      //calculate interpolated angle
      for (int i=0; i<BEAM3K_COLLOCATION_POINTS; i++)
      {
        for(int j=0; j<3; j++)
        {
          theta(j)+=L_i(i)*theta_cp[i](j);
          theta_s(j)+=L_i_xi(i)*theta_cp[i](j)/jacobi_[numgp];
        }
      }
      computestrain(theta,theta_s,K0_[numgp]);

      triad_mat.Clear();
      LARGEROTATIONS::angletotriad(theta,Gref[REFERENCE_NODE],triad_mat);
      SetInitialDynamicClassVariables(numgp, triad_mat, r);

    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

    isinit_ = true;
  }//if(!isinit_)

}//DRT::ELEMENTS::Beam3k::SetUpReferenceGeometry()

/*--------------------------------------------------------------------------------------------*
 |  Set up the reference geometry of the case of a strong Kirchhoff constraint     meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK(const std::vector<LINALG::Matrix<3,1> >& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
   *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
   *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
   *the method and is true the element is initialized another time with respective xrefe and rotrefe;
   *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
   *second initilization in principle (e.g. for periodic boundary conditions*/

  //TODO: It turned out that the approximation of initially curved centerlines via third order Hermite polynomials by simply setting
  //the nodal positions and tangents at the boundary nodes to the corresponding values of the analytical geometry at these nodes yields
  //a worse approximation than that following from a third order Lagrange polynomial interpolation based on nodal positions that coincide
  //with the analytical values at the four element nodes. This worse approximation appears in terms of a larger error in the jacobian,
  //the derivative of the relative angle theta_s and finally of the initial curvature K0_. For strongly curved initial configurations
  //where the error in the initial geometry representation might dominate the discretization error it can be useful to consider strategies
  //which yield in a better Hermite approximation of the initial geometry e.g. by determining the initial nodal values (or the shape
  //function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if(!isinit_ || secondinit)
  {
    const int nnode = 2; //number of nodes

    // set the flag statmechprob_
    statmechprob_ = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(), "STATMECHPROB");

    //Calculate the (initial reference triads) = (initial material triads) at the CPs out of the angles theta0_.
    //So far the initial value for the relative angle is set to zero, i.e.
    //material coordinate system and reference system in the reference configuration coincidence.
    std::vector< LINALG::Matrix<3,3> > Gref(BEAM3K_COLLOCATION_POINTS);
    for (int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node],Gref[node]);
    }

    T0_.resize(2);
    T_.resize(2);
    //write initial nodal tangents in extra vector
    for (int i=0;i<3;i++)
    {
      (T0_[0])(i)=(Gref[0])(i,0);
      (T0_[1])(i)=(Gref[1])(i,0);
      (T_[0])(i)=(Gref[0])(i,0);
      (T_[1])(i)=(Gref[1])(i,0);
    }

    //Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

    //Vector holding angle theta of triads
    std::vector<double > phi_cp;
    phi_cp.resize(BEAM3K_COLLOCATION_POINTS);
    double phi(true);
    double phi_s(true);
    LINALG::Matrix<3,3> triad_mat(true); //material triad at gp

    ResizeClassVariables(gausspoints.nquad);

    //calculate the length of the element via Newton iteration
    Calculate_length(xrefe,T0_,LENGTHCALCNEWTONTOL);

    //Matrices to store the function values of the Lagrange shape functions used to interpolate theta
    LINALG::Matrix<1,BEAM3K_COLLOCATION_POINTS> L_i;
    LINALG::Matrix<1,BEAM3K_COLLOCATION_POINTS> L_i_xi;

    //Matrices to store the (derivative of) the Hermite shape functions
    LINALG::Matrix<1,2*nnode> N_i;
    LINALG::Matrix<1,2*nnode> N_i_xi;
    LINALG::Matrix<1,2*nnode> N_i_xixi;

    //Matrices to store r, dr/dxi and d^2r/dxi^2
    LINALG::Matrix<3,1> r;
    LINALG::Matrix<3,1> r_xi;
    LINALG::Matrix<3,1> r_xixi;
    LINALG::Matrix<3,1> r_s;
    LINALG::Matrix<3,1> r_ss;
    //centerline curvature
    LINALG::Matrix<3,1> kappacl;

    //storage index for collocation point
    int ind=0;

    //Calculate initial material triads at the collocation points
    for (int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)
    {
      //colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;

      //Get values of shape functions
      L_i.Clear();
      N_i_xi.Clear();
      DRT::UTILS::shape_function_1D(L_i,xi,Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

      //Determine storage position for the node colpt
      ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

      //current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();

      for (int i=0; i<3; i++)
      {
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+T0_[0](i)*N_i_xi(1)+T0_[1](i)*N_i_xi(3);
      }

      jacobi_cp_[ind]=r_xi.Norm2();

      //rotate (initial reference triad) = (initial material triad) at the interior CPs on tangential line
      //resulting from the Hermite interpolation. This is necessary for initially curved geometries for which
      //the Hermite interpolation does not deliver the exact tangent values of the analytical representation of the
      //initial curve at the CPs.
      if(ind>1)//only for internal CPs
      {
        LINALG::Matrix<3,3> G_aux(true);
        LARGEROTATIONS::CalculateSRTriads<double>(r_xi,Gref[ind],G_aux);
        //rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind]=G_aux;
        LARGEROTATIONS::triadtoquaternion(G_aux,Qrefconv_[ind]);
        LARGEROTATIONS::quaterniontoangle(Qrefconv_[ind],theta0_[ind]);
      }
      else
      {
        LARGEROTATIONS::triadtoquaternion(Gref[ind],Qrefconv_[ind]);
      }
      //Set initial values for the quaternions describing the orientation of the current Newton iteration
      Qnew_[ind]=Qrefconv_[ind];
      Qrefnew_[ind]=Qrefconv_[ind];
    }//(int node=0;node<BEAM3K_COLLOCATION_POINTS;node++)

    //SETUP INTERPOLATION via calculation of difference angle
    for (int colpt=0; colpt<BEAM3K_COLLOCATION_POINTS; colpt++)
    {
      LINALG::Matrix<3,3> Lambdabarref(true);
      LINALG::Matrix<3,1> tangentref(true);
      LINALG::Matrix<3,1> phivec(true);
      for(int i=0;i<3;i++)
      {
        tangentref(i)=Gref[colpt](i,0);
      }
      LARGEROTATIONS::CalculateSRTriads<double>(tangentref,Gref[REFERENCE_NODE],Lambdabarref);
      LARGEROTATIONS::triadtoangleleft(phivec,Lambdabarref,Gref[colpt]);
      phi_cp[colpt]=0.0;
      for(int i=0;i<3;i++)
      {
        phi_cp[colpt]+=tangentref(i)*phivec(i);
      }
    }

    //Loop through all GPs and computation of all relevant values at each gp
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get values of shape functions
      L_i.Clear();
      L_i_xi.Clear();
      N_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();

      DRT::UTILS::shape_function_1D(L_i,xi,Shape());
      DRT::UTILS::shape_function_1D_deriv1(L_i_xi,xi,Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

      //current value of derivatives at GP (derivatives in xi!)
      r.Clear();
      r_xi.Clear();
      r_xixi.Clear();
      r_s.Clear();
      r_ss.Clear();
      kappacl.Clear();
      phi=0.0;
      phi_s=0.0;

      for (int i=0; i<3; i++)
      {
        r(i)+=xrefe[0](i)*N_i(0)+xrefe[1](i)*N_i(2)+T0_[0](i)*N_i(1)+T0_[1](i)*N_i(3);
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+T0_[0](i)*N_i_xi(1)+T0_[1](i)*N_i_xi(3);
        r_xixi(i)+=xrefe[0](i)*N_i_xixi(0)+xrefe[1](i)*N_i_xixi(2)+T0_[0](i)*N_i_xixi(1)+T0_[1](i)*N_i_xixi(3);
      }

      //calculate jacobi_=||r'_0|| and jacobi2_=r'_0^T r''_0
      jacobi_[numgp]=r_xi.Norm2();
      for(int i=0;i<3;i++)
        jacobi2_[numgp]+=r_xi(i)*r_xixi(i);

      //calculate interpolated angle
      for (int i=0; i<BEAM3K_COLLOCATION_POINTS; i++)
      {
        phi+=L_i(i)*phi_cp[i];
        phi_s+=L_i_xi(i)*phi_cp[i]/jacobi_[numgp];
      }

      //calculate derivatives in s
      r_s=r_xi;
      r_s.Scale(1/jacobi_[numgp]);
      for (int i=0; i<3; i++)
      {
        r_ss(i)=r_xixi(i)/pow(jacobi_[numgp],2.0)-r_xi(i)*jacobi2_[numgp]/pow(jacobi_[numgp],4.0);
      }

      triad_mat.Clear();
      ComputeTriadSK(phi,r_s,Gref[REFERENCE_NODE],triad_mat);
      Calculate_clcurvature(r_s, r_ss, kappacl);

      computestrainSK(phi_s,kappacl,Gref[REFERENCE_NODE],triad_mat,K0_[numgp]);

      SetInitialDynamicClassVariables(numgp, triad_mat, r);

    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

    isinit_ = true;
  }//if(!isinit_)

}//DRT::ELEMENTS::Beam3k::SetUpReferenceGeometrySK()

/*--------------------------------------------------------------------------------------------*
 |  Calculates the element length via a Newton Iteration: f(l)=l-int(|N'd|)dxi=0   meier 01/16|
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Calculate_length(const std::vector<LINALG::Matrix<3,1> >& xrefe,
                                              const std::vector<LINALG::Matrix<3,1> >& trefe,
                                              double tolerance)
{
  const int nnode = 2; //number of nodes
  const int vnode = 2; //interpolated values per node (2: value + derivative of value)

  //Get integration points for exact integration
  //DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::intrule_line_10point);

  //Newton Iteration - Tolerance and residual
  double res=1.0;

  //Integral-value for Gauss Integration
  double int_length=0.0;
  //Derivative value of the length integral for Newton Iteration (=weighted sum over deriv_int, gauss quadrature of: int(d/dl(|N'd|))dxi)
  double deriv_length=0.0;
  //value needed to store the derivative of the integral at the GP: d/dl(|N'd|)
  double deriv_int=0.0;

  //inital value for iteration
  {
    LINALG::Matrix<3,1> tempvec;
    tempvec.Clear();
    for(int i=0; i<3; i++)
    {
      tempvec(i)=xrefe[1](i)-xrefe[0](i);
    }
    length_=tempvec.Norm2();
  }

  //Matrices to store the function values of the shape functions
  LINALG::Matrix<1,nnode*vnode> shapefuncderiv;

  shapefuncderiv.Clear();

  //current value of the derivative at the GP
  LINALG::Matrix<3,1> r_xi;

  while(fabs(res)>tolerance)
  {
    int_length=0;
    deriv_length=0;
    //Loop through all GPs and computation of the length and the derivative of the length
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      deriv_int=0;
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get derivatives of the shape functions
        DRT::UTILS::shape_function_hermite_1D_deriv1(shapefuncderiv,xi,length_,line2);

      //integral of the length
      r_xi.Clear();
      deriv_int=0;
      for (int i=0; i<3; i++)
      {
        r_xi(i)+=xrefe[0](i)*shapefuncderiv(0)+xrefe[1](i)*shapefuncderiv(2)+trefe[0](i)*shapefuncderiv(1)+trefe[1](i)*shapefuncderiv(3);
      }
      int_length+=gausspoints.qwgt[numgp]*r_xi.Norm2();

      //derivative of the integral of the length at GP
      for (int i=0; i<3; i++)
      {
          deriv_int+=(trefe[0](i)*shapefuncderiv(1)/length_+trefe[1](i)*shapefuncderiv(3)/length_)*r_xi(i);
      }
      deriv_length+=gausspoints.qwgt[numgp]*deriv_int/r_xi.Norm2();
    }
    //cout << endl << "Länge:" << length_ << "\tIntegral:" << int_length << "\tAbleitung:" << deriv_length << endl;
    res=length_-int_length;
    //Update
    length_=length_-res/(1-deriv_length); //the derivative of f(l)=l-int(|N'd|)dxi=0 is f'(l)=1-int(d/dl(|N'd|))dxi
  }

//  //*************************************begin: Determine optimal Hermite constant for circle segment******************************
//  if(fabs(trefe[0].Norm2()-1.0)>1.0e-12 or fabs(trefe[1].Norm2()-1.0)>1.0e-12)
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
//  std::cout << std::setprecision(16) << copt << " " << length_ << " " << length_approx << " " << length_analyt << std::endl;
//  //*************************************end: Determine optimal Hermite constant for circle segment******************************

  return;
}

double DRT::ELEMENTS::Beam3k::GetJacobiFacAtXi(const double& xi) const
{
  const int nnode=2;

  // Matrices to store the the Hermite shape function derivative values
  LINALG::Matrix<1,2*nnode> N_i_xi;
  DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

  // jacobi = ds/dxi = ||r'_0||
  LINALG::Matrix<3,1> r_xi;

  for (int dim=0; dim<3; dim++)
  {
    r_xi(dim) += Nodes()[0]->X()[dim]*N_i_xi(0) + Nodes()[1]->X()[dim]*N_i_xi(2)
               + T0_[0](dim)*N_i_xi(1) + T0_[1](dim)*N_i_xi(3);
  }

  return r_xi.Norm2();
}

// Todo remove this redundant method; see Beam3Base::Calc_r for general implementation
/*----------------------------------------------------------------------------------------------------------*
 | Get position vector at xi for given nodal displacements                                        popp 02/16|
 *----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3k::GetPos(const double& xi,
                                                  const LINALG::Matrix<12,1>& disp_totlag_centerline) const
{
  // note: this method expects the absolute ("total Lagrangean") values for positions and tangents of both centerline nodes (local numbering 0 and 1)
  LINALG::Matrix<3,1> r(true);
  LINALG::Matrix<4,1> N_i(true);

  DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

  for (int n=0;n<4;n++)
  {
    for (int i=0;i<3;i++)
    {
      r(i)+=N_i(n)*disp_totlag_centerline(3*n+i);
    }
  }

  return (r);
}

/*----------------------------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::GetPosAtXi(LINALG::Matrix<3,1>&       pos,
                                       const double&              xi,
                                       const std::vector<double>& disp) const
{
  /* we expect the (centerline) displacement state vector (positions and tangents at boundary nodes) here;
   * for flexibility, we also accept complete DoF vector of this element (2*7+1 DoFs) and do the rest automatically here
   * NOTE: the latter option is favorable in case of rotvec_==true:
   *       in case of ROTVEC=true, we need 3 positions, 3 absolute rotvec DOFs AND the length of tangent vector (7th DoF) */
  LINALG::Matrix<12,1> disp_totlag_centerline(true);

  if (disp.size() == 15)
  {
    // in this case, we need to "add" reference values first, because if rotvec_==true,
    // we can extract tangent vectors only from total rotation vectors
    std::vector<double> disp_totlag(disp);
    AddRefValuesDisp(disp_totlag);
    this->ExtractCenterlineDofValues<2,2,double>(disp_totlag,disp_totlag_centerline);
  }
  else if (disp.size() == 12)
  {
    /* in this case, we expect the position and tangent DOF values for both boundary nodes;
     * for rotvec_==true, the tangents are NOT nodal DOFs, so they need to be pre-calculated
     * before calling this method */
    disp_totlag_centerline = LINALG::Matrix<12,1>(&disp[0]);
    AddRefValuesDispCenterline<2,2,double>(disp_totlag_centerline);
  }
  else
  {
    dserror("size mismatch: expected either 12 values for disp_totlag_centerline or 15 for "
        "disp_totlag and got %d",disp.size());
  }

  pos = this->GetPos(xi,disp_totlag_centerline);
}

/*----------------------------------------------------------------------------------------------------------*
 | \brief Get base vectors describing the cross-section orientation and size at a given xi       meier 03/16|
 *----------------------------------------------------------------------------------------------------------*/
std::vector<LINALG::Matrix<3,1> > DRT::ELEMENTS::Beam3k::GetBaseVectors(double& xi) const
{
  std::vector<LINALG::Matrix<3,1> > basevectors(2,LINALG::Matrix<3,1>(true));

  if(weakkirchhoff_)
  {
    //Triads at collocation points
    std::vector<LINALG::Matrix<3,3> > triad_mat_cp(BEAM3K_COLLOCATION_POINTS);  //material triads at collocation points
    std::vector<LINALG::Matrix<3,1> > theta_cp(BEAM3K_COLLOCATION_POINTS);  //relative angle at collocation points

    //Interpolated material triad at xi
    LINALG::Matrix<3,3> triad_mat(true); //vector of material triads at xi
    LINALG::Matrix<3,1> theta(true); //vector of material triads at xi

    //Matrices for individual shape functions and xi-derivatives
    LINALG::Matrix<1,BEAM3K_COLLOCATION_POINTS> L_i;

    //calculate angle at cp (this has to be done in a SEPARATE loop as follows)
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      triad_mat_cp[node].Clear();
      LARGEROTATIONS::quaterniontotriad(Qnew_[node],triad_mat_cp[node]);
    }

    //calculate angle at cp (this has to be done in a SEPARATE loop as follows)
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      theta_cp[node].Clear();
      LARGEROTATIONS::triadtoangleright(theta_cp[node],triad_mat_cp[REFERENCE_NODE],triad_mat_cp[node]);
    }

    //Evaluate shape functions
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      theta.Update(L_i(node),theta_cp[node],1.0);
    }

    //compute material triad at gp
    triad_mat.Clear();
    LARGEROTATIONS::angletotriad(theta,triad_mat_cp[REFERENCE_NODE],triad_mat);

    for(int i=0;i<3;i++)
    {
      basevectors[0](i)=triad_mat(i,1);
      basevectors[1](i)=triad_mat(i,2);
    }

    // ToDo careful, this is a hack for rectangular cross-sections ?!?
    double Ly=pow(pow(12.0*Izz_,3)/(12.0*Iyy_),1.0/8.0);
    double Lz=12.0*Izz_/pow(Ly,3);

    basevectors[0].Scale(0.5*Ly);
    basevectors[1].Scale(0.5*Lz);
  }
  else
  {
    dserror("Method GetBaseVectors() only implemented for the case WK==1 (weak Kirchhoff constraint) so far!");
  }

  return basevectors;
}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValues(const std::vector<T>&                   dofvec,
                                                      LINALG::TMatrix<T,3*vpernode*nnodecl,1>& dofvec_centerline) const
{
  // nnodecl: number of nodes used for interpolation of centerline
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  const int dofperboundarynode = 3*vpernode+1;

  // get current values for position DOFs directly
  for (unsigned int dim=0; dim<3; ++dim)
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      dofvec_centerline(3*vpernode*node+dim) = dofvec[dofperboundarynode*node+dim];
    }

  // Hermite interpolation: get values of tangent DOFs as well
  if(rotvec_==false)
  {
    for (unsigned int dim=0; dim<3; ++dim)
      for (unsigned int node=0; node<nnodecl; ++node)
      {
        // in case of rotvec_==false, tangents are equivalent to nodal DOFs
        dofvec_centerline(3*vpernode*node+3+dim) = dofvec[dofperboundarynode*node+3+dim];
      }
  }
  else
  {
    // values for tangent DOFs must be transformed in case of rotvec_==true
    LINALG::TMatrix<T,3,1> theta(true);
    LINALG::TMatrix<T,3,3> unity(true);
    LINALG::TMatrix<T,3,3> triad(true);

    for (unsigned int dim=0; dim<3; ++dim)
      unity(dim,dim)=1.0;

    for (unsigned int node=0; node<nnodecl; ++node)
    {
      for (unsigned int dim=0; dim<3; ++dim)
      {
        // get values for rotation vector DOFs
        theta(dim)=dofvec[dofperboundarynode*node+3+dim];
      }
      // transform to triad
      LARGEROTATIONS::angletotriad(theta,unity,triad);

      // direction of tangent is equivalent to first base vector of triad; length of tangent is 7th nodal DOF
      for (unsigned int dim=0; dim<3; ++dim)
      {
        dofvec_centerline(3*vpernode*node+3+dim) = triad(dim,0) * (dofvec[dofperboundarynode*node+6]);
      }
    }
  }

}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodecl, unsigned int vpernode, typename T>
void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline(LINALG::TMatrix<T,3*vpernode*nnodecl,1>& dofvec_centerline) const
{
  for (unsigned int dim=0; dim<3; ++dim)
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      dofvec_centerline(3*vpernode*node+dim) += Nodes()[node]->X()[dim];

      // have Hermite interpolation? then update tangent DOFs as well
      if(vpernode==2)
        dofvec_centerline(3*vpernode*node+3+dim) += T0_[node](dim);
    }
}

/*------------------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::AddRefValuesDisp(std::vector<double>& dofvec) const
{
  //For the boundary nodes we have to add the initial values. For the interior nodes we are already
  //done since the initial relative angles are zero: alpha_i(t=0)=0;
  if(rotvec_==false)
  {
    //Calculate total displacements = positions vectors, tangents and relative angles
    for (int node=0;node<2;node++)//loop over boundary nodes
    {
      for (int ndof=0;ndof<7;ndof++)//loop over dofs per node
      {
        if(ndof<3)
        {
          dofvec[7*node+ndof]+=Nodes()[node]->X()[ndof];
        }
        else if (ndof<6)
        {
          dofvec[7*node+ndof]+=(Tref()[node])(ndof-3);
        }
        else
        {
          //nothing to do here: alpha_i(t=0)=0;
        }
      }
    }
  }
  else
  {
    //Calculate total displacements = positions vectors, tangents and relative angles
    for (int node=0;node<2;node++)//loop over boundary nodes
    {
      LINALG::Matrix<3,1> nodalangle(true);
      for (int ndof=0;ndof<7;ndof++)//loop over dofs per node
      {
        if(ndof<3)
        {
          dofvec[7*node+ndof]+=Nodes()[node]->X()[ndof];
        }
        else if (ndof<6)
        {
          //Nothing to do here, rotations are treated below
        }
        else
        {
          //here we have to add the initial length of the tangents at the boundary nodes, i.e. ||r'_i(t=0)||=1:
          dofvec[7*node+ndof]+=1.0;
        }
      }//for (int ndof=0;ndof<7;ndof++)//loop over dofs per node

      LINALG::Matrix<3,1> disptheta(true);
      LINALG::Matrix<3,1> thetanew(true);
      LINALG::Matrix<4,1> deltaQ(true);
      LINALG::Matrix<4,1> Qnew(true);
      for(int i=0; i<3; i++)
      {
        disptheta(i) = dofvec[7*node+3+i];
      }

      LARGEROTATIONS::angletoquaternion(disptheta,deltaQ);

      LINALG::Matrix<4,1> Q0;
      LARGEROTATIONS::angletoquaternion(Theta0()[node],Q0);
      LARGEROTATIONS::quaternionproduct(Q0,deltaQ,Qnew);

      //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
      Qnew.Scale(1.0/Qnew.Norm2());

      //Calculate the new nodal angle thetanew \in ]-PI,PI] -> Here, thetanew \in ]-PI,PI] by quaterniontoangle()
      LARGEROTATIONS::quaterniontoangle(Qnew,thetanew);

      //Finally set rotation values in disp_totlag
      for(int i=0; i<3; i++)
      {
        dofvec[7*node+3+i]=thetanew(i);
      }
    }//for (int node=0;node<2;node++)//loop over boundary nodes
  }
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 01/16|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3kType::Initialize(DRT::Discretization& dis)
{
  //setting up geometric variables for Beam3k elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a Beam3k element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    //if we get so far current element is a Beam3k element and  we get a pointer at it
    DRT::ELEMENTS::Beam3k* currele = dynamic_cast<DRT::ELEMENTS::Beam3k*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3k* failed");

    //reference node position
    std::vector<LINALG::Matrix<3,1> > xrefe;

    const int nnode=currele->NumNode();

    //resize xrefe for the number of nodes to store
    xrefe.resize(nnode);

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int node=0; node<nnode; node++) //element has k nodes
        for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
        {
          xrefe[node](dof) = currele->Nodes()[node]->X()[dof];
        }
    }

    //Set up all geometrical (triads, curvatures, jacobians etc.) quantities describing the (initial) reference geometry
    currele->SetUpReferenceGeometry(xrefe);

  } //for (int num=0; num<dis_.NumMyColElements(); ++num)
  return 0;
}

// explicit template instantiations
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValues<2,2,Sacado::Fad::DFad<double> >
                                              (const std::vector<Sacado::Fad::DFad<double> >&,
                                               LINALG::TMatrix<Sacado::Fad::DFad<double>,12,1>&) const;
template void DRT::ELEMENTS::Beam3k::ExtractCenterlineDofValues<2,2,double>
                                              (const std::vector<double>&,
                                               LINALG::TMatrix<double,12,1>&) const;
template void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline<2,2,Sacado::Fad::DFad<double> >
                                              (LINALG::TMatrix<Sacado::Fad::DFad<double>,12,1>&) const;
template void DRT::ELEMENTS::Beam3k::AddRefValuesDispCenterline<2,2,double>
                                              (LINALG::TMatrix<double,12,1>&) const;
