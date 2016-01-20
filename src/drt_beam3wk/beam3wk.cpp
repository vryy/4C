/*!----------------------------------------------------------------------
\file beam3wk.cpp

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "../drt_beam3wk/beam3wk.H"

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

DRT::ELEMENTS::Beam3wkType DRT::ELEMENTS::Beam3wkType::instance_;

DRT::ELEMENTS::Beam3wkType& DRT::ELEMENTS::Beam3wkType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Beam3wkType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3wk* object = new DRT::ELEMENTS::Beam3wk(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3wkType::Create(const std::string eletype,
                                   const std::string eledistype,
                                 const int id,
                                 const int owner )
{
  if ( eletype=="BEAM3WK" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3wk(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3wkType::Create( const int id, const int owner )

{
  return Teuchos::rcp( new Beam3wk( id, owner ) );
}

void DRT::ELEMENTS::Beam3wkType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
      dserror("method 'NodalBlockInformation' not implemented for element type beam3wk!");
}

void DRT::ELEMENTS::Beam3wkType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  dserror("Function not implemented yet.");
}

void DRT::ELEMENTS::Beam3wkType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3WK"];
  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",6)
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("ROTVEC")
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDouble("IT")
    .AddNamedDouble("IR1")
    .AddNamedDouble("IR2")
    .AddNamedDoubleVector("TRIADS",9)
    ;

  defs["LINE4"]
      .AddIntVector("LINE4",4)
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddNamedDouble("IT")
      .AddNamedDouble("IR1")
      .AddNamedDouble("IR2")
      .AddNamedDoubleVector("TRIADS",12)
      ;

    defs["LIN4"]
      .AddIntVector("LIN4",4)
      .AddNamedInt("ROTVEC")
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddNamedDouble("IT")
      .AddNamedDouble("IR1")
      .AddNamedDouble("IR2")
      .AddNamedDoubleVector("TRIADS",12)
      ;

}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3wk::Beam3wk(int id, int owner) :
DRT::Element(id,owner),
isinit_(false),
crosssec_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
qrefconv_(0),
theta0_(0),
kappa0_(0),
length_(0),
jacobi_(0),
jacobi_cp_(0),
rotvec_(false),
firstcall_(true),
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
DRT::ELEMENTS::Beam3wk::Beam3wk(const DRT::ELEMENTS::Beam3wk& old) :
 DRT::Element(old),
 isinit_(old.isinit_),
 crosssec_(old.crosssec_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 qrefconv_(old.qrefconv_),
 theta0_(old.theta0_),
 kappa0_(old.kappa0_),
 length_(old.length_),
 jacobi_(old.jacobi_),
 jacobi_cp_(old.jacobi_cp_),
 rotvec_(old.rotvec_),
 firstcall_(old.firstcall_),
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
 |  Deep copy this instance of Beam3wk and return pointer to it (public) |
 |                                                                    meier 05/12 |
 *--------------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3wk::Clone() const
{
  DRT::ELEMENTS::Beam3wk* newelement = new DRT::ELEMENTS::Beam3wk(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3wk::~Beam3wk()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::Print(std::ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3wk::Shape() const
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
void DRT::ELEMENTS::Beam3wk::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  //add all class variables of beam2r element

  AddtoPack(data,isinit_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,Irr_);
  AddtoPack<4,1>(data,qrefconv_);
  AddtoPack<3,1>(data,theta0_);
  AddtoPack<3,1>(data,kappa0_);
  AddtoPack(data,length_);
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobi_cp_);
  AddtoPack(data,rotvec_);
  AddtoPack(data,firstcall_);
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
void DRT::ELEMENTS::Beam3wk::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  std::cout<<"We are in Unpack method"<<std::endl;
  //extract all class variables of beam3 element

  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack<4,1>(position,data,qrefconv_);
  ExtractfromPack<3,1>(position,data,theta0_);
  ExtractfromPack<3,1>(position,data,kappa0_);
  ExtractfromPack(position,data,length_);
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobi_cp_);
  rotvec_ = ExtractInt(position,data);
  firstcall_ = ExtractInt(position,data);
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
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3wk::Lines()
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
 | element nodes are already known (public)                   meier 05/12|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Beam3wk::SetUpReferenceGeometry(const std::vector<LINALG::Matrix<3,1> >& xrefe, const bool secondinit)
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
  //the derivative of the relative angle theta_s and finally of the initial curvature kapparef_. For strongly curved initial configurations
  //where the error in the initial geometry representation might dominate the discretization error it can be useful to consider strategies
  //which yield in a better Hermite approximation of the initial geometry e.g. by determining the initial nodal values (or the shape
  //function parameter c=c_{opt} instead of c=l) based on optimization strategies.

  if(!isinit_ || secondinit)
  {
    const int nnode = 2; //number of nodes

    //Calculate the (initial reference triads) = (initial material triads) at the CPs out of the angles theta0_.
    //So far the initial value for the relative angle is set to zero, i.e.
    //material coordinate system and reference system in the reference configuration coincidence.
    std::vector< LINALG::Matrix<3,3> > Gref(COLLOCATION_POINTS);
    for (int node=0;node<COLLOCATION_POINTS;node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node],Gref[node]);
    }

    //write initial nodal tangents in extra vector
    std::vector<LINALG::Matrix<3,1> > trefe(2);
    for (int i=0;i<3;i++)
    {
      (trefe[0])(i)=(Gref[0])(i,0);
      (trefe[1])(i)=(Gref[1])(i,0);
    }

    //Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3wk);

    //Vector holding angle theta of triads
    std::vector<LINALG::Matrix<3,1> > theta_cp;
    theta_cp.resize(COLLOCATION_POINTS);
    LINALG::Matrix<3,1> theta(true);
    LINALG::Matrix<3,1> theta_s(true);
    LINALG::Matrix<3,3> triad_mat(true); //material triad at gp

    //Resize vectors
    qrefconv_.resize(COLLOCATION_POINTS);
    kappa0_.resize(gausspoints.nquad);
    jacobi_.resize(gausspoints.nquad);
    jacobi_cp_.resize(COLLOCATION_POINTS);

    Qconvmass_.resize(gausspoints.nquad);
    Qnewmass_.resize(gausspoints.nquad);
    wconvmass_.resize(gausspoints.nquad);
    wnewmass_.resize(gausspoints.nquad);
    aconvmass_.resize(gausspoints.nquad);
    anewmass_.resize(gausspoints.nquad);
    amodconvmass_.resize(gausspoints.nquad);
    amodnewmass_.resize(gausspoints.nquad);
    rttconvmass_.resize(gausspoints.nquad);
    rttnewmass_.resize(gausspoints.nquad);
    rttmodconvmass_.resize(gausspoints.nquad);
    rttmodnewmass_.resize(gausspoints.nquad);
    rtconvmass_.resize(gausspoints.nquad);
    rtnewmass_.resize(gausspoints.nquad);
    rconvmass_.resize(gausspoints.nquad);
    rnewmass_.resize(gausspoints.nquad);

    //calculate the length of the element via Newton iteration
    calculate_length(xrefe,trefe,LENGTHCALCNEWTONTOL);

    //Matrices to store the function values of the Lagrange shape functions used to interpolate theta
    LINALG::Matrix<1,COLLOCATION_POINTS> L_i;
    LINALG::Matrix<1,COLLOCATION_POINTS> L_i_xi;

    //Matrices to store the (derivative of) the Hermite shape functions
    LINALG::Matrix<1,2*nnode> N_i;
    LINALG::Matrix<1,2*nnode> N_i_xi;

    //Matrices to store r and dr/dxi
    LINALG::Matrix<3,1> r;
    LINALG::Matrix<3,1> r_xi;

    //storage index for collocation point
    int ind=0;

    //Calculate initial material triads at the collocation points
    for (int node=0;node<COLLOCATION_POINTS;node++)
    {
      //colpt=0->xi=-1  colpt=1->xi=0 colpt=2->xi=1
      const double xi=(double)node/(COLLOCATION_POINTS-1)*2-1.0;

      //Get values of shape functions
      L_i.Clear();
      N_i_xi.Clear();
      DRT::UTILS::shape_function_1D(L_i,xi,Shape());
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

      //Determine storage position for the node colpt
      ind=LARGEROTATIONS::NumberingTrafo(node+1, COLLOCATION_POINTS);

      //current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();

      for (int i=0; i<3; i++)
      {
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+trefe[0](i)*N_i_xi(1)+trefe[1](i)*N_i_xi(3);
      }

      jacobi_cp_[ind]=r_xi.Norm2();

      //rotate (initial reference triad) = (initial material triad) at the interior CPs on tangential line
      //resulting from the Hermite interpolation. This is necessary for initially curved geometries for which
      //the Hermite interpolation does not deliver the exact tangent values of the analytical representation of the
      //initial curve at the CPs.
      if(ind>1)//only for internal CPs
      {
        LINALG::Matrix<3,3> G_aux(true);
        CalculateSRTriads<double>(r_xi,Gref[ind],G_aux);
        //rotate also Gref and theta0_ via smallest rotation to get a consistent initial state
        Gref[ind]=G_aux;
        LARGEROTATIONS::triadtoquaternion(G_aux,qrefconv_[ind]);
        LARGEROTATIONS::quaterniontoangle(qrefconv_[ind],theta0_[ind]);
      }
      else
      {
        LARGEROTATIONS::triadtoquaternion(Gref[ind],qrefconv_[ind]);
      }
    }//(int node=0;node<COLLOCATION_POINTS;node++)

    //SETUP INTERPOLATION via calculation of difference angle
    for (int colpt=0; colpt<COLLOCATION_POINTS; colpt++)
    {
      LINALG::Matrix<3,3> Lambda(true);
      theta_cp[colpt].Clear();
      triadtoangle(theta_cp[colpt],Gref[REFERENCE_NODE],Gref[colpt]);
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
        r(i)+=xrefe[0](i)*N_i(0)+xrefe[1](i)*N_i(2)+trefe[0](i)*N_i(1)+trefe[1](i)*N_i(3);
        r_xi(i)+=xrefe[0](i)*N_i_xi(0)+xrefe[1](i)*N_i_xi(2)+trefe[0](i)*N_i_xi(1)+trefe[1](i)*N_i_xi(3);
      }

      //calculate jacobi jacobi_=|r'_0|
      jacobi_[numgp]=r_xi.Norm2();

      //calculate interpolated angle
      for (int i=0; i<COLLOCATION_POINTS; i++)
      {
        for(int j=0; j<3; j++)
        {
          theta(j)+=L_i(i)*theta_cp[i](j);
          theta_s(j)+=L_i_xi(i)*theta_cp[i](j)/jacobi_[numgp];
        }
      }
      computestrain(theta,theta_s,kappa0_[numgp]);

      Qconvmass_[numgp].Clear();
      Qnewmass_[numgp].Clear();
      rconvmass_[numgp].Clear();
      rnewmass_[numgp].Clear();

      wconvmass_[numgp].Clear();
      wnewmass_[numgp].Clear();
      aconvmass_[numgp].Clear();
      anewmass_[numgp].Clear();
      amodconvmass_[numgp].Clear();
      amodnewmass_[numgp].Clear();
      rtconvmass_[numgp].Clear();
      rtnewmass_[numgp].Clear();
      rttconvmass_[numgp].Clear();
      rttnewmass_[numgp].Clear();
      rttmodconvmass_[numgp].Clear();
      rttmodnewmass_[numgp].Clear();

      //compute material triad at gp
      triad_mat.Clear();
      angletotriad(theta,Gref[REFERENCE_NODE],triad_mat);
      LARGEROTATIONS::triadtoquaternion(triad_mat,Qconvmass_[numgp]);
      Qnewmass_[numgp]=Qconvmass_[numgp];
      rconvmass_[numgp]=r;
      rnewmass_[numgp]=r;


    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

    isinit_ = true;
  }//if(!isinit_)

}//DRT::ELEMENTS::Beam3wk::SetUpReferenceGeometry()

//Calculates the length of the element based upon a Newton Iteration
//f(l)=l-int(|N'd|)dxi=0
void DRT::ELEMENTS::Beam3wk::calculate_length(const std::vector<LINALG::Matrix<3,1> >& xrefe,const std::vector<LINALG::Matrix<3,1> >& trefe,double tolerance)
{
  const int nnode = 2; //number of nodes
  const int vnode = 2; //interpolated values per node (2: value + derivative of value)

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

  //Get integration points for exact integration
  //DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3wk);
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::intrule_line_10point);

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

  return;
}

LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3wk::calculate_curvature(const LINALG::Matrix<3,1>& r_s, const LINALG::Matrix<3,1>& r_ss)
{
  LINALG::Matrix<3,1> curvature;

  //spinmatrix Sr' = r'x
  LINALG::Matrix<3,3> Srx;
  LARGEROTATIONS::computespin(Srx,r_s);

  //cross-product r'xr''
  LINALG::Matrix<3,1> Srxrxx;
  Srxrxx.Clear();
  Srxrxx.Multiply(Srx,r_ss);

  for (int i=0;i<3;i++)
  {
    curvature(i)=Srxrxx(i)/pow(r_s.Norm2(),2);
  }

  return curvature;
}

std::vector<double> DRT::ELEMENTS::Beam3wk::calculate_jacobi(const LINALG::Matrix<3,1>& r_xi, const LINALG::Matrix<3,1>& r_xixi, const LINALG::Matrix<3,1>& r_xixixi)
{
  std::vector<double> jacobi_(3);
  double scalar1=0.0;
  double scalar2=0.0;

  //|r'(xi)|
  jacobi_[0]=r_xi.Norm2();

  //r'^T(xi)r''(xi)
  jacobi_[1]=0;
  for (int i=0; i<3; i++)
  {
    jacobi_[1]+=r_xi(i)*r_xixi(i);
  }

  //r''^T(xi)r''^T(xi)+r'^T(xi)r'''^T(xi)
  jacobi_[2]=0;
  for (int i=0; i<3; i++)
  {
    scalar1+=r_xixi(i)*r_xixi(i);
    scalar2+=r_xi(i)*r_xixixi(i);
  }
  jacobi_[2]=scalar1+scalar2;
  return jacobi_;
}

void DRT::ELEMENTS::Beam3wk::triadtoangle(LINALG::Matrix<3,1>& theta,LINALG::Matrix<3,3> triad_ref,LINALG::Matrix<3,3> triad)
{
  LINALG::Matrix<3,3> rotation_matrix;
  LINALG::Matrix<4,1> quaternion;

  rotation_matrix.Clear();
  quaternion.Clear();

  rotation_matrix.MultiplyTN(triad_ref,triad);
  LARGEROTATIONS::triadtoquaternion(rotation_matrix,quaternion);
  LARGEROTATIONS::quaterniontoangle(quaternion,theta);
}

void DRT::ELEMENTS::Beam3wk::triadtoangle(LINALG::TMatrix<FAD,3,1>& theta,LINALG::TMatrix<FAD,3,3> triad_ref,LINALG::TMatrix<FAD,3,3> triad)
{
    TEUCHOS_FUNC_TIME_MONITOR("Beam3wk::triadtoangle");
  LINALG::TMatrix<FAD,3,3> rotation_matrix;
  LINALG::TMatrix<FAD,4,1> quaternion;

  rotation_matrix.Clear();
  quaternion.Clear();

  rotation_matrix.MultiplyTN(triad_ref,triad);
  triadtoquaternion(rotation_matrix,quaternion);
  quaterniontoangle(quaternion,theta);
}


void DRT::ELEMENTS::Beam3wk::angletotriad(LINALG::TMatrix<FAD,3,1> theta,LINALG::TMatrix<FAD,3,3> triad_ref,LINALG::TMatrix<FAD,3,3>& triad)
{
  TEUCHOS_FUNC_TIME_MONITOR("Beam3wk::angletotriad");
  LINALG::TMatrix<FAD,3,3> rotation_matrix;
  LINALG::TMatrix<FAD,4,1> quaternion;

  rotation_matrix.Clear();
  quaternion.Clear();
  triad.Clear();

  angletoquaternion(theta,quaternion);
  quaterniontotriad(quaternion,rotation_matrix);
  triad.Multiply(triad_ref,rotation_matrix);// ???
}

void DRT::ELEMENTS::Beam3wk::angletotriad(LINALG::Matrix<3,1> theta,LINALG::Matrix<3,3> triad_ref,LINALG::Matrix<3,3>& triad)
{
    TEUCHOS_FUNC_TIME_MONITOR("Beam3wk::angletotriad");
  LINALG::Matrix<3,3> rotation_matrix;
  LINALG::Matrix<4,1> quaternion;

  rotation_matrix.Clear();
  quaternion.Clear();
  triad.Clear();

  angletoquaternion(theta,quaternion);
  quaterniontotriad(quaternion,rotation_matrix);
  triad.Multiply(triad_ref,rotation_matrix);// ???
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3wkType::Initialize(DRT::Discretization& dis)
{
  //setting up geometric variables for Beam3wk elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a Beam3wk element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    //if we get so far current element is a Beam3wk element and  we get a pointer at it
    DRT::ELEMENTS::Beam3wk* currele = dynamic_cast<DRT::ELEMENTS::Beam3wk*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3wk* failed");

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

    currele->SetUpReferenceGeometry(xrefe);

  } //for (int num=0; num<dis_.NumMyColElements(); ++num)
  return 0;
}
