/*!----------------------------------------------------------------------
\file beam3eb_anisotrop.cpp

\brief three dimensional nonlinear rod based on a C1 curve

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


*-----------------------------------------------------------------------------------------------------------*/

#include "beam3eb_anisotrop.H"
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


DRT::ELEMENTS::Beam3ebanisotropType DRT::ELEMENTS::Beam3ebanisotropType::instance_;

DRT::ELEMENTS::Beam3ebanisotropType & DRT::ELEMENTS::Beam3ebanisotropType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Beam3ebanisotropType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam3ebanisotrop* object = new DRT::ELEMENTS::Beam3ebanisotrop(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebanisotropType::Create(const std::string eletype,
                                                                       const std::string eledistype,
                                                                       const int id,
                                                                       const int owner )
{
  if ( eletype=="BEAM3EBANISOTROP" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam3ebanisotrop(id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam3ebanisotropType::Create( const int id, const int owner )

{
  return Teuchos::rcp( new Beam3ebanisotrop( id, owner ) );
}

void DRT::ELEMENTS::Beam3ebanisotropType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
      dserror("method 'NodalBlockInformation' not implemented for element type beam3eb_anisotrop!");
//    numdf = 6 + TWISTDOFPN;
//    nv = 6 + TWISTDOFPN;
//    dimns = 4;
}

void DRT::ELEMENTS::Beam3ebanisotropType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  dserror("Function not implemented yet.");
  //DRT::UTILS::ComputeXFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam3ebanisotropType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM3EBANISOTROP"];
  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDoubleVector("TANGENTS",6)
    .AddNamedDoubleVector("CURVATURE",6)
    .AddNamedDouble("OPT")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDoubleVector("TANGENTS",6)
    .AddNamedDoubleVector("CURVATURE",6)
    .AddNamedDouble("OPT")
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDoubleVector("TANGENTS",6)
    .AddNamedDoubleVector("CURVATURE",6)
    .AddNamedDouble("OPT")
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("MOMINY")
    .AddNamedDouble("MOMINZ")
    .AddNamedDouble("MOMINPOL")
    .AddNamedDoubleVector("TANGENTS",6)
    .AddNamedDoubleVector("CURVATURE",6)
    .AddNamedDouble("OPT")
    ;

  defs["LINE4"]
      .AddIntVector("LINE4",4)
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddNamedDoubleVector("TANGENTS",6)
      .AddNamedDoubleVector("CURVATURE",6)
      .AddNamedDouble("OPT")
      ;

    defs["LIN4"]
      .AddIntVector("LIN4",4)
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedDouble("MOMINY")
      .AddNamedDouble("MOMINZ")
      .AddNamedDouble("MOMINPOL")
      .AddNamedDoubleVector("TANGENTS",6)
      .AddNamedDoubleVector("CURVATURE",6)
      .AddNamedDouble("OPT")
      ;

}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebanisotrop::Beam3ebanisotrop(int id, int owner) :
 DRT::ELEMENTS::Beam3Base(id,owner),
isinit_(false),
crosssec_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
int_energy_(0),
firstcall_(true),
timestepcount_(0)
{
  theta_nodes_old_.resize(2,0.0);
  sr_theta_nodes_old_=0.0;
  w0_nodes_.resize(2);
  w0_nodes_[0].Clear();
  w0_nodes_[1].Clear();
  dw0perpdt_nodes_.resize(2);
  dw0perpdt_nodes_[0].Clear();
  dw0perpdt_nodes_[1].Clear();
  dw0paralleldt_nodes_.resize(2,0.0);

  w_nodes_.resize(2);
  w_nodes_[0].Clear();
  w_nodes_[1].Clear();
  dwperpdt_nodes_.resize(2);
  dwperpdt_nodes_[0].Clear();
  dwperpdt_nodes_[1].Clear();
  dwparalleldt_nodes_.resize(2,0.0);
  alphat_nodes_.resize(2,0.0);
  alphatt_nodes_.resize(2,0.0);

  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebanisotrop::Beam3ebanisotrop(const DRT::ELEMENTS::Beam3ebanisotrop& old) :
 DRT::ELEMENTS::Beam3Base(old),
 isinit_(old.isinit_),
 crosssec_(old.crosssec_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),

 t0_(old.t0_),
 n0_(old.n0_),
 b0_(old.b0_),
 dt0ds_(old.dt0ds_),
 db0ds_(old.db0ds_),

 Tref_(old.Tref_),
 G2ref_(old.G2ref_),
 kappa0_(old.kappa0_),
 kappa0g20_(old.kappa0g20_),
 kappa0g30_(old.kappa0g30_),
 tau0_(old.tau0_),
 length_(old.length_),
 copt_(old.copt_),
 jacobi_(old.jacobi_),
 jacobi2_(old.jacobi2_),
 jacobi3_(old.jacobi3_),
 int_energy_(old.int_energy_)
{
  return;
}
/*--------------------------------------------------------------------------------*
 |  Deep copy this instance of Beam3ebanisotrop and return pointer to it (public) |
 |                                                                    meier 05/12 |
 *--------------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3ebanisotrop::Clone() const
{
  DRT::ELEMENTS::Beam3ebanisotrop* newelement = new DRT::ELEMENTS::Beam3ebanisotrop(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ebanisotrop::~Beam3ebanisotrop()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::Print(std::ostream& os) const
{
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3ebanisotrop::Shape() const
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
void DRT::ELEMENTS::Beam3ebanisotrop::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  //add all class variables
  AddtoPack(data,jacobi_);
  AddtoPack<3,1>(data,Tref_);
  AddtoPack<3,1>(data,G2ref_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,isinit_);
  AddtoPack(data,Irr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,jacobi2_);
  AddtoPack(data,jacobi3_);
  AddtoPack(data,length_);
  AddtoPack(data,copt_);
  AddtoPack(data,kappa0_);
  AddtoPack(data,tau0_);
  AddtoPack<3,1>(data,t0_);
  AddtoPack<3,1>(data,n0_);
  AddtoPack<3,1>(data,b0_);
  AddtoPack<3,1>(data,dt0ds_);
  AddtoPack<3,1>(data,db0ds_);
  AddtoPack<3,1>(data,t0_nodes_);
  AddtoPack<3,1>(data,n0_nodes_);
  AddtoPack<3,1>(data,b0_nodes_);
  AddtoPack<3,1>(data,dt0ds_nodes_);
  AddtoPack<3,1>(data,db0ds_nodes_);
  AddtoPack<3,1>(data,w0_nodes_);
  AddtoPack<3,1>(data,dw0perpdt_nodes_);
  AddtoPack(data,dw0paralleldt_nodes_);
  AddtoPack<3,1>(data,w_nodes_);
  AddtoPack<3,1>(data,dwperpdt_nodes_);
  AddtoPack(data,dwparalleldt_nodes_);
  AddtoPack(data,alphat_nodes_);
  AddtoPack(data,alphatt_nodes_);
  AddtoPack(data,firstcall_);
  AddtoPack(data,timestepcount_);
  AddtoPack(data,int_energy_);
  AddtoPack(data,sr_theta_nodes_old_);
  AddtoPack(data,theta_nodes_old_);
  AddtoPack(data,kappa0g20_);
  AddtoPack(data,kappa0g30_);
  AddtoPack<3,1>(data,Tref_);
  AddtoPack<3,1>(data,G2ref_);

  return;
}
/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::Unpack(const std::vector<char>& data)
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

  //extract all class variables of beam3 element
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack<3,1>(position,data,Tref_);
  ExtractfromPack<3,1>(position,data,G2ref_);
  ExtractfromPack(position,data,crosssec_);
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,jacobi2_);
  ExtractfromPack(position,data,jacobi3_);
  ExtractfromPack(position,data,length_);
  ExtractfromPack(position,data,copt_);
  ExtractfromPack(position,data,kappa0_);
  ExtractfromPack(position,data,tau0_);
  ExtractfromPack<3,1>(position,data,t0_);
  ExtractfromPack<3,1>(position,data,n0_);
  ExtractfromPack<3,1>(position,data,b0_);
  ExtractfromPack<3,1>(position,data,dt0ds_);
  ExtractfromPack<3,1>(position,data,db0ds_);
  ExtractfromPack<3,1>(position,data,t0_nodes_);
  ExtractfromPack<3,1>(position,data,n0_nodes_);
  ExtractfromPack<3,1>(position,data,b0_nodes_);
  ExtractfromPack<3,1>(position,data,dt0ds_nodes_);
  ExtractfromPack<3,1>(position,data,db0ds_nodes_);
  ExtractfromPack<3,1>(position,data,w0_nodes_);
  ExtractfromPack<3,1>(position,data,dw0perpdt_nodes_);
  ExtractfromPack(position,data,dw0paralleldt_nodes_);
  ExtractfromPack<3,1>(position,data,w_nodes_);
  ExtractfromPack<3,1>(position,data,dwperpdt_nodes_);
  ExtractfromPack(position,data,dwparalleldt_nodes_);
  ExtractfromPack(position,data,alphat_nodes_);
  ExtractfromPack(position,data,alphatt_nodes_);
  firstcall_ = ExtractInt(position,data);
  ExtractfromPack(position,data,timestepcount_);
  ExtractfromPack(position,data,int_energy_);
  ExtractfromPack(position,data,sr_theta_nodes_old_);
  ExtractfromPack(position,data,theta_nodes_old_);
  ExtractfromPack(position,data,kappa0g20_);
  ExtractfromPack(position,data,kappa0g30_);
  ExtractfromPack<3,1>(position,data,Tref_);
  ExtractfromPack<3,1>(position,data,G2ref_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          meier 05/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam3ebanisotrop::Lines()
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

void DRT::ELEMENTS::Beam3ebanisotrop::SetUpReferenceGeometry(const std::vector<LINALG::Matrix<3,1> >& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
   *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
   *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
   *the method and is true the element is initialized another time with respective xrefe and rotrefe;
   *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
   *second initilization in principle (e.g. for periodic boundary conditions*/

  const int nnode = 2; //number of nodes
  const int vnode = 2; //interpolated values per node (2: value + derivative of value)

  //matrix for current nodal positions and nodal tangents
  std::vector<FAD> disp_totlag(6*nnode + TWISTDOFS, 0.0);

  //set disp_totlag for reference geometry: So far the initial value for the relative angle gamma is set to zero, i.e.
  //material coordinate system and reference system in the reference configuration coincidence: gamma_0(s)=0

  for (int i=0;i<3;i++)
  {
    disp_totlag[i]=xrefe[0](i);
    disp_totlag[7 + i]=xrefe[1](i);
    disp_totlag[3+i]=Tref_[0](i);
    disp_totlag[3+7+i]=Tref_[1](i);
  }

  if(!isinit_ || secondinit)
  {

    std::cout << "SetUpReferenceGeometry!!!" << std::endl;

    //Get integration points for exact integration
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebanisotrop);

    //Resize vectors
    t0_.resize(gausspoints.nquad);
    n0_.resize(gausspoints.nquad);
    b0_.resize(gausspoints.nquad);
    dt0ds_.resize(gausspoints.nquad);
    db0ds_.resize(gausspoints.nquad);
    t0_nodes_.resize(2);
    n0_nodes_.resize(2);
    b0_nodes_.resize(2);
    dt0ds_nodes_.resize(2);
    db0ds_nodes_.resize(2);

    kappa0_.resize(gausspoints.nquad);
    kappa0g20_.resize(gausspoints.nquad);
    kappa0g30_.resize(gausspoints.nquad);
    tau0_.resize(gausspoints.nquad);
    jacobi_.resize(gausspoints.nquad);
    jacobi2_.resize(gausspoints.nquad);
    jacobi3_.resize(gausspoints.nquad);

    //calculate the length of the element via Newton iteration
    calculate_length(xrefe, LENGTHCALCNEWTONTOL);

    //Matrices to store the function values of the shape functions
    LINALG::Matrix<1,vnode*nnode> shapefuncderiv;
    LINALG::Matrix<1,vnode*nnode> shapefuncderiv2;
    LINALG::Matrix<1,vnode*nnode> shapefuncderiv3;

    //Matrices to store r, r' and r'' at gp, where (...)'=d/ds
    LINALG::Matrix<3,1> r_s;
    LINALG::Matrix<3,1> r_ss;
    LINALG::Matrix<3,1> r_sss;

    //Matrices to store dr/dxi, d2r/dxi2 and d3r/dxi3
    LINALG::Matrix<3,1> r_xi;
    LINALG::Matrix<3,1> r_xixi;
    LINALG::Matrix<3,1> r_xixixi;

    //Frenet-Serret curvature vector
    LINALG::Matrix<3,1> kappa_vec;

    //Calculate nodal reference triad and derivative
    for (int node=0;node<2;node++)
    {
      LINALG::TMatrix<FAD,3,3> Stangent;
      LINALG::TMatrix<FAD,3,1> tangent;
      Stangent.Clear();
      tangent.Clear();
      for (int i=0;i<3;i++)
      {
        tangent(i)=Tref_[node](i);
      }
      LARGEROTATIONS::computespin<FAD>(Stangent,tangent);
      //Extract material triads at elements nodes
      for (int i=0;i<3;i++)
      {
        dt0ds_nodes_[node](i)=0.0;
        db0ds_nodes_[node](i)=0.0;
        t0_nodes_[node](i)=Tref_[node](i);
        n0_nodes_[node](i)=G2ref_[node](i);
        b0_nodes_[node](i)=0.0;
        for (int j=0;j<3;j++)
        {
          b0_nodes_[node](i)+=Stangent(i,j).val()*G2ref_[node](j);
        }
      }
    }

    //Quantities which are needed for NSRISR calculation
    //attention: The NSRISR interpolation is applied in order to create a C0-continuous initial triad field also when the
    //pure SR or VP formulation is used.

    //difference angle between triad_ref and triad_bar at the element nodes (For the currently implemented NSRISR formulation,
    //the angle theta_nodes[0] is not needed, but it will for example be needed for the NSRIFS method)
    std::vector<FAD> theta_nodes(2,0.0);
    //material triads at the nodes
    LINALG::TMatrix<FAD,3,6> triads_mat_nodes;
    //reference triads at the nodes
    std::vector<LINALG::TMatrix<FAD,2,3> > triads_ref_nodes;
    //material triad at the gauss points
    LINALG::TMatrix<FAD,2,3> triad_mat_gp;
    //intermediate triad at the gauss points
    LINALG::TMatrix<FAD,2,3> triad_bar_gp;
    //intermediate torsion at the gauss points
    FAD tau_bar_gp;
    //reference torsion at the gauss points
    FAD tau_gp;

    triads_mat_nodes.Clear();
    triads_ref_nodes.resize(2);
    triads_ref_nodes[0].Clear();
    triads_ref_nodes[1].Clear();

    FAD taubar0;
    taubar0=0.0;

    //Calculate the nodal reference triads as well as the difference angle between the reference triad triad_ref_nodes[1]
    // and the intermediate triad (normal_bar, binormal_bar) at the right element node
    DetermineNodalTriads(disp_totlag, triads_mat_nodes, theta_nodes,true, triads_ref_nodes,taubar0);

    //Set initial values for the difference angles
    theta_nodes_old_[0]=theta_nodes[0].val();
    theta_nodes_old_[1]=theta_nodes[1].val();
    sr_theta_nodes_old_=theta_nodes[1].val();

    //Loop through all GPs and computation of all relevant values at each gp
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];

      //Get derivatives of the shape functions
      DRT::UTILS::shape_function_hermite_1D_deriv1(shapefuncderiv,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv2(shapefuncderiv2,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv3(shapefuncderiv3,xi,length_,line2);

      //current value of derivatives at GP (derivatives in xi!)
      r_xi.Clear();
      r_xixi.Clear();
      r_xixixi.Clear();
      kappa_vec.Clear();
      triad_mat_gp.Clear();
      tau_gp = 0.0;
      triad_bar_gp.Clear();
      tau_bar_gp = 0.0;

      for (int i=0; i<3; i++)
      {
        r_xi(i)+=xrefe[0](i)*shapefuncderiv(0)+xrefe[1](i)*shapefuncderiv(2)+Tref_[0](i)*shapefuncderiv(1)+Tref_[1](i)*shapefuncderiv(3);
        r_xixi(i)+=xrefe[0](i)*shapefuncderiv2(0)+xrefe[1](i)*shapefuncderiv2(2)+Tref_[0](i)*shapefuncderiv2(1)+Tref_[1](i)*shapefuncderiv2(3);
        r_xixixi(i)+=xrefe[0](i)*shapefuncderiv3(0)+xrefe[1](i)*shapefuncderiv3(2)+Tref_[0](i)*shapefuncderiv3(1)+Tref_[1](i)*shapefuncderiv3(3);
      }

      //calculate jacobi jacobi_=|r'_0| jacobi2_=(r'_0)^T(r''_0) jacobi3_=(r''_0)^T(r''_0)+(r'_0)^T(r'''_0)
      std::vector<double> jacobi(3);
      jacobi=calculate_jacobi(r_xi,r_xixi,r_xixixi);
      jacobi_[numgp]=jacobi[0];
      jacobi2_[numgp]=jacobi[1];
      jacobi3_[numgp]=jacobi[2];

      //calculate derivatives in s
      r_s=r_xi;
      r_s.Scale(1/jacobi_[numgp]);
      for (int i=0; i<3; i++)
      {
        r_ss(i)=r_xixi(i)/pow(jacobi_[numgp],2.0)-r_xi(i)*jacobi2_[numgp]/pow(jacobi_[numgp],4.0);
        r_sss(i)=r_xixixi(i)/pow(jacobi_[numgp],3.0)-3.0*r_xixi(i)*jacobi2_[numgp]/pow(jacobi_[numgp],5.0)-r_xi(i)*jacobi3_[numgp]/pow(jacobi_[numgp],5)+4.0*r_xi(i)*pow(jacobi2_[numgp],2)/pow(jacobi_[numgp],7.0);
      }

      LINALG::TMatrix<FAD,3,1> r_s_FAD;
      LINALG::TMatrix<FAD,3,1> r_ss_FAD;
      for (int i=0;i<3;i++)
      {
        r_s_FAD(i)=r_s(i);
        r_ss_FAD(i)=r_ss(i);
      }

      FAD abs_r_s = Norm(r_s_FAD);

      //calculate Frenet-Serret triad: triad[0]=tangent triad[1]=normal_fs triad[2]=binormal_fs
      //attention: the vectors normal and binormal are set to zero for almost straight (|r'xr''|<1.0e-12) beams
      std::vector<LINALG::TMatrix<FAD,3,1> > triad_fs(3);
      triad_fs=calculate_fs_triad(r_s_FAD, r_ss_FAD);

      //calculate initial curvature
      kappa_vec=calculate_curvature(r_s, r_ss);
      kappa0_[numgp]=kappa_vec.Norm2();

      //Calculate the intermediate triad triad_bar_gp and the corresponding torsion tau_bar_gp at the gauss points
      CalculateIntermediateTriad(r_s_FAD, r_ss_FAD, numgp, triads_ref_nodes, triad_bar_gp, tau_bar_gp);


//      std::cout << "tau_bar_gp: " << tau_bar_gp.val() << std::endl;
//      tau_bar_gp=(taubar0*(xi+1)*(xi+1)/2.0);
//      std::cout << "tau_bar_gp: " << tau_bar_gp.val() << std::endl;

      //calculate the material triad triad_mat_gp at current gp and the mechanical torsion tau_gp out of the reference system
      //attention: for the initial configuration the angle gamma_0 is set to zero!
      CalculateMaterialTriad(xi, numgp, theta_nodes,triad_bar_gp, tau_bar_gp, 0.0, triad_mat_gp,tau_gp);

      //Set initial values of the class variables describing the initial geometry and the initial reference systems
      tau0_[numgp] = tau_gp.val();
      kappa0g20_[numgp]=0.0;
      kappa0g30_[numgp]=0.0;
      for (int i=0;i<3;i++)
      {
        kappa0g20_[numgp]+=kappa_vec(i)*triad_mat_gp(0,i).val();
        kappa0g30_[numgp]+=kappa_vec(i)*triad_mat_gp(1,i).val();
      }

      for (int i=0;i<3;i++)
      {
        t0_[numgp](i)=triad_fs[0](i).val();
        n0_[numgp](i)=triad_mat_gp(0,i).val();
        b0_[numgp](i)=triad_mat_gp(1,i).val();
        dt0ds_[numgp](i)=kappa0_[numgp]*triad_fs[1](i).val();
        db0ds_[numgp](i)=kappa0g20_[numgp]*t0_[numgp](i) - tau0_[numgp]*triad_mat_gp(0,i).val();
      }

    }//End: Loop through all GPs and computation of all relevant values at each gp

    isinit_ = true;
  }//if(!isinit_)


}//DRT::ELEMENTS::Beam3ebanisotrop::SetUpReferenceGeometry()

//Calculates the length of the element based upon a Newton Iteration
//f(l)=l-int(|N'd|)dxi=0
void DRT::ELEMENTS::Beam3ebanisotrop::calculate_length(const std::vector<LINALG::Matrix<3,1> >& xrefe, double tolerance)
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
        r_xi(i)+=xrefe[0](i)*shapefuncderiv(0)+xrefe[1](i)*shapefuncderiv(2)+Tref_[0](i)*shapefuncderiv(1)+Tref_[1](i)*shapefuncderiv(3);
      }
      int_length+=gausspoints.qwgt[numgp]*r_xi.Norm2();

      //derivative of the integral of the length at GP
      for (int i=0; i<3; i++)
      {
        deriv_int+=(Tref_[0](i)*shapefuncderiv(1)/length_+Tref_[1](i)*shapefuncderiv(3)/length_)*r_xi(i);
      }
      deriv_length+=gausspoints.qwgt[numgp]*deriv_int/r_xi.Norm2();
    }
    //cout << endl << "Länge:" << length_ << "\tIntegral:" << int_length << "\tAbleitung:" << deriv_length << endl;
    res=length_-int_length;
    //Update
    length_=length_-res/(1-deriv_length); //the derivative of f(l)=l-int(|N'd|)dxi=0 is f'(l)=1-int(d/dl(|N'd|))dxi
  }

//  //*************************************begin: Determine optimal Hermite constant for circle segment**********************************************************
//
//  if(fabs(Tref_[0].Norm2()-1.0)>1.0e-12 or fabs(Tref_[1].Norm2()-1.0)>1.0e-12)
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
//  n0.Multiply(Sv0,Tref_[0]);
//  n0.Scale(1.0/n0.Norm2());
//  LARGEROTATIONS::computespin(Sn0,n0);
//  scalarproduct10.MultiplyTN(n0,Tref_[1]);
//  if(scalarproduct10.Norm2()>1.0e-12)
//    dserror("Only plane geometries possible at this point!");
//
//  scalarproduct20.MultiplyTN(v0,Tref_[0]);
//  alpha0 = acos(scalarproduct20(0,0));
//  R0=d0/(2*sin(alpha0));
//  R0vec.Multiply(Sn0,Tref_[0]);
//  R0vec.Scale(1/R0vec.Norm2());
//  R0vec.Scale(R0);
//  r0.Update(1.0,xrefe[0],0.0);
//  r0.Update(-1.0,R0vec,1.0);
//
//  for(int i=0;i<3;i++)
//  {
//    triad0(i,0)=R0vec(i);
//    triad0(i,1)=n0(i);
//    triad0(i,2)=Tref_[0](i);
//  }
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
//      b(i)+=Tref_[0](i)*shapefunc(1)+Tref_[1](i)*shapefunc(3);
//      a(i)+=xrefe[0](i)*shapefunc(0)+xrefe[1](i)*shapefunc(2);
//      r_hermite(i)+=xrefe[0](i)*shapefunc(0)+xrefe[1](i)*shapefunc(2)+Tref_[0](i)*length_*shapefunc(1)+Tref_[1](i)*length_*shapefunc(3);
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
//  //length_=100*3.14159265358979323846/(4.0*NUMELE);
//  //std::cout << "length_: " << length_ << std::endl;
//
//  //*************************************end: Determine optimal Hermite constant for circle segment*********************************************************

  return;
}

LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3ebanisotrop::calculate_curvature(const LINALG::Matrix<3,1>& r_s, const LINALG::Matrix<3,1>& r_ss)
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

LINALG::TMatrix<FAD,3,1> DRT::ELEMENTS::Beam3ebanisotrop::calculate_curvature(LINALG::TMatrix<FAD,3,1>& r_s, LINALG::TMatrix<FAD,3,1>& r_ss)
{
  LINALG::TMatrix<FAD,3,1> curvature;

  //spinmatrix Sr' = r'x
  LINALG::TMatrix<FAD,3,3> Srx;
  LARGEROTATIONS::computespin(Srx,r_s);

  //cross-product r'xr''
  LINALG::TMatrix<FAD,3,1> Srxrxx(true);
  Srxrxx.Clear();
  Srxrxx.Multiply(Srx,r_ss);
  FAD abs_r_s = Norm(r_s);

  for (int i=0;i<3;i++)
  {
    curvature(i)=Srxrxx(i)/pow(abs_r_s,2);
  }

  return curvature;
}

double DRT::ELEMENTS::Beam3ebanisotrop::calculate_fstorsion(const LINALG::Matrix<3,1>& r_s, const LINALG::Matrix<3,1>& r_ss, const LINALG::Matrix<3,1>& r_sss, bool straight)
{
  double torsion;

  if(straight==true)
    torsion=0;
  else
  {
    LINALG::Matrix<3,3> rxrxxrxxx;
    for (int i=0; i<3; i++)
    {
      rxrxxrxxx(i,0)=r_s(i);
      rxrxxrxxx(i,1)=r_ss(i);
      rxrxxrxxx(i,2)=r_sss(i);
    }

    //spinmatrix Sr' = r'x
    LINALG::Matrix<3,3> Srx;
    LARGEROTATIONS::computespin(Srx,r_s);

    //cross-product r'xr''
    LINALG::Matrix<3,1> Srxrxx;
    Srxrxx.Clear();
    Srxrxx.Multiply(Srx,r_ss);

    torsion=rxrxxrxxx.Determinant()/Srxrxx.Norm2();
  }

  return torsion;
}
std::vector<double> DRT::ELEMENTS::Beam3ebanisotrop::calculate_jacobi(const LINALG::Matrix<3,1>& r_xi, const LINALG::Matrix<3,1>& r_xixi, const LINALG::Matrix<3,1>& r_xixixi)
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

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebanisotropType::Initialize(DRT::Discretization& dis)
{
  //setting up geometric variables for Beam3ebanisotrop elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a Beam3ebanisotrop element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    //if we get so far current element is a Beam3ebanisotrop element and  we get a pointer at it
    DRT::ELEMENTS::Beam3ebanisotrop* currele = dynamic_cast<DRT::ELEMENTS::Beam3ebanisotrop*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam3ebanisotrop* failed");

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
