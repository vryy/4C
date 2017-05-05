/*----------------------------------------------------------------------*/
/*!
\file red_air_blood_scatraLine3_impl.cpp

\brief Internal implementation of RedAcinus element

\level 3

\maintainer Lena Yoshihara
*/
/*----------------------------------------------------------------------*/



#include "red_air_blood_scatraLine3_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/air_0d_O2_saturation.H"
#include "../drt_mat/hemoglobin_0d_O2_saturation.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface* DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface::Impl(DRT::ELEMENTS::RedAirBloodScatraLine3* red_acinus)
{
  switch (red_acinus->Shape())
  {
  case DRT::Element::line3:
  {
    static RedAirBloodScatraLine3Impl<DRT::Element::line3>* acinus;
    if (acinus==NULL)
    {
      acinus = new RedAirBloodScatraLine3Impl<DRT::Element::line3>;
    }
    return acinus;
  }
  default:
    dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
  }
  return NULL;
}



/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::RedAirBloodScatraLine3Impl()
{

}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Evaluate(
  RedAirBloodScatraLine3*         ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  Teuchos::RCP<MAT::Material> mat)
{

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Initial(
  RedAirBloodScatraLine3*                     ele,
  Teuchos::ParameterList&                         params,
  DRT::Discretization&                   discretization,
  std::vector<int>&                      lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  Teuchos::RCP<Epetra_Vector> generations   = params.get<Teuchos::RCP<Epetra_Vector> >("generations");

  //--------------------------------------------------------------------
  // get the generation numbers
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  {
    int    gid = ele->Id();
    double val = -2.0;
    generations->ReplaceGlobalValues(1,&val,&gid);
  }

}//RedAirBloodScatraLine3Impl::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Sysmat(
  RedAirBloodScatraLine3*                       ele,
  Epetra_SerialDenseVector&                epnp,
  Epetra_SerialDenseVector&                epn,
  Epetra_SerialDenseVector&                epnm,
  Epetra_SerialDenseMatrix&                sysmat,
  Epetra_SerialDenseVector&                rhs,
  Teuchos::RCP<const MAT::Material>        material,
  Teuchos::ParameterList &                          params,
  double                                   time,
  double                                   dt)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::EvaluateTerminalBC(
  RedAirBloodScatraLine3*           ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Epetra_SerialDenseVector&    rhs,
  Teuchos::RCP<MAT::Material>   material)
{

}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::CalcFlowRates(
  RedAirBloodScatraLine3*           ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    elevec1, //a_volumenp,
  Epetra_SerialDenseVector&    elevec2, //a_volume_strain_np,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)

{

}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::GetCoupledValues(
  RedAirBloodScatraLine3*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{

}

/*----------------------------------------------------------------------*
 |  solve the transport of O2 from air to blood             ismail 02/13|
 |                                                                      |
 | Example of use (1):                                                  |
 |--------------------                                                  |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |     |       ^        |              |                                |
 |     |       |        |              |                                |
 |     |(flow direction)|              |                                |
 |     |                |              |                                |
 |     V                V              |                                |
 |     o=====           o===           |                                |
 |  (node1)  =====   (node3)==         |                                |
 |    or 3      ^ ==== or 1   ==       |                                |
 |              |     ======    ===    |                                |
 |              |           =====  ==  V                                |
 |              |                ======o                                |
 |              |                   (node2)                             |
 |              |                                                       |
 |              |                                                       |
 |    [RED_AIR_BLOOD_SCATRA_LINE3 element]                              |
 |                                                                      |
 | Example of use (2):                                                  |
 |--------------------                                                  |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |     |       ^        |              |                                |
 |     |       |        |              |                                |
 |     |(flow direction)|              |                                |
 |     |                |              |                                |
 |     V                V              |                                |
 |     o================o==============o                                |
 |  (node1)     ^    (node2)        (node3)                             |
 |    or 3      |      or 1           or 1                              |
 |              |                                                       |
 |              |                                                       |
 |    [RED_AIR_BLOOD_SCATRA_LINE3 element]                              |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::SolveBloodAirTransport(
  RedAirBloodScatraLine3*      ele,
  Epetra_SerialDenseVector&    dscatra,
  Epetra_SerialDenseVector&    dvo2,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
}
#if 0
{
  const int   myrank  = discretization.Comm().MyPID();

  // get time-step size
  const double dt = params.get<double>("time step size");

  Teuchos::RCP<Epetra_Vector> volnp    = params.get<Teuchos::RCP<Epetra_Vector> >("volumenp");
  Teuchos::RCP<Epetra_Vector> areanp   = params.get<Teuchos::RCP<Epetra_Vector> >("areanp");
  Teuchos::RCP<Epetra_Vector> scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("scatranp");

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> myvolnp(lm.size());
  DRT::UTILS::ExtractMyValues(*volnp,myvolnp,lm);

  // extract local values from the global vectors
  std::vector<double> myareanp(lm.size());
  DRT::UTILS::ExtractMyValues(*areanp,myareanp,lm);

  //--------------------------------------------------------------------
  // define the nodes connected to an acinus
  std::vector<unsigned int> ai;
  std::vector<unsigned int> ci;

  for (unsigned int i=0; i<lm.size();i++)
  {
    if (ele->Nodes()[i]->GetCondition("RedAirwayScatraAirCond")!=NULL)
    {
      ai.push_back(i);
    }
    else if (ele->Nodes()[i]->GetCondition("RedAirwayScatraHemoglobinCond")!=NULL)
    {
      ci.push_back(i);
    }
    else
    {

    }
  }
#if 0
  ai[0] = 2;
  // define the nodes connected to capillaries

  ci[0] = 0;
  ci[1] = 1;
#endif

  // define an empty pO2 vector
  std::vector<double> pO2(lm.size(),0.0);
  // define an empty vO2 vector
  std::vector<double> vO2(lm.size(),0.0);

  // -------------------------------------------------------------------
  // Convert O2 concentration to PO2
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Get O2 properties in air and blood
  // -------------------------------------------------------------------

  // --------------------------------
  // Get O2 properties in air
  // --------------------------------
  int aid = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
  // check if O2 properties material exists
  if (aid==-1)
  {
    dserror("A material defining O2 properties in air could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* amat = DRT::Problem::Instance()->Materials()->ParameterById(aid);
  const MAT::PAR::Air_0d_O2_saturation* aactmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(amat);

  // get atmospheric pressure
  const double patm = aactmat->atmospheric_p_;
  // get number of O2 moles per unit volume of O2
  const double nO2perVO2a = aactmat->nO2_per_VO2_;

  // --------------------------------
  // Get O2 properties in Blood
  // --------------------------------
  int bid = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_hemoglobin_saturation);

  // check if O2 properties material exists
  if (bid==-1)
  {
    dserror("A material defining O2 properties in blood could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* bmat = DRT::Problem::Instance()->Materials()->ParameterById(bid);
  const MAT::PAR::Hemoglobin_0d_O2_saturation* bactmat = static_cast<const MAT::PAR::Hemoglobin_0d_O2_saturation*>(bmat);

  // how much of blood satisfies this rule
  double per_volume_blood = bactmat->per_volume_blood_;
  // The dominating value of E at zero O2 saturation
  double eo_l = bactmat->eo_lower_;
  // Lower decay constant
  double tau_l = bactmat->tau_lower_;
  // Lower VO2 translator
  double vO2o_l = bactmat->vO2o_lower_;
  // Translating E at maximum O2 saturation
  double eo = bactmat->eo_;
  // The dominating value of E at maximum O2 saturation
  double eo_u = bactmat->eo_upper_;
  // Upper decay constant
  double tau_u = bactmat->tau_upper_;
  // Upper VO2 translator
  double vO2o_u = bactmat->vO2o_upper_;
  // Number of O2 moles per VO2
  double nO2perVO2b = bactmat->nO2_per_VO2_;

  // -------------------------------------------------------------------
  // Evaluate VO2 properties in air
  // -------------------------------------------------------------------
  // get acinar volume
  double vAir = myvolnp[ai[0]];
  {
    // calculate the VO2 at nodes
    vO2[ai[0]] = (myscatranp[ai[0]]*vAir)/nO2perVO2a;
    // calculate pO2
    pO2[ai[0]] = vO2[ai[0]]*patm/vAir;
    cout<<"VO2 air: "<<vO2[ai[0]]<<endl;
    cout<<"PO2 air: "<<pO2[ai[0]]<<endl;
    cout<<"Cair: "<<myscatranp[ai[0]]<<endl;
  }

  // -------------------------------------------------------------------
  // Evaluate VO2 in blood
  // -------------------------------------------------------------------
  // get capillary volume
  double vBlood = 0.0;
  for(unsigned int i=0;i<ci.size();i++)
  {
    vBlood += myvolnp[ci[i]];
  }
  vBlood /= float(ci.size());

  // get the ratio of blood volume to the reference saturation volume
  double alpha = vBlood/per_volume_blood;
  {
    // calculate the VO2 at nodes
    for (unsigned int i = 0; i<ci.size(); i++)
    {
      vO2[ci[i]] = (myscatranp[ci[i]]*vBlood)/nO2perVO2b;
    }

    // calculate PO2 at nodes
    for (unsigned int i = 0; i<ci.size(); i++)
    {
      double E  = (eo_l*exp(tau_l*(vO2[ci[i]]/alpha-vO2o_l))
                   + (eo+eo_u*exp(tau_u*(vO2[ci[i]]/alpha-vO2o_u))));
      pO2[ci[i]] = E*vO2[ci[i]]/alpha;
      cout<<"VO2 blood: "<<vO2[ci[i]]<<endl;
      cout<<"PO2 blood: "<<pO2[ci[i]]<<endl;
      cout<<"Vblood: "<<vBlood<<endl;
      cout<<"Cbld: "<<myscatranp[ci[0]]<<endl;
    }
  }

  // -------------------------------------------------------------------
  // Find the deltaVO2 that flowed from air to blood
  // -------------------------------------------------------------------
  const double tol = 1e-7;
  const unsigned int maxItr = 200;
  double error = 1e7;

  // get Diffussion coefficient
  double D = 0.0;
  ele->getParams("DiffusionCoefficient",D);
  // get wall thickness
  double th = 0.0;
  ele->getParams("WallThickness",th);
  // get percentage of diffusion area
  double percDiffArea = 0.0;
  ele->getParams("PercentageOfDiffusionArea",percDiffArea);
  cout<<"Di: "<<D<<";\t th: "<<th<<" \t D\%: "<<percDiffArea<<endl;
  {
    D *= percDiffArea*(myareanp[ai[0]]/th);
  }
  // define vO2 in air
  double vO2a = vO2[ai[0]];
  // define vO2 in capillary
  double vO2b = 0.5*(vO2[ci[0]]+vO2[ci[1]]);
  // define vO2ab which is the volume of blood flowing from air to blood
  double tempE =     eo_l*exp(tau_l*(vO2b/alpha-vO2o_l))
                 + eo+eo_u*exp(tau_u*(vO2b/alpha-vO2o_u));
  double tempPO2b= tempE*vO2b/alpha;
  double vO2ab = 0.0;//D*(pO2[ai[0]] - tempPO2b)*dt;
  // loop till convergence
  for (unsigned int itr = 0; fabs(error)>tol;itr++)
  {
    // -----------------------------------------------------------------
    // Evaluate all terms associated with flow rate of O2
    // -----------------------------------------------------------------
    // flow rate of O2 from air to blood
    double qO2  = vO2ab/dt;
    // derivative of O2 flowrate w.r.t volume of O2 flowing into blood
    double dqO2_dvO2ab = 1.0/dt;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in air
    // -----------------------------------------------------------------
    // new volume of O2 in air is old volume minus the one flowing into
    // blood
    double vO2ap= vO2a - vO2ab;
    // pO2 in air
    double pO2a = vO2ap*patm/vAir;
    // derivative of pO2a w.r.t volume of O2 flowing into blood
    double dpO2a_dvO2 = -patm/vAir;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in blood
    // -----------------------------------------------------------------
    // volumen of O2 in blood is old volume plus the one flowing inot
    // blood
    double vO2bp= vO2b + vO2ab;
    // pO2 in blood: eb*vO2ab/alpha
    double eb  =     (eo_l*exp(tau_l*(vO2bp/alpha-vO2o_l))
                 + eo+eo_u*exp(tau_u*(vO2bp/alpha-vO2o_u)));
    double pO2b = eb*vO2bp/alpha;
    // derivative of pO2b w.r.t volume of O2 flowing into blood
    double deb_dvO2ab =  (eo_l*exp(tau_l*(vO2bp/alpha-vO2o_l)))*(tau_l/alpha)
                        +(eo_u*exp(tau_u*(vO2bp/alpha-vO2o_u)))*(tau_u/alpha);
    double dpO2b_dvO2ab= (deb_dvO2ab*vO2bp + eb)/alpha;

    // -----------------------------------------------------------------
    // perform the Newton-Raphson step
    // -----------------------------------------------------------------
    // calculate residual
    double f = qO2 - D*(pO2a - pO2b);
    // calculate df/dvO2ab
    double df_dvO2ab = dqO2_dvO2ab - D*(dpO2a_dvO2 - dpO2b_dvO2ab);
    // calculate the corrector
    double dvO2ab = f/(df_dvO2ab);
    // corrector vO2ab
    vO2ab -= dvO2ab;

    // evaluate a dimensionless error
    error = fabs(f*dt);

    // check append iteration step
    if (itr>maxItr)
    {
      dserror("[Warning in ELEMENT(%d)]: solving for VO2(air/blood) did not converge Error: %f\n",ele->Id(),fabs(error));
      exit(1);
    }
  }
  cout<<"+-------------------->>>>>>>>>>>>>>--------------------+"<<endl;
  cout<<"| Diffc: "<<D<<"\t DiffArea: "<<myareanp[ai[0]]<<endl;
  cout<<"| Vo2_a: "<<vO2a<<endl;
  cout<<"| Vo2_b: "<<vO2b<<endl;
  cout<<"| dVO2 : "<<vO2ab<<endl;
  cout<<"+-------------------->>>>>>>>>>>>>>--------------------+"<<endl;
  // -----------------------------------------------------------------
  // update the change of vO2 into dvO2 vector:
  // -----------------------------------------------------------------
  // update dscatra in air
  dscatra[ai[0]] =-vO2ab*nO2perVO2a/vAir;

  // update dscatra in blood
  dscatra[ci[0]] = vO2ab*nO2perVO2b/vBlood;

}
#endif
#if 0
{
  const int   myrank  = discretization.Comm().MyPID();

  // get time-step size
  const double dt = params.get<double>("time step size");

  Teuchos::RCP<Epetra_Vector> volnp    = params.get<Teuchos::RCP<Epetra_Vector> >("volumenp");
  Teuchos::RCP<Epetra_Vector> areanp   = params.get<Teuchos::RCP<Epetra_Vector> >("areanp");
  Teuchos::RCP<Epetra_Vector> scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("scatranp");

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> myvolnp(lm.size());
  DRT::UTILS::ExtractMyValues(*volnp,myvolnp,lm);

  // extract local values from the global vectors
  std::vector<double> myareanp(lm.size());
  DRT::UTILS::ExtractMyValues(*areanp,myareanp,lm);


  //--------------------------------------------------------------------
  // define the nodes connected to an acinus
  std::vector<unsigned int> ai(1,0);
  ai[0] = 2;
  // define the nodes connected to capillaries
  std::vector<unsigned int> ci(lm.size()-1,0);
  ci[0] = 0;
  ci[1] = 1;

  // define an empty pO2 vector
  std::vector<double> pO2(lm.size(),0.0);
  // define an empty vO2 vector
  std::vector<double> vO2(lm.size(),0.0);

  // -------------------------------------------------------------------
  // Convert O2 concentration to PO2
  // -------------------------------------------------------------------
  // get node coordinates and number of elements per node
  const int numnode = iel;
  DRT::Node** nodes = ele->Nodes();
  // get airway length
  LINALG::Matrix<3,iel> xyze;
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // Calculate the length of airway element
  const double length=sqrt(
      pow(xyze(0,ci[0]) - xyze(0,ci[1]),2)
    + pow(xyze(1,ci[0]) - xyze(1,ci[1]),2)
    + pow(xyze(2,ci[0]) - xyze(2,ci[1]),2));

  // print out node coorsd
  cout<<"air node ("<<xyze(0,ai[0])<<","<<xyze(1,ai[0])<<","<<xyze(2,ai[0])<<")"<<endl;
  cout<<"bld node1("<<xyze(0,ci[0])<<","<<xyze(1,ci[0])<<","<<xyze(2,ci[0])<<")"<<endl;
  cout<<"bld node2("<<xyze(0,ci[1])<<","<<xyze(1,ci[1])<<","<<xyze(2,ci[1])<<")"<<endl;

  // -------------------------------------------------------------------
  // Get O2 properties in air and blood
  // -------------------------------------------------------------------

  // --------------------------------
  // Get O2 properties in air
  // --------------------------------
  int aid = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
  // check if O2 properties material exists
  if (aid==-1)
  {
    dserror("A material defining O2 properties in air could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* amat = DRT::Problem::Instance()->Materials()->ParameterById(aid);
  const MAT::PAR::Air_0d_O2_saturation* aactmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(amat);

  // get atmospheric pressure
  const double patm = aactmat->atmospheric_p_;
  // get number of O2 moles per unit volume of O2
  const double nO2perVO2a = aactmat->nO2_per_VO2_;

  // --------------------------------
  // Get O2 properties in Blood
  // --------------------------------
  int bid = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_hemoglobin_saturation);
  cout<<"O2 properties in air: "<<"Patm("<<aactmat->atmospheric_p_<<"),"<<"NO2perVO2("<<aactmat->nO2_per_VO2_<<")"<<endl;

  // check if O2 properties material exists
  if (bid==-1)
  {
    dserror("A material defining O2 properties in blood could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* bmat = DRT::Problem::Instance()->Materials()->ParameterById(bid);
  const MAT::PAR::Hemoglobin_0d_O2_saturation* bactmat = static_cast<const MAT::PAR::Hemoglobin_0d_O2_saturation*>(bmat);
  cout<<"O2 properties in blood: "<<
    "per_volume_blood_("<<bactmat->per_volume_blood_<<")"<<
    "eo_lower("<<bactmat->eo_lower_<<")"<<
    "tau_lower("<<bactmat->tau_lower_<<")"<<
    "vO2o_lower("<<bactmat->vO2o_lower_<<")"<<
    "eo("<<bactmat->eo_<<")"<<
    "eo_upper("<<bactmat->eo_upper_<<")"<<
    "tau_upper("<<bactmat->tau_upper_<<")"<<
    "vO2_upper("<<bactmat->vO2o_upper_<<")"<<
    "nO2perVO2("<<bactmat->nO2_per_VO2_<<")"<<endl;

  // how much of blood satisfies this rule
  double per_volume_blood = bactmat->per_volume_blood_;
  // The dominating value of E at zero O2 saturation
  double eo_l = bactmat->eo_lower_;
  // Lower decay constant
  double tau_l = bactmat->tau_lower_;
  // Lower VO2 translator
  double vO2o_l = bactmat->vO2o_lower_;
  // Translating E at maximum O2 saturation
  double eo = bactmat->eo_;
  // The dominating value of E at maximum O2 saturation
  double eo_u = bactmat->eo_upper_;
  // Upper decay constant
  double tau_u = bactmat->tau_upper_;
  // Upper VO2 translator
  double vO2o_u = bactmat->vO2o_upper_;
  // Number of O2 moles per VO2
  double nO2perVO2b = bactmat->nO2_per_VO2_;

  // -------------------------------------------------------------------
  // Evaluate VO2 properties in air
  // -------------------------------------------------------------------
  // get acinar volume
  double vAir = myvolnp[ai[0]];
  cout<<"volume of air:"<<vAir<<endl;
  {
    // calculate the VO2 at nodes
    vO2[ai[0]] = (myscatranp[ai[0]]*vAir)/nO2perVO2a;
    // calculate pO2
    pO2[ai[0]] = vO2[ai[0]]*patm/vAir;
  }
  cout<<"VO2a: "<<vO2[ai[0]]<<endl;
  cout<<"PO2a: "<<pO2[ai[0]]<<endl;

  // -------------------------------------------------------------------
  // Evaluate VO2 in blood
  // -------------------------------------------------------------------
  // get capillary volume
  cout<<"length of capillaries: "<<length<<endl;
  double vBlood = 0.5*(myareanp[ci[0]]+myareanp[ci[1]])*length;
  cout<<"Area of capillaries: ("<<myareanp[ci[0]]<<","<<myareanp[ci[1]]<<")"<<endl;
  cout<<"length of capillaries:"<<length<<endl;
  cout<<"vBlood: "<<vBlood<<endl;
  // get the ratio of blood volume to the reference saturation volume
  double alpha = vBlood/per_volume_blood;
  cout<<"alpha: "<<alpha<<endl;
  {
    // calculate the VO2 at nodes
    for (unsigned int i = 0; i<ci.size(); i++)
    {
      vO2[ci[i]] = (myscatranp[ci[i]]*vBlood)/nO2perVO2b;
    }

    // calculate PO2 at nodes
    for (unsigned int i = 0; i<ci.size(); i++)
    {
      double E  = (eo_l*exp(tau_l*(vO2[ci[i]]/alpha-vO2o_l))
                   + (eo+eo_u*exp(tau_u*(vO2[ci[i]]/alpha-vO2o_u))));
      pO2[ci[i]] = E*vO2[ci[i]]/alpha;

      cout<<"VO2b_"<<ci[i]<<": "<<vO2[ci[i]]<<endl;
      cout<<"PO2b_"<<ci[i]<<": "<<pO2[ci[i]]<<endl;
    }

  }

  // -------------------------------------------------------------------
  // Find the deltaVO2 that flowed from air to blood
  // -------------------------------------------------------------------
  const double tol = 1e-7;
  const unsigned int maxItr = 200;
  double error = 1e7;

  // get Diffussion coefficient
  double D = 0.0;
  ele->getParams("DiffusionCoefficient",D);
  // get wall thickness
  double th = 0.0;
  ele->getParams("WallThickness",th);
  // get percentage of diffusion area
  double percDiffArea = 0.0;
  ele->getParams("PercentageOfDiffusionArea",percDiffArea);
  {
    D *= percDiffArea*(myareanp[ai[0]]/th);
  }
  // define vO2 in air
  double vO2a = vO2[ai[0]];
  // define vO2 in capillary
  double vO2b = 0.5*(vO2[ci[0]]+vO2[ci[1]]);
  // define vO2ab which is the volume of blood flowing from air to blood
  double tempE =     eo_l*exp(tau_l*(vO2b/alpha-vO2o_l))
                 + eo+eo_u*exp(tau_u*(vO2b/alpha-vO2o_u));
  double tempPO2b= tempE*vO2b/alpha;
  double vO2ab = 0.0;//D*(pO2[ai[0]] - tempPO2b)*dt;
  // loop till convergence
  for (unsigned int itr = 0; fabs(error)>tol;itr++)
  {
    cout<<"+------------------ITR"<<itr<<"-----------------+"<<endl;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with flow rate of O2
    // -----------------------------------------------------------------
    // flow rate of O2 from air to blood
    double qO2  = vO2ab/dt;
    // derivative of O2 flowrate w.r.t volume of O2 flowing into blood
    double dqO2_dvO2ab = 1.0/dt;
    cout<<"Qo2: "<<qO2<<endl;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in air
    // -----------------------------------------------------------------
    // new volume of O2 in air is old volume minus the one flowing into
    // blood
    double vO2ap= vO2a - vO2ab;
    // pO2 in air
    double pO2a = vO2ap*patm/vAir;
    // derivative of pO2a w.r.t volume of O2 flowing into blood
    double dpO2a_dvO2 = -patm/vAir;
    cout<<"vO2a: "<<vO2ap<<endl;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in blood
    // -----------------------------------------------------------------
    // volumen of O2 in blood is old volume plus the one flowing inot
    // blood
    double vO2bp= vO2b + vO2ab;
    // pO2 in blood: eb*vO2ab/alpha
    double eb  =     (eo_l*exp(tau_l*(vO2bp/alpha-vO2o_l))
                 + eo+eo_u*exp(tau_u*(vO2bp/alpha-vO2o_u)));
    double pO2b = eb*vO2bp/alpha;
    // derivative of pO2b w.r.t volume of O2 flowing into blood
    double deb_dvO2ab =  (eo_l*exp(tau_l*(vO2bp/alpha-vO2o_l)))*(tau_l/alpha)
                        +(eo_u*exp(tau_u*(vO2bp/alpha-vO2o_u)))*(tau_u/alpha);
    double dpO2b_dvO2ab= (deb_dvO2ab*vO2bp + eb)/alpha;
    cout<<"vO2b: "<<vO2bp<<endl;

    cout<<"dE/dvO2ab: "<<deb_dvO2ab<<endl;
    cout<<"dPO2/dvO2ab: "<<dpO2b_dvO2ab<<endl;

    // -----------------------------------------------------------------
    // perform the Newton-Raphson step
    // -----------------------------------------------------------------
    // calculate residual
    double f = qO2 - D*(pO2a - pO2b);
    // calculate df/dvO2ab
    double df_dvO2ab = dqO2_dvO2ab - D*(dpO2a_dvO2 - dpO2b_dvO2ab);
    // calculate the corrector
    double dvO2ab = f/(df_dvO2ab);
    // corrector vO2ab
    vO2ab -= dvO2ab;
    cout<<"vO2ab: "<<vO2ab<<endl;
    cout<<"f: "<<f<<endl;
    cout<<"df/dvO2ab: "<<df_dvO2ab<<endl;
    cout<<"del_VO2: "<<dvO2ab<<endl;
    // evaluate a dimensionless error
    error = fabs(f*dt);

    cout<<"+-----------------------------------------------+"<<endl;
    // check append iteration step
    if (itr>maxItr)
    {
      dserror("[Warning in ELEMENT(%d)]: solving for VO2(air/blood) did not converge Error: %f\n",ele->Id(),fabs(error));
      exit(1);
//      cout<<"[Warning in ELEMENT("<<ele->Id()<<")]: solving for VO2(air/blood) did not converge (Error: "<<fabs(error)<<")"<<endl;
    }
  }
  cout<<"Change in O2 is: "<<vO2ab<<endl;

  // -----------------------------------------------------------------
  // update the change of vO2 into dvO2 vector:
  //  In this step we need to take care, that in blood, a small change
  //  in VO2 can result in an unrealistic PO2 values. This is due to
  //  the exponential dominating relationship between PO2 and VO2
  //  near the fully saturated Hemoglobin.
  //
  //  To resolve this issue we do the following:
  //  1 - Update VO2 in air and blood.
  //  2 - Evaluate the PO2 in air and blood.
  //  3 - if PO2b upstream > PO2a then PO2b upstream = PO2a
  //  4 - if (3) is true add the remaining vO2ab to vO2b downstream
  // -----------------------------------------------------------------

  // ---------------------------------------
  // (step 1) define the new vO2 vector
  // ---------------------------------------
  std::vector<double> nvO2(vO2);
  // update the volume of O2 in air
  for (unsigned int i = 0; i<ai.size();i++)
  {
    nvO2[ai[i]] = vO2[ai[i]]- vO2ab;
  }
  // update the volume of O2 in blood
  for (unsigned int i = 0; i<ci.size();i++)
  {
    nvO2[ci[i]] = vO2[ci[i]]+ vO2ab;
  }

  // ---------------------------------------
  // (step 2-1) evaluate the new pO2
  // ---------------------------------------
  std::vector<double> npO2(lm.size(),0.0);
  // Evaluate new pO2 in air
  npO2[ai[0]] = nvO2[ai[0]]*patm/vAir;
  cout<<"New PO2 air: "<<npO2[ai[0]]<<endl;
  // ---------------------------------------
  // (step 2-2) Evaluate PO2 in blood
  // ---------------------------------------
  for (unsigned int i = 0; i<ci.size(); i++)
  {
    double E  = (eo_l*exp(tau_l*(nvO2[ci[i]]/alpha-vO2o_l))
                 + (eo+eo_u*exp(tau_u*(nvO2[ci[i]]/alpha-vO2o_u))));
    npO2[ci[i]] = E*nvO2[ci[i]]/alpha;
    cout<<"New PO2 bld"<<i<<": "<<npO2[ci[i]]<<endl;
  }


  // ---------------------------------------
  // (step 3) check if PO2b upstream > PO2a
  //          if true: correct PO2b upstream
  // ---------------------------------------

  // check which node is upstream and which is downstream
  std::vector<unsigned int> nci(ci);
  if (npO2[ci[0]] > npO2[ci[1]])
  {
    cout<<"True!!!! Flipping..."<<endl;
    // invert the numbering of the nodes
    for (unsigned int i = 0; i<ci.size(); i++)
    {
      nci[ci.size()-1-i] = ci[i];
    }
  }

  // check if PO2b upstream > PO2a
  // if true: set PO2b upstream = PO2a
  if (npO2[nci[1]]>npO2[ai[0]])
    npO2[nci[1]]=npO2[ai[0]];

  // ---------------------------------------
  // (step 4) if (step 3) is true add the
  // remaining vO2ab to vO2b downstream
  // ---------------------------------------

  // find the new vO2b upstream
  double vO2bu =  vO2[nci[1]];
  double pO2bu = npO2[nci[1]];
  error = 1e7;
  cout<<"+-----------------------------------------------+"<<endl;
  for (unsigned int itr = 0; fabs(error)>tol;itr++)
  {
    // evaluate Eb upstream
    double eb=  (eo_l*exp(tau_l*(vO2bu/alpha-vO2o_l))
              + (eo+eo_u*exp(tau_u*(vO2bu/alpha-vO2o_u))));

    // derivative of eb w.r.t upstream volume of O2 in blood
    double deb_dvO2bu =  (eo_l*exp(tau_l*(vO2bu/alpha-vO2o_l)))*(tau_l/alpha)
                        +(eo_u*exp(tau_u*(vO2bu/alpha-vO2o_u)))*(tau_u/alpha);

    // --------------------------------
    // perform the Newton-Raphson step
    // --------------------------------
    // calculate residual
    double f = pO2bu - eb*vO2bu/alpha;
    // calculate the residual derivative w.r.t vO2bu
    double df_dvO2bu =  -(deb_dvO2bu*vO2bu/alpha  + eb/alpha);
    // calculate the corrector
    double dvO2bu = f/(df_dvO2bu);
    // corrector vO2ab
    vO2bu -= dvO2bu;
    // evaluate a dimensionless error
    // error = f/vBlood;
    error = (pO2bu/eb*alpha - vO2bu)/vBlood;

    // check append iteration step

    if (itr>maxItr)
    {
      dserror("[Warning in ELEMENT(%d)]: solving for VO2(blood) did not converge Error: %f\n",ele->Id(),fabs(error));
      exit(1);
    }
  }
  cout<<"+-----------------------------------------------+"<<endl;
  // update vO2 upstream
  nvO2[nci[1]] = vO2bu;

  // correct vO2b downstream
  nvO2[nci[0]] +=  vO2ab - (nvO2[nci[1]] - vO2[nci[1]]);


  for (unsigned int i = 0; i<ci.size(); i++)
  {
    double E  = (eo_l*exp(tau_l*(nvO2[ci[i]]/alpha-vO2o_l))
                 + (eo+eo_u*exp(tau_u*(nvO2[ci[i]]/alpha-vO2o_u))));
    cout<<"Corrected PO2 bld"<<i<<": "<<E*nvO2[ci[i]]/alpha<<endl;
    cout<<"Corrected VO2 bld"<<i<<": "<<nvO2[ci[i]]<<endl;
  }


  // -------------------------------------------------------------------
  // convert the new vO2 in blood and air into a concentration
  // and update the change of concentration into dscatra
  // -------------------------------------------------------------------
  // update the concentration of O2 in air
  dscatra[ai[0]] = (-vO2ab)*nO2perVO2a/vAir;

  cout<<"dVO2 (air/bld): "<<vO2ab<<endl;
  cout<<"dVO2 air "<<nvO2[ai[0]]-vO2[ai[0]]<<endl;
  // update the concentration of O2 in blood
  for (unsigned int i = 0; i<ci.size(); i++)
  {
    dscatra[ci[i]]          =(nvO2[ci[i]]-vO2[ci[i]])*nO2perVO2b/vBlood;
    cout<<"dVO2 bld"<<i<<" "<<nvO2[ci[i]]-vO2[ci[i]]<<endl;
    cout<<"new Scatra"<<i<<" "<<myscatranp[ci[i]]+dscatra[ci[i]]<<endl;
  }
}
#endif
