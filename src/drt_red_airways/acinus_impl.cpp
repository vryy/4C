/*----------------------------------------------------------------------*/
/*!
\file acinus_impl.cpp

\brief Internal implementation of RedAcinus element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/




#include "acinus_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/maxwell_0d_acinus.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinusImplInterface* DRT::ELEMENTS::RedAcinusImplInterface::Impl(DRT::ELEMENTS::RedAcinus* red_acinus)
{
  switch (red_acinus->Shape())
  {
  case DRT::Element::line2:
  {
    static AcinusImpl<DRT::Element::line2>* acinus;
    if (acinus==NULL)
    {
      acinus = new AcinusImpl<DRT::Element::line2>;
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
DRT::ELEMENTS::AcinusImpl<distype>::AcinusImpl()
{

}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcinusImpl<distype>::Evaluate(
  RedAcinus*                 ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RefCountPtr<MAT::Material> mat)
{
  //  const int   myrank  = discretization.Comm().MyPID();
  
  //  const int numnode = iel;
  const int elemVecdim = elevec1_epetra.Length () ;
  vector<int>::iterator it_vcr;

  // construct views
  //  LINALG::Matrix<1*iel,1*iel> elemat1(elemat1_epetra.A(),true);
  //  LINALG::Matrix<1*iel,    1> elevec1(elevec1_epetra.A(),true);
  // elemat2, elevec2, and elevec3 are never used anyway

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get time-step size
  const double dt = params.get<double>("time step size");

  // get time
  const double time = params.get<double>("total time");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------


  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow, pressure,
  // ---------------------------------------------------------------------

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  RefCountPtr<const Epetra_Vector> pn   = discretization.GetState("pn");
  RefCountPtr<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  RefCountPtr<Epetra_Vector> acinar_vnp  = params.get<RefCountPtr<Epetra_Vector> >("acinar_vnp");
  RefCountPtr<Epetra_Vector> acinar_vn   = params.get<RefCountPtr<Epetra_Vector> >("acinar_vn");


  RefCountPtr<Epetra_Vector> qin_nm  = params.get<RefCountPtr<Epetra_Vector> >("qin_nm");
  RefCountPtr<Epetra_Vector> qin_n   = params.get<RefCountPtr<Epetra_Vector> >("qin_n");
  RefCountPtr<Epetra_Vector> qin_np  = params.get<RefCountPtr<Epetra_Vector> >("qin_np");

  RefCountPtr<Epetra_Vector> qout_np = params.get<RefCountPtr<Epetra_Vector> >("qout_np");
  RefCountPtr<Epetra_Vector> qout_n  = params.get<RefCountPtr<Epetra_Vector> >("qout_n");
  RefCountPtr<Epetra_Vector> qout_nm = params.get<RefCountPtr<Epetra_Vector> >("qout_nm");
  
  RefCountPtr<Epetra_Vector> sysmat_iad = params.get<RefCountPtr<Epetra_Vector> >("sysmat_iad");

 if (pnp==null || pn==null || pnm==null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // extract local values from the global vectors
  vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  // extract local values from the global vectors
  vector<double> mypnm(lm.size());
  DRT::UTILS::ExtractMyValues(*pnm,mypnm,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(elemVecdim);
  Epetra_SerialDenseVector epn (elemVecdim);
  Epetra_SerialDenseVector epnm(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)   = mypnp[i];
    epn(i)    = mypn[i];
    epnm(i)   = mypnm[i];
  }

  double e_acin_e_vnp;
  double e_acin_e_vn;

  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    e_acin_e_vnp = (*acinar_vnp)[ele->LID()];
    e_acin_e_vn  = (*acinar_vn )[ele->LID()];
  }

  // get the volumetric flow rate from the previous time step
  ParameterList elem_params;
  elem_params.set<double>("qout_np",(*qout_np)[ele->LID()]);
  elem_params.set<double>("qout_n" ,(*qout_n )[ele->LID()]);
  elem_params.set<double>("qout_nm",(*qout_nm)[ele->LID()]);
  elem_params.set<double>("qin_np" ,(*qin_np )[ele->LID()]);
  elem_params.set<double>("qin_n"  ,(*qin_n  )[ele->LID()]);
  elem_params.set<double>("qin_nm" ,(*qin_nm )[ele->LID()]);

  elem_params.set<double>("acin_vnp" ,e_acin_e_vnp);
  elem_params.set<double>("acin_vn"  ,e_acin_e_vn );

  elem_params.set<double>("lungVolume_np",params.get<double>("lungVolume_np"));
  elem_params.set<double>("lungVolume_n",params.get<double>("lungVolume_n"));
  elem_params.set<double>("lungVolume_nm",params.get<double>("lungVolume_nm"));
  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         epnp,
         epn,
         epnm,
         elemat1_epetra,
         elevec1_epetra,
         mat,
         elem_params,
         time,
         dt);

  for (unsigned int i=0;i<lm.size();i++)
  {
    int    gid = lm[i];
    double val = elemat1_epetra(i,i);
    sysmat_iad->ReplaceGlobalValues(1,&val,&gid);
  }

#if 0
  cout<<">>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"sysmat: "<<elemat1_epetra<<endl;
  cout<<">>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"rhs: "<<elevec1_epetra<<endl;
  cout<<">>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
#endif

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::Initial(
  RedAcinus*                             ele,
  ParameterList&                         params,
  DRT::Discretization&                   discretization,
  vector<int>&                           lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  const int   myrank  = discretization.Comm().MyPID();

  RCP<Epetra_Vector> p0np    = params.get<RCP<Epetra_Vector> >("p0np");
  RCP<Epetra_Vector> p0n     = params.get<RCP<Epetra_Vector> >("p0n");
  RCP<Epetra_Vector> p0nm    = params.get<RCP<Epetra_Vector> >("p0nm");

  RCP<Epetra_Vector> radii   = params.get<RCP<Epetra_Vector> >("radii");

  RCP<Epetra_Vector> generations   = params.get<RCP<Epetra_Vector> >("generations");
  RCP<Epetra_Vector> a_bc          = params.get<RCP<Epetra_Vector> >("acini_bc");

  //  RCP<Epetra_Vector> a_volume      = params.get<RCP<Epetra_Vector> >("acini_volume");
  RCP<Epetra_Vector> a_e_volume    = params.get<RCP<Epetra_Vector> >("acini_e_volume");

  //  vector<int>::iterator it = lm.begin();

  //vector<int> lmowner;
  vector<int> lmstride;
  RCP<vector<int> > lmowner = Teuchos::rcp(new vector<int>);
  ele->LocationVector(discretization,lm,*lmowner,lmstride);

  //--------------------------------------------------------------------
  // Initialize the pressure vectors
  //--------------------------------------------------------------------
  if(myrank == (*lmowner)[0])
  {
    int    gid = lm[0];
    double val = 0.0;
    p0np->ReplaceGlobalValues(1,&val,&gid);
    p0n ->ReplaceGlobalValues(1,&val,&gid);
    p0nm->ReplaceGlobalValues(1,&val,&gid);

  }
  {
    
    // find the volume of an acinus condition
    int    gid2 = ele->Id();
    double acin_vol = 0.0;
    ele->getParams("AcinusVolume",acin_vol);
    a_e_volume->ReplaceGlobalValues(1,&acin_vol,&gid2);
  }

  //--------------------------------------------------------------------
  // get the generation numbers
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  for (int i = 0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("RedLungAcinusCond")||ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
    {
      // find the acinus condition
      int    gid = ele->Id();
      double val = 1.0;
      a_bc->ReplaceGlobalValues(1,&val,&gid);
    }
  }
  {
    int    gid = ele->Id();
    int    generation = -1;
    double val = double(generation);
    generations->ReplaceGlobalValues(1,&val,&gid);
  }

#if 0
  cout<<"!!!!!!!!!!!!-----------------------!!!!!!!!!!!!"<<endl;
  cout<<"!!!!!!!!!!!!------ BRILLIANT ------!!!!!!!!!!!!"<<endl;
  cout<<"!!!!!!!!!!!!-----------------------!!!!!!!!!!!!"<<endl;
#endif

}//AcinusImpl::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::Sysmat(
  RedAcinus*                               ele,
  Epetra_SerialDenseVector&                epnp,
  Epetra_SerialDenseVector&                epn,
  Epetra_SerialDenseVector&                epnm,
  Epetra_SerialDenseMatrix&                sysmat,
  Epetra_SerialDenseVector&                rhs,
  Teuchos::RCP<const MAT::Material>        material,
  ParameterList &                          params,
  double                                   time,
  double                                   dt)
{
  //  const int elemVecdim = epnp.Length () ;
  double E1 = 0.0;
  double E2 = 0.0;
  double Rt = 0.0;
  double Ra = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_0d_maxwell_acinus)
  {
    // get actual material
    const MAT::Maxwell_0d_acinus* actmat = static_cast<const MAT::Maxwell_0d_acinus*>(material.get());
    E1 = actmat->Stiffness1();
    E2 = actmat->Stiffness2();
    Rt = actmat->Viscosity1();
    Ra = actmat->Viscosity2();
#if 0
    cout<<"E1: "<<E1<<endl;
    cout<<"E2: "<<E2<<endl;
    cout<<"Rt: "<<Rt<<endl;
    cout<<"Ra: "<<Ra<<endl;
#endif
  }
  else
  {
    dserror("Material law is not a Newtonia fluid");
    exit(1);
  }

#if 0
  // set element data
  const int numnode = iel;
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();

  LINALG::Matrix<3,iel> xyze;
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

  // Calculate the length of acinus element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));
#endif


  // check here, if we really have an acinus !!


  rhs.Scale(0.0);
  sysmat.Scale(0.0);


  double q_out    = params.get<double>("qout_n");
  //  double qout_np  = params.get<double>("qout_np");

  //  double q_in    = params.get<double>("qin_n");
  //  double qin_np  = params.get<double>("qin_np");

  // get the generation number
  int generation = 0;
  ele->getParams("Generation",generation);


  /********** START HERE ************/

  // -------------------------------------------------------------------
  // Adding the acinus model if Prescribed
  // -------------------------------------------------------------------
  
  double acin_vnp = params.get<double>("acin_vnp");
  double acin_vn  = params.get<double>("acin_vn");
  
  
  //------------------------------------------------------------------
  // Evaluate the number of acini on the end of an outlet.
  // This is hard coded for now, but will be fixed later
  //
  // For now Schroters model is assumed to be Vmodel = 1mm^3
  // Thus we should be given the total-acini-volume/total-t-bronchi-area
  // to calculate the number of Schroters acini at each outlet
  //------------------------------------------------------------------
  
  double qnp = params.get<double>("qin_np");
  double qn  = params.get<double>("qin_n");
  double qnm = params.get<double>("qin_nm");
 
  
  //----------------------------------------------------------------
  // Read in the material information
  //----------------------------------------------------------------
  double VolAcinus;
  ele->getParams("AcinusVolume",VolAcinus);
  double volAlvDuct;
  ele->getParams("AlveolarDuctVolume",volAlvDuct);
  //  cout<<"Acinus Vol: "<<VolAcinus<<endl;
  //  cout<<"AlvDuc Vol: "<<volAlvDuct<<endl;
  const double NumOfAcini = double(floor(VolAcinus/volAlvDuct));
  
  if (NumOfAcini < 1.0)
  {
    dserror("Acinus condition at node (%d) has zero acini",ele->Id());
  }
  
  
  double p1nm = epnm(0);
  double p1n  = epn(0);

  double p2nm = epnm(1);
  double p2n  = epn(1);
  
  double vnp= acin_vnp;
  double vn = acin_vn;  

  if (ele->Type() ==  "NeoHookean")
  {
    const double Kp_np = 1.0/(E1*dt);
    const double Kp_n  = 1.0/(E1*dt);
    
    sysmat(0,0) = -1.0*(Kp_np)*NumOfAcini; sysmat(0,1) =  1.0*(Kp_np)*NumOfAcini;
    sysmat(1,0) =  1.0*(Kp_np)*NumOfAcini; sysmat(1,1) = -1.0*(Kp_np)*NumOfAcini;
    
    rhs(0)      = -1.0*(Kp_n*(p1n-p2n))*NumOfAcini;
    rhs(1)      =  1.0*(Kp_n*(p1n-p2n))*NumOfAcini;
  }
  else if (ele->Type() ==  "KelvinVoigt")
  {
    const double Kp_np = 2.0/(E1*dt) + 1.0/Rt;
    const double Kp_n  = 2.0/(E1*dt) - 1.0/Rt;
    
    sysmat(0,0)  = -1.0*(Kp_np*NumOfAcini); sysmat(0,1)  =  1.0*(Kp_np*NumOfAcini);
    sysmat(1,0)  =  1.0*(Kp_np*NumOfAcini); sysmat(1,1)  = -1.0*(Kp_np*NumOfAcini);
    
    rhs(0)       = -1.0*(q_out + p1n*NumOfAcini * Kp_n);
    rhs(1)       =  1.0*(q_out + p1n*NumOfAcini * Kp_n);
  }
  else if (ele->Type() ==  "ViscoElastic_2dof")
  {
    const double Kp_np = E2/dt + Rt/(dt*dt);
    const double Kp_n  = E2/dt + 2.0*Rt/(dt*dt);
    const double Kp_nm = -Rt/(dt*dt);
    const double Kq_np = E1*E2 + Rt*(E1 + E2)/dt;
    const double Kq_n  = -Rt*(E1+E2)/dt;
    
    sysmat(0,0) = -1.0*( Kp_np/Kq_np)*NumOfAcini; sysmat(0,1) =  1.0*( Kp_np/Kq_np)*NumOfAcini;
    sysmat(1,0) =  1.0*( Kp_np/Kq_np)*NumOfAcini; sysmat(1,1) = -1.0*( Kp_np/Kq_np)*NumOfAcini;
    
    rhs(0)      = -1.0*((Kp_n*(p1n-p2n) + Kp_nm*(p1nm-p2nm))*NumOfAcini/Kq_np + Kq_n/Kq_np*qn);
    rhs(1)      =  1.0*((Kp_n*(p1n-p2n) + Kp_nm*(p1nm-p2nm))*NumOfAcini/Kq_np + Kq_n/Kq_np*qn);
    
  }
  else if (ele->Type() ==  "ViscoElastic_3dof")
  {
    const double Kp_np = Rt/(dt*dt) + E2/dt;
    const double Kp_n  = -2.0*Rt/(dt*dt) - E2/dt;
    const double Kp_nm = Rt/(dt*dt);
    
    const double Kq_np = Ra*Rt/(dt*dt) + (E1*Rt + E2*Ra + E2*Rt)/dt + E1*E2;
    const double Kq_n  = -2.0*Ra*Rt/(dt*dt) - (E1*Rt + E2*Ra + E2*Rt)/dt;
    const double Kq_nm = Ra*Rt/(dt*dt);
    
    sysmat(0,0) = -1.0*( Kp_np/Kq_np)*NumOfAcini; sysmat(0,1) =  1.0*( Kp_np/Kq_np)*NumOfAcini;
    sysmat(1,0) =  1.0*( Kp_np/Kq_np)*NumOfAcini; sysmat(1,1) = -1.0*( Kp_np/Kq_np)*NumOfAcini;
    
    rhs(0)      = -1.0*(-(Kp_n*(p1n-p2n) + Kp_nm*(p1nm-p2nm))*NumOfAcini/Kq_np + (Kq_n*qn + Kq_nm*qnm)/Kq_np);
    rhs(1)      =  1.0*(-(Kp_n*(p1n-p2n) + Kp_nm*(p1nm-p2nm))*NumOfAcini/Kq_np + (Kq_n*qn + Kq_nm*qnm)/Kq_np);
    
  }
  else if (ele->Type() ==  "Exponential")
  {
    const double Vo  = volAlvDuct;
    double dvnp= (vnp/NumOfAcini)- Vo;
    double dvn = (vn /NumOfAcini)- Vo; 
    
    //------------------------------------------------------------
    // V  = A + B*exp(-K*P)
    //
    // The P-V curve is fitted to create the following
    // P1 = E1.(V-Vo)
    // 
    // E1 = a + b.(V-Vo) + c.exp(d.(V-Vo))
    //------------------------------------------------------------
    
    double kp_np = Rt/(E2*dt)+1;
    double kp_n  =-Rt/(E2*dt);
    double kq_np = Rt*Ra/(E2*dt) + (Ra+Rt);
    double kq_n  =-Rt*Ra/(E2*dt);
    
    double term_nonlin = 0.0;

    //------------------------------------------------------------
    // for now the (a,b,c,d) components are not read from the 
    // input file
    //------------------------------------------------------------
    double a = 6449.0 ;
    double b = 33557.7;
    double c = 6.5158;
    double d = 47.9892;
    
    //------------------------------------------------------------
    // get the terms assosciated with the nonlinear behavior of 
    // E1
    //------------------------------------------------------------
    double pnpi = 0.0;
    double pnpi2= 0.0;
    double dpnpi_dt = 0.0;
    double dpnpi2_dt= 0.0;
    
    // componets of linearized E1
    pnpi      = (a + b*dvnp + c*exp(d*dvnp))*dvnp;
    pnpi2     = (a + 2*b*dvnp + c*exp(d*dvnp)*(d*dvnp+1));
    
    // componets of linearized d(E1)/dt
    dpnpi_dt  = (a+2*b*dvnp+c*exp(d*dvnp)*(1+d*dvnp))*(dvnp-dvn)/dt;
    dpnpi2_dt = (2*b+d*c*exp(d*dvnp)*(1+d*dvnp) + c*d*exp(d*dvnp))*(dvnp-dvn)/dt + (a+2*b*dvnp+c*exp(d*dvnp)*(1+d*dvnp))/dt;
    
#if 0
    cout<<"+------+++++++ DEBUG +++++++------+"<<endl;
    cout<<"pnpi: "<<pnpi<<endl;
    cout<<"dvnp: "<<dvnp<<endl;
    cout<<"pnpi2: "<<pnpi2<<endl;
    cout<<"dpnpi_dt: "<<dpnpi_dt<<endl;
    cout<<"dvn: "<<dvn<<endl;
    cout<<"dt: "<<dt<<endl;
    cout<<"dpnpi2_dt: "<<dpnpi2_dt<<endl;
    cout<<"qn: "<<qn<<endl;
    cout<<"+------+++++++ DEBUG +++++++------+"<<endl;
#endif
    
    term_nonlin = pnpi + pnpi2*(-(dvnp) +(qn/NumOfAcini)*dt/2 + dvn);
    kq_np = kq_np + pnpi2/2*dt;
    term_nonlin = term_nonlin + dpnpi_dt*Rt/E2  + dpnpi2_dt*Rt/E2 *(-(dvnp)+(qnp/NumOfAcini)*dt/2 + dvn);
    kq_np = kq_np + dpnpi2_dt*Rt/E2/2*dt;
    
#if 0
    cout<<"+------------- DEBUG -------------+"<<endl;
    cout<<"NumOfAcini: "<<NumOfAcini<<endl;
    cout<<"kp_np: "<< kp_np<<"\t"<<"kp_n: "<<kp_n<<endl;
    cout<<"kq_np: "<< kq_np<<"\t"<<"kq_n: "<<kq_n<<"\t"<<"term_nonlin: "<<term_nonlin<<endl;
    cout<<"+------------- DEBUG -------------+"<<endl;
#endif


    sysmat(0,0) = -1.0*( kp_np/kq_np)*NumOfAcini;   sysmat(0,1) =  1.0*( kp_np/kq_np)*NumOfAcini; 
    sysmat(1,0) =  1.0*( kp_np/kq_np)*NumOfAcini;   sysmat(1,1) = -1.0*( kp_np/kq_np)*NumOfAcini; 

    rhs(0)      = -1.0*(-(kp_n*(p1n-p2n) - term_nonlin)*NumOfAcini/kq_np +( kq_n*qn)/kq_np);
    rhs(1)      =  1.0*(-(kp_n*(p1n-p2n) - term_nonlin)*NumOfAcini/kq_np +( kq_n*qn)/kq_np);
   
  }
  else if (ele->Type() ==  "DoubleExponential")
  {
    const double Vo  = volAlvDuct;
    double dvnp= (vnp/NumOfAcini)- Vo;
    double dvn = (vn /NumOfAcini)- Vo;
    
    //------------------------------------------------------------
    // V  = A + B*exp(-K*P)
    //
    // The P-V curve is fitted to create the following
    // P1 = E1.(V-Vo)
    //
    // E1 = a + b.(V-Vo) + c.exp(d.(V-Vo))
    //------------------------------------------------------------
    double kp_np = Rt/(E2*dt)+1.0;
    double kp_n  =-Rt/(E2*dt);
    double kq_np = Rt*Ra/(E2*dt) + (Ra+Rt);
    double kq_n  =-Rt*Ra/(E2*dt);
    
    double term_nonlin = 0.0;
    //------------------------------------------------------------
    // for now the (a,b,c,d) components are not read from the
    // input file
    //------------------------------------------------------------
    // OLD VALUES
    //        double a = 6449.0 ;
    //        double b = 33557.7;
    //        double c = 6.5158;
    //        double d = 47.9892;
    
    // NEW VALUES
    double a = 6510.99;
    double b = 3.5228E04;
    double c = 6.97154E-06;
    double d = 144.716;
        
    double a2= 0.0;
    double b2= 0.0;
    double c2= 38000.0*1.4;
    double d2=-90.0000;
    
    //------------------------------------------------------------
    // get the terms assosciated with the nonlinear behavior of
    // E1
    //------------------------------------------------------------
    double pnpi = 0.0;
    double pnpi2= 0.0;
    double dpnpi_dt = 0.0;
    double dpnpi2_dt= 0.0;

    // componets of linearized E1
    pnpi      = (a + b *dvnp + c *exp(d *dvnp))*dvnp;
    pnpi     += (a2+ b2*dvnp + c2*exp(d2*dvnp))*dvnp;
    
    pnpi2     = (a + 2.0*b *dvnp + c *exp(d *dvnp)*(d *dvnp+1.0));
    pnpi2    += (a2+ 2.0*b2*dvnp + c2*exp(d2*dvnp)*(d2*dvnp+1.0));

    // componets of linearized d(E1)/dt
    dpnpi_dt  = (a +2.0*b *dvnp+c *exp(d *dvnp)*(1.0+d *dvnp))*(dvnp-dvn)/dt;
    dpnpi_dt += (a2+2.0*b2*dvnp+c2*exp(d2*dvnp)*(1.0+d2*dvnp))*(dvnp-dvn)/dt;
    
    dpnpi2_dt = (2.0*b +d *c *exp(d *dvnp)*(1.0+d *dvnp) + c *d *exp(d *dvnp))*(dvnp-dvn)/dt + (a +2.0*b *dvnp+c *exp(d *dvnp)*(1.0+d *dvnp))/dt;
    dpnpi2_dt+= (2.0*b2+d2*c2*exp(d2*dvnp)*(1.0+d2*dvnp) + c2*d2*exp(d2*dvnp))*(dvnp-dvn)/dt + (a2+2.0*b2*dvnp+c2*exp(d2*dvnp)*(1.0+d2*dvnp))/dt;
    
    // Add up the nonlinear terms
    

    
    //---------------
    term_nonlin = pnpi + pnpi2*(-(dvnp) +(qn/NumOfAcini)*dt/2.0 + dvn);
    kq_np = kq_np + pnpi2/2.0*dt;
    term_nonlin = term_nonlin + dpnpi_dt*Rt/E2  + dpnpi2_dt*Rt/E2 *(-(dvnp)+(qnp/NumOfAcini)*dt/2.0 + dvn);
    kq_np = kq_np + dpnpi2_dt*Rt/E2/2.0*dt;
    
    sysmat(0,0) = -1.0*( kp_np/kq_np)*NumOfAcini; sysmat(0,1) =  1.0*( kp_np/kq_np)*NumOfAcini;
    sysmat(1,0) =  1.0*( kp_np/kq_np)*NumOfAcini; sysmat(1,1) = -1.0*( kp_np/kq_np)*NumOfAcini;
    
    rhs(0)      = -1.0*((- kp_n*(p1n-p2n) + term_nonlin)*NumOfAcini/kq_np +(kq_n*qn)/kq_np);
    rhs(1)      =  1.0*((- kp_n*(p1n-p2n) + term_nonlin)*NumOfAcini/kq_np +(kq_n*qn)/kq_np);
  }
  else
  {
    dserror("[%s] is not defined as a reduced dimensional lung acinus material",ele->Type().c_str());
    exit(1);
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvaluateTerminalBC(
  RedAcinus*                   ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  vector<int>&                 lm,
  Epetra_SerialDenseVector&    rhs,
  RefCountPtr<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get total time
  const double time = params.get<double>("total time");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // the number of nodes
  const int numnode = lm.size();
  vector<int>::iterator it_vcr;

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BCs
  // ---------------------------------------------------------------------------------

  for(int i = 0; i<ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner()== myrank)
    {
      if(ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond") || ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond") || ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond"))
      {
        string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = *(condition->Get<string>("boundarycond"));
          
        
          const  vector<int>*    curve  = condition->Get<vector<int>    >("curve");
          double curvefac = 1.0;
          const  vector<double>* vals   = condition->Get<vector<double> >("val");
          const vector<int>*     functions = condition->Get<vector<int> >("funct");

          // -----------------------------------------------------------------
          // Read in the value of the applied BC
          // -----------------------------------------------------------------
          if((*curve)[0]>=0)
          {
            curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
            BCin = (*vals)[0]*curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
            exit(1);
          }
          
          int functnum = -1;
          if (functions) functnum = (*functions)[0];
          else functnum = -1;
          
          double functionfac = 0.0;
          if(functnum>0)
          {
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,(ele->Nodes()[i])->X(),time,NULL);
          }
          BCin += functionfac;

          // -----------------------------------------------------------------------------
          // get the local id of the node to whome the bc is prescribed
          // -----------------------------------------------------------------------------
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i]->Id(),discretization.Comm().MyPID());
            exit(1);
          }
        }
        else if (ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
        {
          const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");

          RCP<ParameterList> CoupledTo3DParams  =
            params.get<RCP<ParameterList > >("coupling with 3D fluid params");
          // -----------------------------------------------------------------
          // If the parameter list is empty, then something is wrong!
          // -----------------------------------------------------------------
          if (CoupledTo3DParams.get()==NULL)
          {
            dserror("Cannot prescribe a boundary condition from 3D to reduced D, if the parameters passed don't exist");
            exit(1);
          }
          
          // -----------------------------------------------------------------
          // Read in Condition type
          // -----------------------------------------------------------------
          //        Type = *(condition->Get<string>("CouplingType"));
          // -----------------------------------------------------------------
          // Read in coupling variable rescribed by the 3D simulation
          //
          //     In this case a map called map3D has the following form:
          //     +-----------------------------------------------------------+
          //     |           map< string               ,  double        >    |
          //     |     +------------------------------------------------+    |
          //     |     |  ID  | coupling variable name | variable value |    |
          //     |     +------------------------------------------------+    |
          //     |     |  1   |   flow1                |     0.12116    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  2   |   pressure2            |    10.23400    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     .  .   .   ....                 .     .......    .    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  N   |   variableN            |    value(N)    |    |
          //     |     +------+------------------------+----------------+    |
          //     +-----------------------------------------------------------+
          // -----------------------------------------------------------------
          
          int ID = condition->GetInt("ConditionID");
          RCP<map<string,double> > map3D;
          map3D   = CoupledTo3DParams->get<RCP<map<string,double > > >("3D map of values");
          
          // find the applied boundary variable
          std::stringstream stringID;
          stringID<< "_"<<ID;
          for (map<string,double>::iterator itr = map3D->begin(); itr!=map3D->end(); itr++)
          {
            string VariableWithId = itr->first;
            size_t found;
            found= VariableWithId.rfind(stringID.str());
            if (found!=string::npos)
            {
              Bc   = string(VariableWithId,0,found);
              BCin = itr->second;
              break;
            }
          }

        }
        else if (ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond");
          // Get the type of prescribed bc
          Bc  = *(condition->Get<string>("phase1"));
        
          double period  = condition->GetDouble("period");
          double period1 = condition->GetDouble("phase1_period");

          unsigned int phase_number = 0;

          if (fmod(time,period) > period1)
          {
            phase_number = 1;
            Bc = *(condition->Get<string>("phase2"));
          }

          const  vector<int>*    curve  = condition->Get<vector<int> >("curve");
          double curvefac = 1.0;
          const  vector<double>* vals   = condition->Get<vector<double> >("val");

          // -----------------------------------------------------------------
          // Read in the value of the applied BC
          // -----------------------------------------------------------------
          if((*curve)[phase_number]>=0)
          {
            curvefac = DRT::Problem::Instance()->Curve((*curve)[phase_number]).f(time);
            BCin = (*vals)[phase_number]*curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
            exit(1);
          }

          // -----------------------------------------------------------------------------
          // get the local id of the node to whome the bc is prescribed
          // -----------------------------------------------------------------------------
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i]->Id(),discretization.Comm().MyPID());
            exit(1);
          } 
        }
        else
        {

        }

        if (Bc == "pressure" || Bc == "ExponentialPleuralPressure")
        {
          RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
          RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");
          
          if (bcval==null||dbctog==null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }
          
          if (Bc == "ExponentialPleuralPressure")
          {
            const double ap =  -977.203;
            const double bp = -3338.290;
            const double cp =    -7.686;
            const double dp =  2034.470;
            const double VFR   = 1240000.0;
            const double TLC   = 4760000.0 - VFR;
            
            const double lungVolumenp = params.get<double>("lungVolume_np") - VFR;

            const double TLCnp= lungVolumenp/TLC;

            double Pp_np = ap + bp*exp(cp*TLCnp) + dp*TLCnp;
            
            BCin += Pp_np;
          }          
          
          // set pressure at node i
          int    gid; 
          double val; 
          
          gid = lm[i];
          val = BCin;
          bcval->ReplaceGlobalValues(1,&val,&gid);
      
          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);

        }
        else if (Bc == "flow")
        {
          // ----------------------------------------------------------
          // Since a node might belong to multiple elements then the 
          // flow might be added to the rhs multiple time.
          // To fix this the flow is devided by the number of elements
          // (which is the number of branches). Thus the sum of the
          // final added values is the actual prescribed flow.
          // ----------------------------------------------------------
          int numOfElems = (ele->Nodes()[i])->NumElement();
          BCin /= double(numOfElems);

          // get rhs
          //          RefCountPtr<Epetra_Vector> rhs  = params.get<RCP<Epetra_Vector> >("rhs");
          //          if (rhs==null)
          //          {
          //            dserror("Cannot get state vector 'rhs'");
          //            exit(1);
          //          }
          
          // set pressure at node i
          //          int    gid; 
          //          double val; 
          
          //          gid =  lm[i];
          //          cout<<"FLOW in: "<<BCin<<" with old rhs: "<<rhs(i)<<" With "<<numOfElems<<" elements"<<endl;
          rhs(i) += -BCin + rhs(i);

          //          rhs->ReplaceGlobalValues(1,&val,&gid);
        }
        else
        {
          dserror("precribed [%s] is not defined for reduced acinuss",Bc.c_str());
          exit(1);
        }
        
      }
      
      else if(ele->Nodes()[i]->GetCondition("RedLungAcinusCond"))
      {
        #if 0
        int gid;
        double val;
        RefCountPtr<Epetra_Vector> abc  = params.get<RCP<Epetra_Vector> >("abc");
        gid = lm[i];
        val = 1;
        abc->ReplaceGlobalValues(1,&val,&gid);
        #endif
        // ---------------------------------------------------------------
        // If the node is conected to a reduced dimesnional acinus
        // ---------------------------------------------------------------
        
        // At this state do nothing, since this boundary is resolved
        // during the assembly of Sysmat and RHS
      }
      else
      {
        #if 1
        // ---------------------------------------------------------------
        // If the node is a terminal node, but no b.c is prescribed to it
        // then a zero output pressure is assumed
        // ---------------------------------------------------------------
        if (ele->Nodes()[i]->NumElement() == 1)
        {
          // -------------------------------------------------------------
          // get the local id of the node to whome the bc is prescribed
          // -------------------------------------------------------------

          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
            exit(1);
          }
          
          RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
          RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");
          
          if (bcval==null||dbctog==null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }        
          
          
          // set pressure at node i
          int    gid; 
          double val; 
          
          gid = lm[i];
          val = 0.0;
          bcval->ReplaceGlobalValues(1,&val,&gid);
          
          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);

          //          const double* X = ele->Nodes()[i]->X();
          //          printf("WARNING: node(%d) is free on [%f,%f,%f] \n",gid+1,X[0],X[1],X[2]);
        }
        #endif

        #if 0
        // ---------------------------------------------------------------
        // If the node is a terminal node, but no b.c is prescribed to it
        // then a zero output pressure is assumed
        // ---------------------------------------------------------------
        if (ele->Nodes()[i]->NumElement() == 1)
        {
          RefCountPtr<Epetra_Vector> rhs  = params.get<RCP<Epetra_Vector> >("rhs");
          if (rhs==null)
          {
            dserror("Cannot get state vector 'rhs'");
            exit(1);
          }
          
          // set pressure at node i
          int    gid; 
          double val; 
          
          gid =  lm[i];
          val =  0.0;
          //rhs->ReplaceGlobalValues(1,&val,&gid);
        }
        #endif
      } // END of if there is no BC but the node still is at the terminal
      
    } // END of if node is available on this processor
  } // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::CalcFlowRates(
  RedAcinus*                   ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    elevec1, //a_volumenp,
  Epetra_SerialDenseVector&    elevec2, //a_volume_strain_np,
  vector<int>&                 lm,
  RefCountPtr<MAT::Material>   material)

#if 1
{
  
  //  const int numnode = iel;
  const int elemVecdim = elevec1.Length () ;
  vector<int>::iterator it_vcr;

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get time-step size
  const double dt = params.get<double>("time step size");

  // get time
  const double time = params.get<double>("total time");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------


  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // ---------------------------------------------------------------------
  // get all general state vectors: flow, pressure,
  // ---------------------------------------------------------------------

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  RefCountPtr<const Epetra_Vector> pn   = discretization.GetState("pn");
  RefCountPtr<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  RefCountPtr<Epetra_Vector> qin_nm  = params.get<RefCountPtr<Epetra_Vector> >("qin_nm");
  RefCountPtr<Epetra_Vector> qin_n   = params.get<RefCountPtr<Epetra_Vector> >("qin_n");
  RefCountPtr<Epetra_Vector> qin_np  = params.get<RefCountPtr<Epetra_Vector> >("qin_np");

  RefCountPtr<Epetra_Vector> qout_np = params.get<RefCountPtr<Epetra_Vector> >("qout_np");
  RefCountPtr<Epetra_Vector> qout_n  = params.get<RefCountPtr<Epetra_Vector> >("qout_n");
  RefCountPtr<Epetra_Vector> qout_nm = params.get<RefCountPtr<Epetra_Vector> >("qout_nm");

  RefCountPtr<Epetra_Vector> acinar_vn          = params.get<RefCountPtr<Epetra_Vector> >("acinar_vn");
  RefCountPtr<Epetra_Vector> acinar_vnp         = params.get<RefCountPtr<Epetra_Vector> >("acinar_vnp");
  RefCountPtr<Epetra_Vector> a_volume_strain_np = params.get<RefCountPtr<Epetra_Vector> >("acinar_vnp_strain");


 if (pnp==null || pn==null || pnm==null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // extract local values from the global vectors
  vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  // extract local values from the global vectors
  vector<double> mypnm(lm.size());
  DRT::UTILS::ExtractMyValues(*pnm,mypnm,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(elemVecdim);
  Epetra_SerialDenseVector epn (elemVecdim);
  Epetra_SerialDenseVector epnm(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)   = mypnp[i];
    epn(i)    = mypn[i];
    epnm(i)   = mypnm[i];
  }

  double e_acin_vnp = 0.0;
  double e_acin_vn = 0.0;

  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    e_acin_vnp = (*acinar_vnp)[ele->LID()];
    e_acin_vn  = (*acinar_vn )[ele->LID()];
  }


  // get the volumetric flow rate from the previous time step
  ParameterList elem_params;
  elem_params.set<double>("qout_np",(*qout_np)[ele->LID()]);
  elem_params.set<double>("qout_n" ,(*qout_n )[ele->LID()]);
  elem_params.set<double>("qout_nm",(*qout_nm)[ele->LID()]);
  elem_params.set<double>("qin_np" ,(*qin_np )[ele->LID()]);
  elem_params.set<double>("qin_n"  ,(*qin_n  )[ele->LID()]);
  elem_params.set<double>("qin_nm" ,(*qin_nm )[ele->LID()]);

  elem_params.set<double>("acin_vnp" ,e_acin_vnp);
  elem_params.set<double>("acin_vn"  ,e_acin_vn );

  Epetra_SerialDenseMatrix sysmat (elemVecdim, elemVecdim,true);
  Epetra_SerialDenseVector  rhs (elemVecdim);


  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         epnp,
         epn,
         epnm,
         sysmat,
         rhs,
         material,
         elem_params,
         time,
         dt);


  double qn = (*qin_n  )[ele->LID()];
  double qnp= -1.0*(sysmat(0,0)*epnp(0) + sysmat(0,1)*epnp(1) - rhs(0));
#if 0
  cout<<"Calculating FLOWRATE (Q): "<<qnp<<endl;
  cout<<"["<<ele->LID()<<"] Q calc: -("<<sysmat(0,0)<<"*"<<epnp(0)<<"+"<<sysmat(0,1)<<"*"<<epnp(1)<<"-"<<rhs(0)<<")"<<endl;
  cout<<"p1: "<<epnp(0)<<" \tp2: "<<epnp(1)<<endl;
  cout<<"rhs=("<<rhs(0)<<","<<rhs(1)<<")"<<endl;
  cout<<"Q>>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"Q sysmat: "<<sysmat<<endl;
  cout<<"Q>>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"Q rhs: "<<rhs<<endl;
  cout<<"Q>>>>>>>>>>>>>>>>>===================>>>>>>>>>>>>>>>>>"<<endl;
#endif


  //  ele->setVars("flow_in",qin);
  //  ele->setVars("flow_out",qout);
  int gid = ele->Id();

  qin_np  -> ReplaceGlobalValues(1,&qnp,&gid);
  qout_np -> ReplaceGlobalValues(1,&qnp,&gid);
  {
    double acinus_volume = e_acin_vn;
    acinus_volume +=  0.5*(qnp+qn)*dt;
    acinar_vnp-> ReplaceGlobalValues(1,& acinus_volume,&gid);
    
    //----------------------------------------------------------------
    // Read in the material information
    //----------------------------------------------------------------
    double vo = 0.0;
    ele->getParams("AcinusVolume",vo);
    double avs_np = (acinus_volume - vo)/vo;
    a_volume_strain_np -> ReplaceGlobalValues(1,&avs_np,&gid);
  }


}
#else
{
  //  const int   myrank  = discretization.Comm().MyPID();

  //  const int numnode = iel;

  RefCountPtr<const Epetra_Vector> pnp   = discretization.GetState("pnp");
  RefCountPtr<const Epetra_Vector> pn    = discretization.GetState("pn");
  RefCountPtr<const Epetra_Vector> acinar_vn = discretization.GetState("acinar_vn");

  RefCountPtr<Epetra_Vector> qin_np      = params.get<RefCountPtr<Epetra_Vector> >("qin_np");
  RefCountPtr<Epetra_Vector> qout_np     = params.get<RefCountPtr<Epetra_Vector> >("qout_np");
  RefCountPtr<Epetra_Vector> qin_n       = params.get<RefCountPtr<Epetra_Vector> >("qin_n");
  RefCountPtr<Epetra_Vector> qout_n      = params.get<RefCountPtr<Epetra_Vector> >("qout_n");

  // get time-step size
  const double dt = params.get<double>("time step size");

  double dens = 0.0;
  double visc = 0.0;

  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    // get actual material
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

    // get density
    dens = actmat->Density();
    
    // get dynamic viscosity
    visc = actmat->Viscosity();
  }
  else
  {
    dserror("Material law is not a Newtonia fluid");
    exit(1);
  }

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  vector<double> mypn (lm.size());
  vector<double> myacinar_vn (lm.size());

  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);
  DRT::UTILS::ExtractMyValues(*pn ,mypn ,lm);
  DRT::UTILS::ExtractMyValues(*acinar_vn,myacinar_vn,lm);

  const int numnode = lm.size();

  // create objects for element arrays
  //  LINALG::Matrix<numnode,1> epnp;
  Epetra_SerialDenseVector epnp(numnode);
  Epetra_SerialDenseVector epn(numnode);
  Epetra_SerialDenseVector a_volumen(numnode);
  for (unsigned int i=0;i<lm.size();++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)      = mypnp[i];
    epn(i)       = mypn[i];
    a_volumen(i) = myacinar_vn[i];    
  }


  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();

  LINALG::Matrix<3,iel> xyze;
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }


  // Calculate the length of acinus element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));


  //--------------------------------------------------------------
  //               Get the elements flow in/out
  //--------------------------------------------------------------
  double eqin_n   = (*qin_n  )[ele->LID()];
  double eqout_n  = (*qout_n )[ele->LID()];
  double eqin_np  = (*qin_np )[ele->LID()];
  double eqout_np = (*qout_np)[ele->LID()];

  // get the generation number
  int generation = 0;
  ele->getParams("Generation",generation);

  double R = -1.0;

  // get element information
  double A;
  ele->getParams("Area",A);
  // evaluate Poiseuille resistance
  double Rp = 8.0*PI*visc*L/(pow(A,2));
  // evaluate the Reynolds number
  const double Re = 2.0*fabs(eqout_np)/(visc/dens*sqrt(A*PI));
  if (ele->Resistance() == "Poiseuille")
  {
    R = Rp;
  }
  else if (ele->Resistance() == "Pedley")
  {
    //-----------------------------------------------------------------
    // resistance evaluated using Pedley's model from :
    // Pedley et al (1970)
    //-----------------------------------------------------------------
    double gamma = 0.327;
    R  = gamma* (sqrt(Re * 2.0*sqrt(A/PI)/L)) * Rp;
    
    //-----------------------------------------------------------------
    // Correct any resistance smaller than Poiseuille's one
    //-----------------------------------------------------------------
    //  if (R < Rp)
    //    {
    //      R = Rp;
    //    }
    double alfa = sqrt(2.0*sqrt(A/PI)/L);
    double Rep  = 1.0/((gamma*alfa)*(gamma*alfa));
    double k    = 0.50;
    double st   = 1.0/(1.0+exp(-2*k*(Re-Rep)));
    
    R = R*st + Rp*(1.0-st);

  }
  else if (ele->Resistance() == "Generation_Dependent_Pedley")
  {
    //-----------------------------------------------------------------
    // Gamma is taken from Ertbruggen et al
    //-----------------------------------------------------------------
    double gamma = 0.327;
    switch(generation)
    {
    case 0: 
      gamma = 0.162;
      break;
    case 1: 
      gamma = 0.239;
      break;
    case 2: 
      gamma = 0.244;
      break;
    case 3: 
      gamma = 0.295;
      break;
    case 4: 
      gamma = 0.175;
      break;
    case 5: 
      gamma = 0.303;
      break;
    case 6: 
      gamma = 0.356;
      break;
    case 7: 
      gamma = 0.566;
      break;
    default:
      gamma = 0.327;
    }
    //-----------------------------------------------------------------
    // resistance evaluated using Pedley's model from :
    // Pedley et al (1970)
    //-----------------------------------------------------------------
    R  = gamma* (sqrt(Re * 2.0*sqrt(A/PI)/L)) * Rp;
    
    //-----------------------------------------------------------------
    // Correct any resistance smaller than Poiseuille's one
    //-----------------------------------------------------------------
    if (R < Rp)
    {
      R = Rp;
    }
  }
  else if (ele->Resistance() == "Reynolds")
  {
    R = Rp *(3.4+2.1e-3 *Re);
  }
  else
  {
    dserror("[%s] is not a defined resistance model",ele->Resistance().c_str());
  }
  
  // ------------------------------------------------------------------
  // Find the acinus type
  // ------------------------------------------------------------------
  if(ele->Type() == "Resistive")
  {
    eqin_np = eqout_np = (epnp(0)-epnp(1))/R;
  }
  else if(ele->Type() == "InductoResistive")
  {
    // get element information
    double  Ew, tw;
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);
  }
  else if(ele->Type() == "ComplientResistive")
  {
    // get element information
    double Ew, tw;
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

 
    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // get element information from the previous time step
    double qcn =  eqin_n - eqout_n;

    // -----------------------------------------------------------
    // calculate capacitance flow at time step n+1
    // -----------------------------------------------------------
    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - qcn;

    eqout_np = (epnp(0)-epnp(1))/R;

    eqin_np  = eqout_np + qcnp_val;
    
  }
  else if(ele->Type() == "RLC")
  {
    // get element information
    double Ew, tw;
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // find Inductance I
    //    const double I = dens*L/A;

    // get element information from the previous time step
    double qcn,qln;
    qcn = eqin_n - eqout_n;

    // -----------------------------------------------------------
    // calculate the inductance flow at time step n and n+1
    // qln = qRn (qRn: flow in resistor at n)
    // -----------------------------------------------------------
    double qlnp = 0.0;
    qln  = eqout_n;
    qlnp = (epnp(2) - epnp(1))/R;

    // -----------------------------------------------------------
    // calculate capacitance flow at time step n+1
    // -----------------------------------------------------------
#if 1
    double qcnp_val = (C/dt)*(epnp(0)-epn(0));
#else
    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - qcn;
#endif

    
    eqout_np = qlnp;
    eqin_np  = eqout_np + qcnp_val;
  }
  else if(ele->Type() == "SUKI")
  {
    
  }
  else
  {
    dserror("[%s] is not an implimented elements yet",(ele->Type()).c_str());
    exit(1);
  }

  //  ele->setVars("flow_in",qin);
  //  ele->setVars("flow_out",qout);
  int gid = ele->Id();

  qin_np  -> ReplaceGlobalValues(1,&eqin_np ,&gid);
  qout_np -> ReplaceGlobalValues(1,&eqout_np,&gid);

  for (int i = 0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("RedLungAcinusCond"))
    {
      double acinus_volume = a_volumen(i);
      acinus_volume +=  0.5*(eqout_np+eqout_n)*dt;
      a_volumenp(i) = acinus_volume;

      //----------------------------------------------------------------
      // Read in the material information
      //----------------------------------------------------------------
      const double VolPerArea = (ele->Nodes()[i]->GetCondition("RedLungAcinusCond"))->GetDouble("VolumePerArea");
      const double VolAcinus  = (ele->Nodes()[i]->GetCondition("RedLungAcinusCond"))->GetDouble("Acinus_Volume");
      double vo =  VolPerArea*A;
      a_volume_strain_np(i) =  (acinus_volume - vo)/vo;
    }
  }

}//CalcFlowRates
#endif

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::GetCoupledValues(
  RedAcinus*                   ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  vector<int>&                 lm,
  RefCountPtr<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get total time
  //  const double time = params.get<double>("total time");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // the number of nodes
  const int numnode = lm.size();
  vector<int>::iterator it_vcr;

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BCs
  // ---------------------------------------------------------------------------------
  for(int i = 0; i<ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner()== myrank)
    {
      if(ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
      {

          const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");        
          RCP<ParameterList> CoupledTo3DParams  =
            params.get<RCP<ParameterList > >("coupling with 3D fluid params");
          // -----------------------------------------------------------------
          // If the parameter list is empty, then something is wrong!
          // -----------------------------------------------------------------
          if (CoupledTo3DParams.get()==NULL)
          {
            dserror("Cannot prescribe a boundary condition from 3D to reduced D, if the parameters passed don't exist");
            exit(1);
          }


        // -----------------------------------------------------------------
        // Compute the variable solved by the reduced D simulation to be
        // passed to the 3D simulation
        //
        //     In this case a map called map1D has the following form:
        //     +-----------------------------------------------------------+
        //     |              map< string            ,  double        > >  |
        //     |     +------------------------------------------------+    |
        //     |     |  ID  | coupling variable name | variable value |    |
        //     |     +------------------------------------------------+    |
        //     |     |  1   |   flow1                |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  2   |   pressure2            |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     .  .   .   ....                 .     .......    .    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  N   |   variable(N)          | trash value(N) |    |
        //     |     +------+------------------------+----------------+    |
        //     +-----------------------------------------------------------+
        // -----------------------------------------------------------------

        int ID = condition->GetInt("ConditionID");
        RCP<map<string,double> >  map1D;
        map1D   = CoupledTo3DParams->get<RCP<map<string,double> > >("reducedD map of values");

        string returnedBC = *(condition->Get<string>("ReturnedVariable"));

        double BC3d = 0.0;
        if (returnedBC  == "flow")
        {
          // MUST BE DONE
        }
        else if (returnedBC == "pressure")
        {
          BC3d     = epnp(i);
        }
        else
        {
          string str = (*condition->Get<string>("ReturnedVariable"));
          dserror("%s, is an unimplimented type of coupling",str.c_str());
          exit(1);
        }
        std::stringstream returnedBCwithId;
        returnedBCwithId << returnedBC <<"_" << ID;
        
        //        cout<<"Return ["<<returnedBC<<"] form 1D problem to 3D SURFACE of ID["<<ID<<"]: "<<BC3d<<endl;

        // -----------------------------------------------------------------
        // Check whether the coupling wrapper has already initialized this
        // map else wise we will have problems with parallelization, that's
        // because of the preassumption that the map is filled and sorted
        // Thus we can use parallel addition
        // -----------------------------------------------------------------
        
        map<string,double>::iterator itrMap1D;
        itrMap1D = map1D->find(returnedBCwithId.str());
        if (itrMap1D == map1D->end())
        {
          dserror("The 3D map for (1D - 3D coupling) has no variable (%s) for ID [%d]",returnedBC.c_str(),ID );
          exit(1);
        }
        
        // update the 1D map
        (*map1D)[returnedBCwithId.str()] = BC3d;
      }   
    } // END of if node is available on this processor
  } // End of node i has a condition
}

