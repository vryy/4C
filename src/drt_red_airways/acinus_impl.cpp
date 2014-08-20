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
#include "../drt_mat/air_0d_O2_saturation.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matlist.H"
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
  //  const int   myrank  = discretization.Comm().MyPID();

  //  const int numnode = iel;
  const int elemVecdim = elevec1_epetra.Length () ;
  std::vector<int>::iterator it_vcr;

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

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  Teuchos::RCP<const Epetra_Vector> pn   = discretization.GetState("pn");
  Teuchos::RCP<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  Teuchos::RCP<const Epetra_Vector> ial  = discretization.GetState("intr_ac_link");

  Teuchos::RCP<Epetra_Vector> acinar_vnp  = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> acinar_vn   = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");


  Teuchos::RCP<Epetra_Vector> qin_nm  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_nm");
  Teuchos::RCP<Epetra_Vector> qin_n   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_n");
  Teuchos::RCP<Epetra_Vector> qin_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  Teuchos::RCP<Epetra_Vector> qout_n  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_n");
  Teuchos::RCP<Epetra_Vector> qout_nm = params.get<Teuchos::RCP<Epetra_Vector> >("qout_nm");


  if (pnp==Teuchos::null || pn==Teuchos::null || pnm==Teuchos::null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  // extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
  DRT::UTILS::ExtractMyValues(*pnm,mypnm,lm);

  // extract local values from the global vectors
  std::vector<double> myial(lm.size());
  DRT::UTILS::ExtractMyValues(*ial,myial,lm);

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
  Teuchos::ParameterList elem_params;
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


  double Ao = 0.0;
  ele->getParams("Area",Ao);
  if (myial[1] > 0.0)
  {
    elemat1_epetra(1,0)=0.0;
    elemat1_epetra(1,1)=0.0;
    elevec1_epetra(1)  =0.0;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::Initial(
  RedAcinus*                             ele,
  Teuchos::ParameterList&                params,
  DRT::Discretization&                   discretization,
  std::vector<int>&                      lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  const int   myrank  = discretization.Comm().MyPID();

  Teuchos::RCP<Epetra_Vector> p0np    = params.get<Teuchos::RCP<Epetra_Vector> >("p0np");
  Teuchos::RCP<Epetra_Vector> p0n     = params.get<Teuchos::RCP<Epetra_Vector> >("p0n");
  Teuchos::RCP<Epetra_Vector> p0nm    = params.get<Teuchos::RCP<Epetra_Vector> >("p0nm");

  Teuchos::RCP<Epetra_Vector> generations   = params.get<Teuchos::RCP<Epetra_Vector> >("generations");
  Teuchos::RCP<Epetra_Vector> a_bc          = params.get<Teuchos::RCP<Epetra_Vector> >("acini_bc");

  //  Teuchos::RCP<Epetra_Vector> a_volume      = params.get<Teuchos::RCP<Epetra_Vector> >("acini_volume");
  Teuchos::RCP<Epetra_Vector> a_e_volume    = params.get<Teuchos::RCP<Epetra_Vector> >("acini_e_volume");

  //  std::vector<int>::iterator it = lm.begin();

  //std::vector<int> lmowner;
  std::vector<int> lmstride;
  Teuchos::RCP<std::vector<int> > lmowner = Teuchos::rcp(new std::vector<int>);
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

  // find the volume of an acinus conditions
  {
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
    if(ele->Nodes()[i]->GetCondition("RedAirwayEvalLungVolCond"))
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

  bool solveScatra = params.get<bool>("solveScatra");
  Teuchos::RCP<Epetra_Vector>  junVolMix_Corrector;
  Teuchos::RCP<Epetra_Vector>  scatranp;
  Teuchos::RCP<Epetra_Vector>  e1scatranp;
  Teuchos::RCP<Epetra_Vector>  e2scatranp;

  if (solveScatra)
  {
    junVolMix_Corrector  = params.get<Teuchos::RCP<Epetra_Vector> >("junVolMix_Corrector");
    scatranp             = params.get<Teuchos::RCP<Epetra_Vector> >("scatranp");
    e1scatranp           = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
    e2scatranp           = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");
  }

  if (solveScatra)
  {
    double A=0.0;
    double V=0.0;
    ele->getParams("AcinusVolume",V);
    ele->getParams("Area",A);
    int    gid = lm[1];
    junVolMix_Corrector->ReplaceGlobalValues(1,&A,&gid);

    for (int sci=0;sci<iel;sci++)
    {
      int sgid = lm[sci];
      int    esgid = ele->Id();
      // -------------------------------------------------------------
      //
      // -------------------------------------------------------------
      if (ele->Nodes()[sci]->GetCondition("RedAirwayScatraAirCond"))
      {
        double intSat = ele->Nodes()[sci]->GetCondition("RedAirwayScatraAirCond")->GetDouble("INITIAL_CONCENTRATION");
        int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
        // check if O2 properties material exists
        if (id==-1)
        {
          dserror("A material defining O2 properties in air could not be found");
          exit(1);
        }
        const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
        const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

        // get atmospheric pressure
        double patm = actmat->atmospheric_p_;
        // get number of O2 moles per unit volume of O2
        double nO2perVO2 = actmat->nO2_per_VO2_;

        // calculate the PO2 at nodes
        double pO2 = intSat*patm;

        // calculate VO2
        double vO2 = V*(pO2/patm);
        // evaluate initial concentration
        double intConc = nO2perVO2*vO2/V;

        scatranp->ReplaceGlobalValues(1,&intConc,&sgid);
        e1scatranp->ReplaceGlobalValues(1,&intConc,&esgid);
        e2scatranp->ReplaceGlobalValues(1,&intConc,&esgid);
      }
      else
      {
        dserror("0D Acinus scatra must be predefined as \"air\" only");
        exit(1);
      }
    }
  }
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
  Teuchos::ParameterList &                          params,
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
  }
  else
  {
    dserror("Material law is not a Newtonia fluid");
    exit(1);
  }

  // check here, if we really have an acinus !!

  rhs.Scale(0.0);
  sysmat.Scale(0.0);


  double q_out    = params.get<double>("qout_n");

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
  //  std::cout<<"Acinus Vol: "<<VolAcinus<<std::endl;
  //  std::cout<<"AlvDuc Vol: "<<volAlvDuct<<std::endl;
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

    ele->getParams("E1_0",a);
    ele->getParams("E1_LIN",b);
    ele->getParams("E1_EXP",c);
    ele->getParams("TAU",d);

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

    term_nonlin = pnpi + pnpi2*(-(dvnp) +(qn/NumOfAcini)*dt/2 + dvn);
    kq_np = kq_np + pnpi2/2*dt;
    term_nonlin = term_nonlin + dpnpi_dt*Rt/E2  + dpnpi2_dt*Rt/E2 *(-(dvnp)+(qnp/NumOfAcini)*dt/2 + dvn);
    kq_np = kq_np + dpnpi2_dt*Rt/E2/2*dt;

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

    ele->getParams("E1_01"  ,a);
    ele->getParams("E1_LIN1",b);
    ele->getParams("E1_EXP1",c);
    ele->getParams("TAU1"   ,d);

    ele->getParams("E1_02"  ,a2);
    ele->getParams("E1_LIN2",b2);
    ele->getParams("E1_EXP2",c2);
    ele->getParams("TAU2"   ,d2);

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
  /* Acinus Type "VolumetricOgden": continuum mechanics derivation of cauchy stress (=hydrostatic pressure)
     for Ogden material for purely volumetric deformation                                  (croth 08/2014)*/
  else if (ele->Type() ==  "VolumetricOgden")
  {
    //Variables for acinus
    const double Vo  = volAlvDuct;
    double vi_n  = (acin_vn/NumOfAcini);
    double qi_n  = (qn /NumOfAcini);
    double qi_np = (qnp/NumOfAcini);

    //Parameters kappa and beta
    double kappa = 0.0;
    double beta  = 0.0;
    ele->getParams("kappa",kappa);
    ele->getParams("beta",beta);

    //Linear branches of the Maxwell model (E2, B=R_t, B_a=R_a), notation according to interacinar dependency paper
    double Kp_np  =  Rt/(E2*dt)+1.0;
    double Kp_n   = -Rt/(E2*dt);
    double Kq_np  =  Rt*Ra/(E2*dt)+Rt+Ra;
    double Kq_n   = -Rt*Ra/(E2*dt);
    double rhsLin = -Kp_n*(p1n-p2n) + Kq_n*qi_n;

    //Branch E_1 of the Maxwell model: Hydrostatic pressure (=Cauchy stress) for Ogden material
    // P_1  = P_c + P_d
    // where P_c = (kappa/beta)*(lambda^(-3))
    //       P_d =-(kappa/beta)*(lambda^(-3-3*beta))
    //       \lambda is the volumetric strain ratio, \lambda = (V/Vo)^(1/3)
    double vi_np = qi_np*dt+vi_n;
    double Kq_npNL = (Rt/E2) * (-kappa*Vo/(pow(vi_np,2.0)*beta) +
                                   (beta+1.0)*kappa*(pow(Vo/(vi_np),beta+1.0))/((vi_np)*beta));
    double rhsNL   = (Vo/vi_n) * (kappa/beta) * (1-pow((Vo/vi_n),beta));

    //To Do: FD-Check

    //add linearisation part to system matrix
    Kq_np += Kq_npNL;

    //Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
    sysmat(0,0) = -1.0*( Kp_np/Kq_np)*NumOfAcini;    sysmat(0,1) =  1.0*( Kp_np/Kq_np)*NumOfAcini;
    sysmat(1,0) =  1.0*( Kp_np/Kq_np)*NumOfAcini;    sysmat(1,1) = -1.0*( Kp_np/Kq_np)*NumOfAcini;

    //Build the corresponding right hand side
    rhs(0)      = -1.0*((rhsLin+rhsNL)*NumOfAcini/Kq_np);
    rhs(1)      =  1.0*((rhsLin+rhsNL)*NumOfAcini/Kq_np);
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
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Epetra_SerialDenseVector&    rhs,
  Teuchos::RCP<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get total time
  const double time = params.get<double>("total time");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // the number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==Teuchos::null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
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
        std::string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = *(condition->Get<std::string>("boundarycond"));

          const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
          double curvefac = 1.0;
          const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");
          const std::vector<int>*     functions = condition->Get<std::vector<int> >("funct");

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

          // get curve2
          int curve2num = -1;
          double curve2fac = 1.0;
          if (curve) curve2num = (*curve)[1];
          if (curve2num>=0 )
            curve2fac = DRT::Problem::Instance()->Curve(curve2num).f(time);

          BCin += functionfac*curve2fac;

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

          Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams  =
            params.get<Teuchos::RCP<Teuchos::ParameterList > >("coupling with 3D fluid params");
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
          //        Type = *(condition->Get<std::string>("CouplingType"));
          // -----------------------------------------------------------------
          // Read in coupling variable rescribed by the 3D simulation
          //
          //     In this case a map called map3D has the following form:
          //     +-----------------------------------------------------------+
          //     |           std::map< string               ,  double        >    |
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
          Teuchos::RCP<std::map<std::string,double> > map3D;
          map3D   = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string,double > > >("3D map of values");

          // find the applied boundary variable
          std::stringstream stringID;
          stringID<< "_"<<ID;
          for (std::map<std::string,double>::iterator itr = map3D->begin(); itr!=map3D->end(); itr++)
          {
            std::string VariableWithId = itr->first;
            size_t found;
            found= VariableWithId.rfind(stringID.str());
            if (found!=std::string::npos)
            {
              Bc   = std::string(VariableWithId,0,found);
              BCin = itr->second;
              break;
            }
          }

        }
        else if (ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond");
          // Get the type of prescribed bc
          Bc  = *(condition->Get<std::string>("phase1"));

          double period  = condition->GetDouble("period");
          double period1 = condition->GetDouble("phase1_period");

          unsigned int phase_number = 0;

          if (fmod(time,period) > period1)
          {
            phase_number = 1;
            Bc = *(condition->Get<std::string>("phase2"));
          }

          const  std::vector<int>*    curve  = condition->Get<std::vector<int> >("curve");
          double curvefac = 1.0;
          const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

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

        if (Bc == "pressure" || Bc == "VolumeDependentPleuralPressure")
        {
          Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

          if (bcval==Teuchos::null||dbctog==Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          if (Bc == "VolumeDependentPleuralPressure")
          {
            DRT::Condition * pplCond = ele->Nodes()[i]->GetCondition("RedAirwayVolDependentPleuralPressureCond");
            double Pp_np = 0.0;
            if (pplCond)
            {
              const  std::vector<int>*    curve  = pplCond->Get<std::vector<int>    >("curve");
              double curvefac = 1.0;
              const  std::vector<double>* vals   = pplCond->Get<std::vector<double> >("val");

              // -----------------------------------------------------------------
              // Read in the value of the applied BC
              //
              // -----------------------------------------------------------------
              if((*curve)[0]>=0)
              {
                curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
              }

              std::string ppl_Type = *(pplCond->Get<std::string>("TYPE"));

              double ap = pplCond->GetDouble("P_PLEURAL_0");
              double bp = pplCond->GetDouble("P_PLEURAL_LIN");
              double cp = pplCond->GetDouble("P_PLEURAL_NONLIN");
              double dp = pplCond->GetDouble("TAU");
              double RV    = pplCond->GetDouble("VFR");
              double TLC   = pplCond->GetDouble("TLC");
              const double lungVolumenp = params.get<double>("lungVolume_n") - RV;
              const double TLCnp= lungVolumenp/(TLC-RV);


              if (ppl_Type == "Exponential")
              {
                Pp_np = ap +  bp*TLCnp+ cp*exp(dp*TLCnp);
              }
              else if (ppl_Type == "Polynomial")
              {
                Pp_np = -pow(1.0/(TLCnp+RV/(TLC-RV)),dp)*cp + bp*TLCnp + ap;
              }
              else
              {
                dserror("Unknown volume pleural pressure type: %s",ppl_Type.c_str());
              }
              Pp_np *= curvefac*((*vals)[0]);
            }
            else
            {
              dserror("No volume dependent pleural pressure condition was defined");
            }
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
          //          Teuchos::RCP<Epetra_Vector> rhs  = params.get<Teuchos::RCP<Epetra_Vector> >("rhs");
          //          if (rhs==Teuchos::null)
          //          {
          //            dserror("Cannot get state vector 'rhs'");
          //            exit(1);
          //          }

          // set pressure at node i
          //          int    gid;
          //          double val;

          //          gid =  lm[i];
          //          std::cout<<"FLOW in: "<<BCin<<" with old rhs: "<<rhs(i)<<" With "<<numOfElems<<" elements"<<std::endl;
          rhs(i) += -BCin + rhs(i);

          //          rhs->ReplaceGlobalValues(1,&val,&gid);
        }
        else
        {
          dserror("precribed [%s] is not defined for reduced acinuss",Bc.c_str());
          exit(1);
        }

      }
      else
      {
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

          Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

          if (bcval==Teuchos::null||dbctog==Teuchos::null)
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
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)

{

  //  const int numnode = iel;
  const int elemVecdim = lm.size() ;
//  std::vector<int>::iterator it_vcr;

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

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  Teuchos::RCP<const Epetra_Vector> pn   = discretization.GetState("pn");
  Teuchos::RCP<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  Teuchos::RCP<Epetra_Vector> qin_nm  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_nm");
  Teuchos::RCP<Epetra_Vector> qin_n   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_n");
  Teuchos::RCP<Epetra_Vector> qin_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  Teuchos::RCP<Epetra_Vector> qout_n  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_n");
  Teuchos::RCP<Epetra_Vector> qout_nm = params.get<Teuchos::RCP<Epetra_Vector> >("qout_nm");

  Teuchos::RCP<Epetra_Vector> acinar_vn          = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> a_volume_strain_np = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp_strain");


  if (pnp==Teuchos::null || pn==Teuchos::null || pnm==Teuchos::null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  // extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
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
  Teuchos::ParameterList elem_params;
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


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/13|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::CalcElemVolume(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)

{
  // get all essential vector variables
  Teuchos::RCP<Epetra_Vector> acinar_vnp   = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> elemVolumenp = params.get<Teuchos::RCP<Epetra_Vector> >("elemVolumenp");

  // get acinus size
  double evolnp = (*elemVolumenp)[ele->LID()];

  // get element global ID
  int gid = ele->Id();

  // update elem
  elemVolumenp->ReplaceGlobalValues(1,& evolnp,&gid);

}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::GetCoupledValues(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get total time
  //  const double time = params.get<double>("total time");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // the number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==Teuchos::null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
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
          Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams  =
            params.get<Teuchos::RCP<Teuchos::ParameterList > >("coupling with 3D fluid params");
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
        //     |              std::map< string            ,  double        > >  |
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
        Teuchos::RCP<std::map<std::string,double> >  map1D;
        map1D   = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string,double> > >("reducedD map of values");

        std::string returnedBC = *(condition->Get<std::string>("ReturnedVariable"));

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
          std::string str = (*condition->Get<std::string>("ReturnedVariable"));
          dserror("%s, is an unimplimented type of coupling",str.c_str());
          exit(1);
        }
        std::stringstream returnedBCwithId;
        returnedBCwithId << returnedBC <<"_" << ID;

        //        std::cout<<"Return ["<<returnedBC<<"] form 1D problem to 3D SURFACE of ID["<<ID<<"]: "<<BC3d<<std::endl;

        // -----------------------------------------------------------------
        // Check whether the coupling wrapper has already initialized this
        // map else wise we will have problems with parallelization, that's
        // because of the preassumption that the map is filled and sorted
        // Thus we can use parallel addition
        // -----------------------------------------------------------------

        std::map<std::string,double>::iterator itrMap1D;
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


/*----------------------------------------------------------------------*
 |  calculate the ammount of fluid mixing inside a          ismail 02/13|
 |  junction                                                            |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::GetJunctionVolumeMix(RedAcinus*                   ele,
                                                              Teuchos::ParameterList&      params,
                                                              DRT::Discretization&         discretization,
                                                              Epetra_SerialDenseVector&    volumeMix_np,
                                                              std::vector<int>&            lm,
                                                              Teuchos::RCP<MAT::Material>           material)
{
  // get flow rate out at time step n+1
  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  // get acinar volumes at time step n+1
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");

  // get the element qout
  double q_out = (*qout_np)[ele->LID()];

  // if transport is flowing into the acinus
  if( q_out >= 0.0)
  {
    ele->getParams("Area",volumeMix_np(1));
  }
  // els if transport is flowing out of the acinus
  else
  {
    ele->getParams("Area",volumeMix_np(0));
    ele->getParams("Area",volumeMix_np(1));
  }

  // extra treatment if an acinus is not connected to anything else
  for (int i=0;i<iel;i++)
  {
    if( ele->Nodes()[i]->NumElement()==1)
      ele->getParams("Area",volumeMix_np(i));
  }
}


/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::SolveScatra(RedAcinus*                   ele,
                                                     Teuchos::ParameterList&      params,
                                                     DRT::Discretization&         discretization,
                                                     Epetra_SerialDenseVector&    scatranp,
                                                     Epetra_SerialDenseVector&    volumeMix_np,
                                                     std::vector<int>&            lm,
                                                     Teuchos::RCP<MAT::Material>           material)
{
  const int   myrank  = discretization.Comm().MyPID();
  Teuchos::RCP<Epetra_Vector> qin_np   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");

  Teuchos::RCP<Epetra_Vector> e1scatran = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatran");
  Teuchos::RCP<Epetra_Vector> e2scatran = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatran");

  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  Teuchos::RCP<Epetra_Vector> acinar_vn          = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");

  double volumenp = (*acinar_vnp)[ele->LID()];
  double volumen = (*acinar_vn )[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);
  // get area
  double area = myvolmix[1];

  // get the elements Qin and Qout
  double q_out = (*qout_np)[ele->LID()];
  double q_in  = (*qin_np)[ele->LID()];
  double e1s   = (*e1scatran)[ele->LID()];
  double e2s   = (*e2scatran)[ele->LID()];

  // get time step size
  // const double dt = params.get<double>("time step size");

  // get time
  const double time = params.get<double>("total time");

  //--------------------------------------------------------------------
  // get element length
  //--------------------------------------------------------------------

  // evaluate velocity at nodes (1) and (2)
  double vel1 = q_in /area;
  double vel2 = q_out/area;

  LINALG::Matrix<2,1> velv;
  velv(0,0)=vel1;
  velv(1,0)=vel2;
  // get average velocity
  double vel = vel2;

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------

  if (vel >=0.0)
  {
    // extrapolate the analytical solution
    int gid = ele->Id();
    double scnp = 0.0;
    scnp = (e2s*volumen + e1s*(volumenp-volumen))/(volumenp);
    e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
  else
  {
    // extrapolate the analytical solution
    double scnp = 0.0;
    int gid = ele->Id();
    scnp = (e2s*volumen + e2s*(volumenp-volumen))/(volumenp);
    {
      e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
      e1scatranp->ReplaceGlobalValues(1,&scnp,&gid);
    }
  }

  for (int i=0;i<2;i++)
  {
    if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedScatraCond") &&     myrank == ele->Nodes()[i]->Owner())
    {
      double scnp =0.0;
      DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedScatraCond");
      // Get the type of prescribed bc

      const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
      double curvefac = 1.0;
      const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

      // -----------------------------------------------------------------
      // Read in the value of the applied BC
      // -----------------------------------------------------------------
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      if (curvenum>=0 )
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

      scnp = (*vals)[0]*curvefac;

      const std::vector<int>*    functions = condition->Get<std::vector<int> >("funct");
      int functnum = -1;
      if (functions) functnum = (*functions)[0];
      else functnum = -1;

      double functionfac = 0.0;
      if(functnum>0)
      {
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,(ele->Nodes()[i])->X(),time,NULL);
      }
      scnp += functionfac;

      // ----------------------------------------------------
      // convert O2 saturation to O2 concentration
      // ----------------------------------------------------
      // get O2 properties in air
      int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
      // check if O2 properties material exists
      if (id==-1)
      {
        dserror("A material defining O2 properties in air could not be found");
        exit(1);
      }
      const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

      // get atmospheric pressure
      double patm = actmat->atmospheric_p_;
      // get number of O2 moles per unit volume of O2
      double nO2perVO2 = actmat->nO2_per_VO2_;
      // calculate the PO2 at nodes
      double pO2 = scnp*patm;
      // calculate VO2
      double vO2 = volumenp*(pO2/patm);
      // evaluate initial concentration
      scnp = nO2perVO2*vO2/volumenp;
      //------------
      if(i==0)
      {
        int    gid = ele->Id();
        double val = scnp;
        if (vel<0.0)
          val = (*e1scatranp)[ele->LID()];
//        if (ele->Owner()==myrank)
        {
          e1scatranp->ReplaceGlobalValues(1,&val,&gid);
        }
        scatranp(0) = val*area;
      }
      else
      {
        int    gid = ele->Id();
        double val = scnp;
        if (vel>=0.0)
          val = (*e2scatranp)[ele->LID()];
//        if (ele->Owner()==myrank)
        {
          e2scatranp->ReplaceGlobalValues(1,&val,&gid);
        }
        scatranp(1) = val*area;
      }
    }
  }


  {
    scatranp(1) = (*e2scatranp)[ele->LID()]*area;
  }
  if (vel <  0.0)
  {
    scatranp(0) = (*e1scatranp)[ele->LID()]*area;
  }
}//SolveScatra



/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::SolveScatraBifurcations(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    scatranp,
  Epetra_SerialDenseVector&    volumeMix_np,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
//  const int   myrank  = discretization.Comm().MyPID();

  Teuchos::RCP<Epetra_Vector> qin_np   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");

  //  Teuchos::RCP<Epetra_Vector> scatran = params.get<Teuchos::RCP<Epetra_Vector> >("scatran");
  Teuchos::RCP<const Epetra_Vector> scatran  = discretization.GetState("scatranp");

  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);
  // get area
  double area = myvolmix[1];

  // get the elements Qin and Qout
  double q_out = (*qout_np)[ele->LID()];
  double q_in  = (*qin_np)[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> myscatran(lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myscatran,lm);

  // evaluate velocity at nodes (1) and (2)
  double vel1 = q_in /area;
  double vel2 = q_out/area;

  LINALG::Matrix<2,1> velv;
  velv(0,0)=vel1;
  velv(1,0)=vel2;
  // get average velocity
  double vel = 0.5*(vel1+vel2);

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------
  if (vel >=0.0)
  {
    // extrapolate the analytical solution
    double scnp = myscatran[0];
    int gid = ele->Id();
    e1scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
  else
  {
    // extrapolate the analytical solution
    double scnp = myscatran[1];
    int gid = ele->Id();
    e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
}//SolveScatraBifurcations


/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::UpdateScatra(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{

  const int   myrank  = discretization.Comm().MyPID();

  Teuchos::RCP<const Epetra_Vector> dscatranp  = discretization.GetState("dscatranp");
  Teuchos::RCP<Epetra_Vector> dscatranp_m = params.get<Teuchos::RCP<Epetra_Vector> >("dscatranp");
  Teuchos::RCP<Epetra_Vector> qin_np     = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  // get flowrate
  double qin = (*qin_np)[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> mydscatra(lm.size());
  DRT::UTILS::ExtractMyValues(*dscatranp,mydscatra,lm);

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------
  if (qin<0.0)
  {
    int gid = lm[1];
    double val = mydscatra[1];
    if (myrank == ele->Nodes()[1]->Owner())
    {
      dscatranp_m->ReplaceGlobalValues(1,&val,&gid);
    }
  }
}//UpdateScatra


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::UpdateElem12Scatra(
  RedAcinus*                  ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.GetState("scatranp");
  Teuchos::RCP<const Epetra_Vector> dscatranp = discretization.GetState("dscatranp");
  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  Teuchos::RCP<Epetra_Vector> qin_np     = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");
  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> mydscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*dscatranp,mydscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);

  // get flowrate
  double qin = (*qin_np)[ele->LID()];
  // Get the average concentration

  // ---------------------------------------------------------------------
  // element scatra must be updated only at the capillary nodes.
  // ---------------------------------------------------------------------
  //  double e2s = (*e2scatranp)[ele->LID()] + mydscatranp[1]*myvolmix[1];
  double e2s = myscatranp[1];

  int gid = ele->Id();
  e2scatranp->ReplaceGlobalValues(1,&e2s,&gid);
  if (qin<0.0)
  {
    e1scatranp->ReplaceGlobalValues(1,&e2s,&gid);
  }
}



/*----------------------------------------------------------------------*
 |  calculate PO2 from concentration                        ismail 06/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvalPO2FromScatra(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get Po2 vector
  Teuchos::RCP<Epetra_Vector> po2        = params.get<Teuchos::RCP<Epetra_Vector> >("PO2");
  // get Po2 vector
  Teuchos::RCP<const Epetra_Vector> scatran = discretization.GetState("scatranp");
  // get acinar volume
  Teuchos::RCP<Epetra_Vector> acinar_vnp = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  // -------------------------------------------------------------------
  // extract scatra values
  // -------------------------------------------------------------------
  // extract local values from the global vectors
  std::vector<double> myscatran(lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myscatran,lm);

  // -------------------------------------------------------------------
  // find out if the material type is Air or Blood
  // -------------------------------------------------------------------
  std::string fluidType = "none";
  // if RedAirwayScatraAirCond then material type is air
  if (ele->Nodes()[0]->GetCondition("RedAirwayScatraAirCond")!=NULL
    && ele->Nodes()[1]->GetCondition("RedAirwayScatraAirCond")!=NULL)
  {
    fluidType = "air";
  }
  else
  {
    dserror("A scalar transport element must be defined either as \"air\"");
    exit(1);
  }

  // define a empty pO2 vector
  double pO2 = 0.0;

  // -------------------------------------------------------------------
  // Get O2 properties in air
  // -------------------------------------------------------------------
  if (fluidType == "air")
  {
    // -----------------------------------------------------------------
    // Get O2 properties in air
    // -----------------------------------------------------------------

    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
    // check if O2 properties material exists
    if (id==-1)
    {
      dserror("A material defining O2 properties in air could not be found");
      exit(1);
    }
    const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

    // get atmospheric pressure
    double patm = actmat->atmospheric_p_;
    // get number of O2 moles per unit volume of O2
    double nO2perVO2 = actmat->nO2_per_VO2_;

    // -----------------------------------------------------------------
    // Calculate Vo2 in air
    // -----------------------------------------------------------------
    // get airway volume
    double vAir = (*acinar_vnp)[ele->LID()];
    // calculate the VO2 at nodes
    double vO2 = (vAir*myscatran[lm.size()-1])/nO2perVO2;
    // calculate PO2 at nodes
    pO2 = patm*vO2/vAir;
  }
  else
  {
    dserror("A scalar transport element must be defined either as \"air\" or \"blood\"");
    exit(1);
  }

  // -------------------------------------------------------------------
  // Set element pO2 to PO2 vector
  // -------------------------------------------------------------------
  int    gid = lm[lm.size()-1];
  double val = pO2;
  if(myrank == ele->Nodes()[lm.size()-1]->Owner())
  {
    po2->ReplaceGlobalValues(1,&val,&gid);
  }

}// EvalPO2FromScatra


/*----------------------------------------------------------------------*
 |  calculate essential nodal values                        ismail 06/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvalNodalEssentialValues(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    nodal_surface,
  Epetra_SerialDenseVector&    nodal_volume,
  Epetra_SerialDenseVector&    nodal_avg_scatra,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  // ---------------------------------------------------------------------
  // get all general state vectors: flow, pressure,
  // ---------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> acinar_e_v  = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_v");
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.GetState("scatranp");

  // -------------------------------------------------------------------
  // extract scatra values
  // -------------------------------------------------------------------
  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  //----------------------------------------------------------------
  // find the volume of an acinus
  //----------------------------------------------------------------
  // get the current acinar volume
  double volAcinus = (*acinar_e_v)[ele->LID()];
  // set nodal volume
  nodal_volume[1] = volAcinus;

  //----------------------------------------------------------------
  // find the average scalar transport concentration
  //----------------------------------------------------------------
  // set nodal flowrate
  nodal_avg_scatra[0] = myscatranp[1];
  nodal_avg_scatra[1] = myscatranp[1];

  //----------------------------------------------------------------
  // find the total gas exchange surface inside an acinus
  //----------------------------------------------------------------
  // get the initial volume of an acinus
  double volAcinus0;
  ele->getParams("AcinusVolume",volAcinus0);
  // get the initial volume of an alveolar duct
  double volAlvDuct0;
  ele->getParams("AlveolarDuctVolume",volAlvDuct0);
  // find the  number of alveolar duct
  const double numOfAlvDucts = double(floor(volAcinus0/volAlvDuct0));
  // define the number of alveoli per alveolar duct
  const double nAlveoliPerAlveolarDuct = 36.0;
  // define the number of alveoli per duct
  const double nAlveoliPerDuct = 4.0;
  // find the volume of one alveolus
  const double volAlveolus = volAcinus/numOfAlvDucts/nAlveoliPerAlveolarDuct;
  // find the surface of one alveolus
  const double surfAlveolus = (6.0+12.0*sqrt(3.0))*pow(volAlveolus/(8.0*sqrt(2.0)),2.0/3.0);
  // get the length of an edge of an alveolus
  const double al = pow(volAlveolus/(8.0*sqrt(2.0)),1.0/3.0)/3.0;
  // find the surface of an alveolar duct
  const double surfAlveolarDuct = (nAlveoliPerAlveolarDuct - 2.0*nAlveoliPerDuct)*surfAlveolus + 6.0*(al*al);
  // find the surface of an acinus
  const double surfAcinus = surfAlveolarDuct*numOfAlvDucts;
  // set nodal surface area
  nodal_surface[1] = surfAcinus;
}
