/*----------------------------------------------------------------------*/
/*!
\file airway_impl.cpp

\brief Internal implementation of RedAirway element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/


#ifdef D_RED_AIRWAYS
#ifdef CCADISCRET

#include "airway_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirwayImplInterface* DRT::ELEMENTS::RedAirwayImplInterface::Impl(DRT::ELEMENTS::RedAirway* red_airway)
{
  switch (red_airway->Shape())
  {
  case DRT::Element::line2:
  {
    static AirwayImpl<DRT::Element::line2>* airway;
    if (airway==NULL)
    {
      airway = new AirwayImpl<DRT::Element::line2>;
    }
    return airway;
  }
  default:
    dserror("shape %d (%d nodes) not supported", red_airway->Shape(), red_airway->NumNode());
  }
  return NULL;
}



/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AirwayImpl<distype>::AirwayImpl()
{

}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AirwayImpl<distype>::Evaluate(
  RedAirway*                 ele,
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
  RefCountPtr<const Epetra_Vector> pnm   = discretization.GetState("pnm");

  RefCountPtr<Epetra_Vector> qin_n   = params.get<RefCountPtr<Epetra_Vector> >("qin_n");
  RefCountPtr<Epetra_Vector> qin_np  = params.get<RefCountPtr<Epetra_Vector> >("qin_np");
  RefCountPtr<Epetra_Vector> qout_n  = params.get<RefCountPtr<Epetra_Vector> >("qout_n");
  RefCountPtr<Epetra_Vector> qout_np = params.get<RefCountPtr<Epetra_Vector> >("qout_np");

  if (pnp==null || pn==null || pnm==null )
    dserror("Cannot get state vectors 'pnp', 'pn', an/or 'pnm''");

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
    epnm(i)    = mypnm[i];
  }


  // get the volumetric flow rate from the previous time step
  ParameterList elem_params;
  elem_params.set<double>("qout_n" ,(*qout_n )[ele->LID()]);
  elem_params.set<double>("qout_np",(*qout_np)[ele->LID()]);
  elem_params.set<double>("qin_n"  ,(*qin_n  )[ele->LID()]);
  elem_params.set<double>("qin_np" ,(*qin_np )[ele->LID()]);

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


  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::Initial(
  RedAirway*                             ele,
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

  vector<int>::iterator it = lm.begin();

  //vector<int> lmowner;
  vector<int> lmstride;
  RCP<vector<int> > lmowner = rcp(new vector<int>);
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


    if (gid == 0)
    {
      double val2 = 0.0;
      double A;
      ele->getParams("Area",A);

      val2 = sqrt(A/M_PI);
      radii->ReplaceGlobalValues(1,&val2,&gid);
    }
  }
  if(myrank == (*lmowner)[1])
  {
    int    gid = lm[1];
    double val = 0.0;
    p0np->ReplaceGlobalValues(1,&val,&gid);
    p0n ->ReplaceGlobalValues(1,&val,&gid);
    p0nm->ReplaceGlobalValues(1,&val,&gid);

    double A;
    ele->getParams("Area",A);

    val = sqrt(A/M_PI);
    radii->ReplaceGlobalValues(1,&val,&gid);
  }

  //--------------------------------------------------------------------
  // initialize the volumetric flow rate vectors
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  //  {
  //    int    gid = ele->Id();
  //  double val = 0.0;
  //}

}//AirwayImpl::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::Sysmat(
  RedAirway*                               ele,
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


  double dens = 0.0;
  double visc = 0.0;

  //  const int elemVecdim = epnp.Length () ;

  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    // get actual material
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

    // get density
    dens = actmat->Density();
    
    // get kinetic viscosity
    visc = actmat->Viscosity();
  }
  else
  {
    dserror("Material law is not a Newtonia fluid");
    exit(1);
  }

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

  rhs.Scale(0.0);
  sysmat.Scale(0.0);

#if 0
  if (ele->Owner() != myrank)
  {
    return;
  }
#endif
  // check here, if we really have an airway !!

  // Calculate the length of airway element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));

  double q_out    = params.get<double>("qout_n");
  double qout_np  = params.get<double>("qout_np");

  double q_in    = params.get<double>("qin_n");
  //  double qin_np  = params.get<double>("qin_np");


  if(ele->type() == "PoiseuilleResistive")
  {
    // get element information
    double A;
    ele->getParams("Area",A);

    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    const double Re= 2.0*fabs(qout_np)/(visc*sqrt(A*PI));
    const double R = 8.0*PI*visc*dens*L/(pow(A,2))*(3.4 + 2.1*0.001*Re);

    sysmat(0,0) = -1.0/R  ; sysmat(0,1) =  1.0/R ;
    sysmat(1,0) =  1.0/R  ; sysmat(1,1) = -1.0/R ;

    rhs(0) = 0.0;
    rhs(1) = 0.0;

  }
  else if(ele->type() == "InductoResistive")
  {

  }
  else if(ele->type() == "ComplientResistive")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    const double Re= 2.0*fabs(qout_np)/(visc*sqrt(A*PI));
    const double R = 8.0*PI*visc*dens*L/(pow(A,2))*(3.4 + 2.1*0.001*Re);

    // find Capacitance C
    const double C = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));


    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    sysmat(0,0) = -1.0/R -2.0*C/dt ; sysmat(0,1) =  1.0/R ;
    sysmat(1,0) =  1.0/R           ; sysmat(1,1) = -1.0/R ;

    // get element information from the previous time step
    double qcn = q_in - q_out;
    //  ele->getVars("capacitor_flow",qcn);

    //------------------------------------------------------------
    //               Calculate the right hand side
    //------------------------------------------------------------
    rhs(0) = -epn(0)*2.0*C/dt - qcn;
    rhs(1) = 0.0;

    // calculate out flow at the current time step
    //    q_out = (epnp(0)-epnp(1))/R;
    
  }
  else if(ele->type() == "RLC")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    const double Re= 2.0*fabs(qout_np)/(visc*sqrt(A*PI));
    const double R = 8.0*PI*visc*dens*L/(pow(A,2))*(3.4 + 2.1*0.001*Re);

    // find Capacitance C
    const double C = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // find Inductance I
    const double I = dens*L/A;


    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    //    sysmat(0,0) = -2.0*C/dt -dt/(2.0*I); sysmat(0,1) =  dt/(2.0*I)       ; sysmat(0,2) =  0.0   ;
    //    sysmat(1,0) =  dt/(2.0*I)          ; sysmat(1,1) = -dt/(2.0*I)-1.0/R ; sysmat(1,2) =  1.0/R ;
    //    sysmat(2,0) =  0.0                 ; sysmat(2,1) =  1.0/R            ; sysmat(2,2) = -1.0/R ;


    sysmat(0,0) = -2.0*C/dt -dt/(2.0*I); sysmat(0,1) =  0.0   ; sysmat(0,2) =  dt/(2.0*I)       ;
    sysmat(1,0) =  0.0                 ; sysmat(1,1) = -1.0/R ; sysmat(1,2) =  1.0/R            ;
    sysmat(2,0) =  dt/(2.0*I)          ; sysmat(2,1) =  1.0/R ; sysmat(2,2) = -dt/(2.0*I)-1.0/R ;



    // get element information from the previous time step
    double qcn, qln;
    //  ele->getVars("capacitor_flow",qcn);
    //  ele->getVars("inductor_flow",qln);
    qcn = q_in - q_out;
    qln =  (epn(2)-epn(1))/R;

    //------------------------------------------------------------
    //               Calculate the right hand side
    //------------------------------------------------------------
    rhs(0) = -qcn - epn(0)*2.0*C/dt + qln + dt/(2.0*I)*(epn(0)-epn(2));
    rhs(1) = 0.0;
    rhs(2) = -qln - dt/(2.0*I)*(epn(0)-epn(2));

    // calculate out flow at the current time step
    //    q_out = (epnp(2)-epnp(1))/R;
  }
  else if(ele->type() == "SUKI")
  {

  }
  else
  {
    dserror("[%s] is not an implimented elements yet",(ele->type()).c_str());
    exit(1);
  }

  // -------------------------------------------------------------------
  // Adding the acinus model if Prescribed
  // -------------------------------------------------------------------

  // Check for the simple piston connected to a spring
  for(int i = 0; i<ele->NumNode(); i++)
  {

    //------------------------------------------------------------------
    // Evaluate the number of acini on the end of an outlet.
    // This is hard coded for now, but will be fixed later
    //
    // For now Schroters model is assumed to be Vmodel = 1mm^3
    // Thus we should be given the total-acini-volume/total-t-bronchi-area
    // to calculate the number of Schroters acini at each outlet
    //------------------------------------------------------------------

    DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedLungAcinusCond");
    if(condition)
    {
      double Area  = 0.0;
      ele->getParams("Area",Area);

      //----------------------------------------------------------------
      // Read in the material information
      //----------------------------------------------------------------
      const double VolPerArea = condition->GetDouble("VolumePerArea");
      const double VolAcinus  = condition->GetDouble("Acinus_Volume");
      const double E1 = condition->GetDouble("Stiffness1") / VolAcinus;
      const double E2 = condition->GetDouble("Stiffness2") / VolAcinus;
      const double B  = condition->GetDouble("Viscosity")  / VolAcinus;
      const double NumOfAcini = double(floor(VolPerArea*Area/VolAcinus));
      if (NumOfAcini < 1.0)
      {
        dserror("Acinus condition at node (%d) has zero acini",ele->Nodes()[i]->Id());
      }
      
      // -------------------------------------------------------------
      // Read in the pleural pressure
      // -------------------------------------------------------------
      const  vector<int>*    curve  = condition->Get<vector<int>    >("curve");
      const  vector<double>* vals   = condition->Get<vector<double> >("val");
      
      double Pp_nm = 0.0;
      double Pp_n  = 0.0;
      double Pp_np = 0.0;
      double qn    = q_out;
      
      double pnm = epnm(i);
      double pn  = epn(i);
      
      // evaluate the pleural pressure at (t - dt), (t), and (t + dt)
      if((*curve)[0]>=0)
      {
        Pp_nm = DRT::Problem::Instance()->Curve((*curve)[0]).f(time - dt);
        Pp_nm *= (*vals)[0];
        Pp_n  = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
        Pp_n  *= (*vals)[0];
        Pp_np = DRT::Problem::Instance()->Curve((*curve)[0]).f(time + dt);
        Pp_np *= (*vals)[0];
        
        // Evaluate the pleural contribution
        pn  -= Pp_n;
        pnm -= Pp_nm;
        //        cout<<"pns : "<<pn<<pnm<<p
      }
      
      if (i==0)
      {
        dserror("SMTHG IS WRONG with Node (%d) in elem(%d)",ele->Nodes()[i]->Id(),ele->Id());
        exit(1);
      }
      string MatType = *(condition->Get<string>("materialType"));

      if (MatType == "NeoHookean")
      {
        const double K = condition->GetDouble("Stiffness1");
        
        sysmat(i,i) += pow(-1.0,i)*(2.0/(K*dt));
        rhs(i)      += pow(-1.0,i)*(q_out + 2.0/(K*dt)*epn(i));
      }
      else if (MatType == "KelvinVoigt")
      {
        const double R  = B;
        const double Kp_np = 2.0/(E1*dt) + 1.0/R;
        const double Kp_n  = 2.0/(E1*dt) - 1.0/R;

        //        cout<<"KS: "<<Kp_np<<"  "<<Kp_n<<"  "<<R<<endl;
        sysmat(i,i) += pow(-1.0,i)*(Kp_np*NumOfAcini);
        //        cout<<"SYSmat: "<<pow(-1.0,i)*(Kp_np*NumOfAcini)<<endl;
        rhs(i)      += pow(-1.0,i)*(q_out + pn*NumOfAcini * Kp_n + Pp_np*NumOfAcini * Kp_np);
        //        cout<<"rhs: "<<pow(-1.0,i)*(q_out + epn(i)*NumOfAcini * Kp_n + Pp_np*NumOfAcini * Kp_np)<<endl;
      }
      else if (MatType == "ViscoElastic_2dof")
      {

#if 0
        const double K = E1*E2/2.0 + (E2*B+B*E1)/dt;
        const double Kp= E2/dt + B/pow(dt,2);
        double Qeq_n, Qeq_nm;
        double pnm = epnm(i);

        Qeq_nm = -(pnm*B/pow(dt,2))/K;
        Qeq_n  = (E2*epn(i)/dt +2.0*B*epn(i)/pow(dt,2) + E1*E2/2.0*q_out - (E2*B+B*E1)*q_out/dt)/K;
        
        sysmat(i,i) += pow(-1.0,i)*(Kp/K);
        rhs(i)      += pow(-1.0,i)*(Qeq_n + Qeq_nm);
#endif
        const double R  = B;
        const double Kp_np =  (E1+E2)/dt + R/(dt*dt);
        const double Kp_n  = -(E1+E2)/dt - 2.0*R/(dt*dt);
        const double Kp_nm =  (R/(dt*dt));
        const double Kq_np =  (E1*E2 + R*E1/dt);
        const double Kq_n  =  (-R*E1/dt);

        double Qeq_n, Qeq_nm;

        Qeq_nm = (-Kp_nm* pnm *NumOfAcini)/Kq_np;
        Qeq_n  = (-Kp_n * pn  *NumOfAcini + Kq_n*qn)/Kq_np;
        
        // -------------------------------------------------------------
        // Evaluate sysmat and rhs
        // -------------------------------------------------------------
        sysmat(i,i) += pow(-1.0,i)*(Kp_np/Kq_np)*NumOfAcini;
        rhs(i)      += pow(-1.0,i)*(Qeq_n + Qeq_nm + Kp_np/Kq_np*Pp_np*NumOfAcini);
      }
      else
      {
        dserror("[%s] is not defined as a reduced dimensional lung acinus material",MatType.c_str());
        exit(1);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::EvaluateTerminalBC(
  RedAirway*                   ele,
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
      if(ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond") || ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
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

          cout<<"Return ["<<Bc<<"] form 3D problem to 1D POINT of ID["<<ID<<"]: "<<BCin<<endl;
        }
        else
        {

        }

        if (Bc == "pressure")
        {
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
          dserror("precribed [%s] is not defined for reduced airways",Bc.c_str());
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

          // const double* X = ele->Nodes()[i]->X();
          // printf("WARNING: node(%d) is free on [%f,%f,%f] \n",gid+1,X[0],X[1],X[2]);
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
void DRT::ELEMENTS::AirwayImpl<distype>::CalcFlowRates(
  RedAirway*                   ele,
  ParameterList&               params,
  DRT::Discretization&         discretization,
  vector<int>&                 lm,
  RefCountPtr<MAT::Material>   material)
{
  //  const int   myrank  = discretization.Comm().MyPID();

  //  const int numnode = iel;

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  RefCountPtr<Epetra_Vector> qin_np  = params.get<RefCountPtr<Epetra_Vector> >("qin_np");
  RefCountPtr<Epetra_Vector> qout_np = params.get<RefCountPtr<Epetra_Vector> >("qout_np");
  RefCountPtr<Epetra_Vector> qin_n  = params.get<RefCountPtr<Epetra_Vector> >("qin_n");
  RefCountPtr<Epetra_Vector> qout_n = params.get<RefCountPtr<Epetra_Vector> >("qout_n");


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
    
    // get kinetic viscosity
    visc = actmat->Viscosity();
  }
  else
  {
    dserror("Material law is not a Newtonia fluid");
    exit(1);
  }



  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  const int numnode = lm.size();

  // create objects for element arrays
  //  LINALG::Matrix<numnode,1> epnp;
  Epetra_SerialDenseVector epnp(numnode);
  for (unsigned int i=0;i<lm.size();++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
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

  // check here, if we really have an airway !!

  // Calculate the length of airway element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));


  //--------------------------------------------------------------
  //               Calculate flowrates
  //--------------------------------------------------------------
  double qin = (*qin_n )[ele->LID()];
  double qout= (*qout_n)[ele->LID()];
  double Qin_np  = (*qin_np )[ele->LID()];
  double Qout_np = (*qout_np)[ele->LID()];
  if(ele->type() == "PoiseuilleResistive")
  {
    // get element information
    double A;
    ele->getParams("Area",A);

    if(params.get<string> ("solver type") == "Nonlinear")
    {
      const int maxitr= 10;
      for (int itr =0; itr<=maxitr; itr++)
      {
        //        const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
        const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
        const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
        const double R = K1 + K2*fabs(Qout_np);
        
        const double df_dq = K1 + 2.0*K2*fabs(Qout_np);
        const double     f = - (epnp(0)-epnp(1)) + R*Qout_np;
        qin = -f/df_dq + Qout_np;
        //      qin = (epnp(0)-epnp(1))/R;
        
        if (fabs(Qin_np - qin) <= fabs(1e-6*qin))
        {
          Qin_np = Qout_np = qout= qin;
          break;
        }
        if (itr == maxitr)
        {
          cout<<"Warning didn't converge: "<< fabs(Qin_np - qin) <<" Q = "<<qout<<endl;
        }
        Qin_np = Qout_np = qout= qin;
      }     
    }
    if(params.get<string> ("solver type") == "Linear")
    {
      //      const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
      const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
      const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
      const double R = K1 + K2*fabs(Qout_np);

      qin = (epnp(0)-epnp(1))/R;
      Qin_np = Qout_np = qout= qin;
    }

  }
  else if(ele->type() == "InductoResistive")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

  }
  else if(ele->type() == "ComplientResistive")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    //    const double R    = 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
 
    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));


    // get element information from the previous time step
    double qcn =  qin - qout;
    //    ele->getVars("capacitor_flow",qcn);

    RefCountPtr<const Epetra_Vector> pn   = discretization.GetState("pn");
    //    RefCountPtr<const Epetra_Vector> qcn  = discretization.GetState("qcn");
    //    RCP<Epetra_Vector>               qcnp = params.get<RCP<Epetra_Vector> >("qcnp");


    // extract local values from the global vectors
    vector<double> mypn(lm.size());
    DRT::UTILS::ExtractMyValues(*pn,mypn,lm);
    
    // create objects for element arrays
    LINALG::Matrix<iel,1> epn;
    for (int i=0;i<numnode;++i)
    {
      // split area and volumetric flow rate, insert into element arrays
      epn(i)    = mypn[i];
    }

    // -----------------------------------------------------------
    // calculate capacitance flow at time step n+1
    // -----------------------------------------------------------
    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - qcn;
    ele->setVars("capacitor_flow",qcnp_val);
    //    int    gid = ele->LID();
    //    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - (*qcn)[gid];
    //    qcnp->ReplaceGlobalValues(1,&qcnp_val,&gid);


    if(params.get<string> ("solver type") == "Nonlinear")
    {
      const int maxitr= 10;
      for (int itr =0; itr<=maxitr; itr++)
      {
        //        const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
        const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
        const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
        const double R = K1 + K2*fabs(Qout_np);
        
        const double df_dq = K1 + 2.0*K2*fabs(Qout_np);
        const double     f = - (epnp(0)-epnp(1)) + R*Qout_np;
        qout = -f/df_dq + Qout_np;
        //      qin = (epnp(0)-epnp(1))/R;
      
        if (fabs(Qout_np - qout) <= fabs(1e-6*qout))
        {
          Qout_np = qout;
          break;
        }
        if (itr == maxitr)
        {
          cout<<"Warning didn't converge: "<< fabs(Qout_np - qout) <<" Q = "<<qout<<endl;
        }
        Qout_np = qout;
      }     
    }
    if(params.get<string> ("solver type") == "Linear")
    {
      //      const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
      const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
      const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
      const double R = K1 + K2*fabs(Qout_np);

      qout = (epnp(0)-epnp(1))/R;
      Qout_np = qout;
    }

   Qin_np = qin  = qout + qcnp_val;
    
  }
  else if(ele->type() == "RLC")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    const double R    = 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
 
    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // find Inductance I
    //    const double I = dens*L/A;

    // get element information from the previous time step
    double qcn,qln;
    qcn = qin - qout;
    //    ele->getVars("capacitor_flow",qcn);
    

    // get pressure values at n+1
    RefCountPtr<const Epetra_Vector> pn   = discretization.GetState("pn");

    // extract local values from the global vectors
    vector<double> mypn(lm.size());
    DRT::UTILS::ExtractMyValues(*pn,mypn,lm);
    
    // create objects for element arrays
    Epetra_SerialDenseVector epn(numnode);
    //  LINALG::Matrix<numnode,1> epn;
    for (int i=0;i<numnode;++i)
    {
      // split area and volumetric flow rate, insert into element arrays
      epn(i)    = mypn[i];
    }

    // -----------------------------------------------------------
    // calculate the inductance flow at time step n and n+1
    // qln = qRn (qRn: flow in resistor at n)
    // -----------------------------------------------------------
    double qlnp = 0.0;
    qln  = ( epn(2) -  epn(1))/R;
    qlnp = (epnp(2) - epnp(1))/R;//((epnp(0) - epnp(2))+(epn(0) - epn(2)))*dt/(2.0*I) + qln;

    // -----------------------------------------------------------
    // calculate capacitance flow at time step n+1
    // -----------------------------------------------------------
    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - qcn;

    //    int    gid = ele->LID();
    //    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - (*qcn)[gid];
    //    qcnp->ReplaceGlobalValues(1,&qcnp_val,&gid);

    if(params.get<string> ("solver type") == "Nonlinear")
    {
      const int maxitr= 10;
      for (int itr =0; itr<=maxitr; itr++)
      {
        //        const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
        const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
        const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
        const double R = K1 + K2*fabs(Qout_np);
        
        const double df_dq = K1 + 2.0*K2*fabs(Qout_np);
        const double     f = - (epnp(0)-epnp(1)) + R*Qout_np;
        qout = -f/df_dq + Qout_np;

  
        if (fabs(Qout_np - qout) <= fabs(1e-6*qout))
        {
          Qout_np = qout;
          break;
        }
        if (itr == maxitr)
        {
          cout<<"Warning didn't converge: "<< fabs(Qout_np - qout) <<" Q = "<<qout<<endl;
        }
        Qout_np = qout;
      }     
    }
    if(params.get<string> ("solver type") == "Linear")
    {
      //      const double Re= 2.0*fabs(Qout_np)/(visc*sqrt(A*PI))*1.0;
      const double K1= 8.0*PI*visc*dens*L/(pow(A,2))*3.4;
      const double K2= 8.0*PI*visc*dens*L/(pow(A,2))*2.1*0.001*2.0/(visc*sqrt(A*PI));
      const double R = K1 + K2*fabs(Qout_np);

      qout = (epnp(0)-epnp(1))/R;
      Qout_np = qout;
    }

    Qin_np = qin  = qout + qcnp_val;
   
    qout = qlnp;
    qin  = qout + qcnp_val;
  }
  else if(ele->type() == "SUKI")
  {
    
  }
  else
  {
    dserror("[%s] is not an implimented elements yet",(ele->type()).c_str());
    exit(1);
  }

  //  ele->setVars("flow_in",qin);
  //  ele->setVars("flow_out",qout);
  int gid = ele->Id();

  qin_np  -> ReplaceGlobalValues(1,&qin,&gid);
  qout_np -> ReplaceGlobalValues(1,&qout,&gid);


}//CalcFlowRates


/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::GetCoupledValues(
  RedAirway*                   ele,
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
          if (i==0)
            ele->getVars("flow_in", BC3d);
          else
            ele->getVars("flow_out", BC3d);
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
        
        cout<<"Return ["<<returnedBC<<"] form 1D problem to 3D SURFACE of ID["<<ID<<"]: "<<BC3d<<endl;

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

#endif
#endif
