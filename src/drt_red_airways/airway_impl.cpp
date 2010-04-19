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

  if (pnp==null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  // extract local values from the global vectors
  vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epn(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epn(i)    = mypn[i];
  }

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         epnp,
         epn,
         elemat1_epetra,
         elevec1_epetra,
         mat,
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
  RCP<vector<int> > lmowner = rcp(new vector<int>);
  ele->LocationVector(discretization,lm,*lmowner);

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
    double val2 = 0.0;
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
  if(myrank == ele->Owner())
  {
    int    gid = ele->Id();
    double val = 0.0;
  }

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
  Epetra_SerialDenseMatrix&                sysmat,
  Epetra_SerialDenseVector&                rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{

  //    cout<<">>>>>>>>>>>>>>>>>Hello Sysmat"<<endl;
  double dens = 0.0;
  double visc = 0.0;

  const int elemVecdim = epnp.Length () ;

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


  // check here, if we really have an airway !!

  // Calculate the length of airway element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));

  double q_out = 0.0;
  //  cout<<"hello Element type"<<endl;
  if(ele->type() == "PoiseuilleResistive")
  {
    // get element information
    double A;
    ele->getParams("Area",A);

    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    const double R = 8.0*PI*visc*dens*L/(pow(A,2));
    sysmat(0,0) = -1.0/R  ; sysmat(0,1) =  1.0/R ;
    sysmat(1,0) =  1.0/R  ; sysmat(1,1) = -1.0/R ;
    
    rhs(0) = 0.0;
    rhs(1) = 0.0;

    q_out = (epnp(0)-epnp(1))/R;

  }
  else if(ele->type() == "InductoResistive")
  {

  }
  else if(ele->type() == "ComplientResistive")
  {
    //    cout<<"Hello ComplientRes"<<endl;
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    const double R = 8.0*PI*visc*dens*L/(pow(A,2));

    // find Capacitance C
    const double C = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));


    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    sysmat(0,0) = -1.0/R -2.0*C/dt ; sysmat(0,1) =  1.0/R ;
    sysmat(1,0) =  1.0/R           ; sysmat(1,1) = -1.0/R ;

    // get element information from the previous time step
    double qcn;
    ele->getVars("capacitor_flow",qcn);

    //------------------------------------------------------------
    //               Calculate the right hand side
    //------------------------------------------------------------
    rhs(0) = -epnp(0)*2.0*C/dt - qcn;
    rhs(1) = 0.0;

    // calculate out flow at the current time step
    q_out = (epnp(0)-epnp(1))/R;
    
    //    cout<<"bye bye RC"<<endl;
  }
  else if(ele->type() == "RLC")
  {
    //    cout<<"RLC"<<endl;
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);
    //    cout<<"A:  "<<A<<endl;
    //    cout<<"Ew: "<<Ew<<endl;
    //    cout<<"tw: "<<tw<<endl;
    // find Resistance 
    const double R = 8.0*PI*visc*dens*L/(pow(A,2));

    // find Capacitance C
    const double C = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // find Inductance I
    const double I = dens*L/A;

    //  cout<<"R: "<<R<<endl;
    //  cout<<"C: "<<C<<endl;
    //  cout<<"L: "<<I<<endl;
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
    ele->getVars("capacitor_flow",qcn);
    ele->getVars("inductor_flow",qln);
    qln =  (epnp(2)-epnp(1))/R;
    //    cout<<"qcn: "<<qcn<<endl;
    //    cout<<"qln: "<<qln<<endl;
    //------------------------------------------------------------
    //               Calculate the right hand side
    //------------------------------------------------------------
    rhs(0) = -qcn - epnp(0)*2.0*C/dt + qln + dt/(2.0*I)*(epnp(0)-epnp(2));
    rhs(1) = 0.0;
    rhs(2) = -qln - dt/(2.0*I)*(epnp(0)-epnp(2));

    //    cout<<sysmat<<endl;
    //    cout<<rhs<<endl;
    // calculate out flow at the current time step
    q_out = (epnp(2)-epnp(1))/R;
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

    DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedLungAcinusCond");
    if(condition)
    {
      if (i==0)
      {
        dserror("SMTHG IS WRONG with Node (%d) in elem(%d)",ele->Nodes()[i]->Id(),ele->Id());
        exit(1);
      }
      string MatType = *(condition->Get<string>("materialType"));

      if (MatType == "NeoHookean")
      {
        const double K = condition->GetDouble("Stiffness1");
        
        sysmat(i,i) += pow(-1.0,i)*(2.0*K/(dt));
        rhs(i)      += pow(-1.0,i)*(q_out + 2.0*K/(dt)*epnp(i));
      }
      else if (MatType == "ViscoElastic_2dof")
      {
        const double E1 = condition->GetDouble("Stiffness1");
        const double E2 = condition->GetDouble("Stiffness2");
        const double B  = condition->GetDouble("Viscosity");

        const double K = E1*E2/2.0 + (E2*B-B*E1)/dt;
        const double Kp= E2/dt + B/pow(dt,2);
        double Qeq_n, Qeq_nm;
        double pnm = epn(i);

        Qeq_nm = -(pnm*B/pow(dt,2))*q_out/K;
        Qeq_n  = (E2*epnp(i)/dt +2.0*epnp(i)/pow(dt,2) + E1*E2/2.0*q_out - (E2*B-B*E1)*q_out)/K;
        
        sysmat(i,i) += pow(-1.0,i)*(Kp/K);
        rhs(i)      += pow(-1.0,i)*(Qeq_n + Qeq_nm);
        
      }
      else
      {
        dserror("[%s] is not defined as a reduced dimensional lung acinus material",MatType.c_str());
        exit(1);
      }
    }
  }
  //    cout<<"bye bye sysmat"<<endl;
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
      
      if(ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
      {
        DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
        // Get the type of prescribed bc
        string Bc = *(condition->Get<string>("boundarycond"));
        
        // double get bc value
        double BCin = 0.0;
        
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
          val = -BCin + (*rhs)[gid];
          rhs->ReplaceGlobalValues(1,&val,&gid);
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
  //  const int numnode = iel;

  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");

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
  for (int i=0;i<lm.size();++i)
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
  double qin = 0.0;
  double qout= 0.0;
  if(ele->type() == "PoiseuilleResistive")
  {
    // get element information
    double A;
    ele->getParams("Area",A);

    const double R = 8.0*PI*visc*dens*L/(pow(A,2));

    qin = (epnp(0)-epnp(1))/R;
    qout= qin;
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
    const double R    = 8.0*PI*visc*dens*L/(pow(A,2));
 
    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));


    // get element information from the previous time step
    double qcn;
    ele->getVars("capacitor_flow",qcn);

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


    qout = (epnp(0)-epnp(1))/R;
    qin  = qout + qcnp_val;
    
  }
  else if(ele->type() == "RLC")
  {
    // get element information
    double A, Ew, tw;
    ele->getParams("Area",A);
    ele->getParams("WallCompliance",Ew);
    ele->getParams("WallThickness",tw);

    // find Resistance 
    const double R    = 8.0*PI*visc*dens*L/(pow(A,2));
 
    // find Capacitance C
    const double C    = 2.0*pow(A,1.5)*L/(Ew*tw*sqrt(M_PI));

    // find Inductance I
    const double I = dens*L/A;

    // get element information from the previous time step
    double qcn,qln;
    ele->getVars("capacitor_flow",qcn);
    

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
    ele->setVars("inductor_flow",qlnp);    

    // -----------------------------------------------------------
    // calculate capacitance flow at time step n+1
    // -----------------------------------------------------------
    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - qcn;
    ele->setVars("capacitor_flow",qcnp_val);
    ele->setVars("inductor_flow",qcnp_val);
    //    int    gid = ele->LID();
    //    double qcnp_val = (2.0*C/dt)*(epnp(0)-epn(0)) - (*qcn)[gid];
    //    qcnp->ReplaceGlobalValues(1,&qcnp_val,&gid);


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

  ele->setVars("flow_in",qin);
  ele->setVars("flow_out",qout);

}

#endif
#endif
