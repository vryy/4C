/*----------------------------------------------------------------------*/
/*!
\file airway_impl.cpp

\brief Internal implementation of artery_lin_exp element

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

#include "../drt_mat/cnst_1d_art.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AirwayImplInterface* DRT::ELEMENTS::AirwayImplInterface::Impl(DRT::ELEMENTS::RedAirway* red_airway)
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
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AirwayImpl<distype>::AirwayImpl()
  : qn_(),
    pn_()
{

}

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

  // the number of nodes
  const int numnode = iel;
  vector<int>::iterator it_vcr;

  // construct views
  LINALG::Matrix<1*iel,1*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<1*iel,    1> elevec1(elevec1_epetra.A(),true);
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

  RefCountPtr<const Epetra_Vector> qnp  = discretization.GetState("qnp");
  RefCountPtr<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (qnp==null)
    dserror("Cannot get state vectors 'qnp'");
  if (pnp==null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  vector<double> myqnp(lm.size());
  DRT::UTILS::ExtractMyValues(*qnp,myqnp,lm);

  vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> epnp;
  LINALG::Matrix<numnode,1> eqnp;
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqnp[i];
    epnp(i)    = mypnp[i];
  }
  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         eqnp,
         epnp,
         elemat1,
         elevec1,
         mat,
         dt);


  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 07/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::Initial(
  RedAirway*                             ele,
  ParameterList&                         params,
  DRT::Discretization&                   discretization,
  vector<int>&                           lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  RCP<Epetra_Vector> q0    = params.get<RCP<Epetra_Vector> >("q0");
  RCP<Epetra_Vector> p0    = params.get<RCP<Epetra_Vector> >("p0");
  vector<int>        lmown = *(params.get<RCP<vector<int> > >("lmowner"));
  int myrank  = discretization.Comm().MyPID();

  vector<int>::iterator it = lm.begin();

  if(myrank == lmown[0])
  {
    int gid = lm[0];
    double val = 0.0;
    q0->ReplaceGlobalValues(1,&val,&gid);
    p0->ReplaceGlobalValues(1,&val,&gid);
  }
  if(myrank == lmown[1])
  {
    int gid = lm[1];
    double val = 0.0;
    q0->ReplaceGlobalValues(1,&val,&gid);
    p0->ReplaceGlobalValues(1,&val,&gid);
  }


}//ArteryLinExp::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 07/09|
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                               ______                                 |
 |                     _____-----      -----_____                       |
 |           _______---                          ---______              |
 | ->       / \                                         / \   ->        |
 | -->     |   |                                       |   |  -->       |
 | ---->  |     |                                     |     | ---->     |
 | ---->  |     |                                     |     | ---->     |
 | ---->  |     |                                     |     | ---->     |
 | -->     |   |                                       |   |  -->       |
 | ->       \_/_____                                ____\_/   ->        |
 |                  ---_____                _____---                    |
 |                          -----______-----                            |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AirwayImpl<distype>::Sysmat(
  RedAirway*                               ele,
  const LINALG::Matrix<iel,1>&             eqnp,
  const LINALG::Matrix<iel,1>&             epnp,
  LINALG::Matrix<1*iel,1*iel>&             sysmat,
  LINALG::Matrix<1*iel,    1>&             rhs,
  Teuchos::RCP<const MAT::Material>        material,
  double                                   dt)
{

  LINALG::Matrix<1*iel,1> qn;
  LINALG::Matrix<1*iel,1> pn;
  for(int i=0; i<iel; i++)
  {
    qn(i,0) = eqnp(i);
    pn(i,0)   = epnp(i);

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

  rhs.Clear();
  sysmat.Clear();


  // check here, if we really have an airway !!

  // Calculate the length of artery element
  const double L=sqrt(
            pow(xyze(0,0) - xyze(0,1),2)
          + pow(xyze(1,0) - xyze(1,1),2)
          + pow(xyze(2,0) - xyze(2,1),2));

  if(ele->type() == "PoiseuilleResistive")
  {
    //------------------------------------------------------------
    //               Calculate the System Matrix
    //------------------------------------------------------------
    const double R = 1.0;
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
  }
  else if(ele->type() == "RLC")
  {
  }
  else if(ele->type() == "SUKI")
  {

  }
  else
  {
    dserror("[%s] is not an implimented elements yet",(ele->type()).c_str());
    exit(1);
  }


#if 0
  cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<endl;
  cout<<"+------------------------ THE FINAL R-LHS------------------------+"<<endl;
  cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<endl;
  cout<<"rhs is: "<<rhs<<endl;
  cout<<"lhs is: "<<sysmat<<endl;
  cout<<"With L= "<<L<<endl;
  cout<<"|+++++++++++++++++++++++++!!!!      !!!!-------------------------|"<<endl;
  cout<<"+-------------------------!!!!!!!!!!!!!!-------------------------+"<<endl;
#endif

}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 07/09|
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

#if 0
  RCP<Epetra_Vector>   Wfnp  = params.get<RCP<Epetra_Vector> >("Wfnp");
  RCP<Epetra_Vector>   Wbnp  = params.get<RCP<Epetra_Vector> >("Wbnp");

  // get time-step size
  const double dt = params.get<double>("time step size");

  // check here, if we really have an artery !!
  // Define Geometric variables
  double Ao1 = 0.0;
  double Ao2 = 0.0;
  // Define blood material variables
  double visc=0.0;
  double dens=0.0;
  // Define artery's material variables
  double t1 = 0.0;
  double t2 = 0.0;
  double E1 = 0.0;
  double E2 = 0.0;
  double nue= 0.0;
  // Define artery's external forces
  double pext1 =0.0;
  double pext2 =0.0;
  // check here, if we really have an artery !!
  if( material->MaterialType() == INPAR::MAT::m_cnst_art)
  {
    const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());
    // Read in initial cross-sectional area at node 1
    Ao1    = M_PI*pow(actmat->Diam()/2,2);
    // Read in initial cross-sectional area at node 2
    Ao2    = Ao1;
    // Read in blood density
    dens   = actmat->Density();
    // Read in blodd viscosity
    visc   = actmat->Viscosity();
    // Read in artery's thickness at node 1
    t1     = actmat->Th();
    // Read in artery's thickness at node 2
    t2     = t1;
    // Read in artery's Youngs modulus of elasticity thickness at node 1
    E1     = actmat->Young();
    // Read in artery's Youngs modulus of elasticity thickness at node 2
    E2     = E1;
    // Read in artery's Poisson's ratio
    nue    = actmat->Nue();
    // Read in artery's external forces at node 1
    pext1  = actmat->pext(0);
    // Read in artery's external forces at node 1
    pext2  = actmat->pext(1);
    
    // Set up all the needed vectors for furthur calculations
    area0_(0,0) = Ao1;
    area0_(1,0) = Ao2;
    th_(0,0)    = t1;
    th_(1,0)    = t2;
    young_(0,0) = E1;
    young_(1,0) = E2;
    pext_(0,0)  = pext1;
    pext_(1,0)  = pext2;
  }
  else
    dserror("Material law is not an artery");


  // the number of nodes
  const int numnode = iel;
  vector<int>::iterator it_vcr;

  RefCountPtr<const Epetra_Vector> qanp  = discretization.GetState("qanp");

  if (qanp==null)
    dserror("Cannot get state vectors 'qanp'");

  // extract local values from the global vectors
  vector<double> myqanp(lm.size());
  DRT::UTILS::ExtractMyValues(*qanp,myqanp,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eareanp;
  LINALG::Matrix<numnode,1> eqnp;

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    eqnp(i)    = myqanp[1+(i*2)];
    qn_(i)     = myqanp[1+(i*2)];
    eareanp(i) = myqanp[0+(i*2)];
    an_(i)     = myqanp[0+(i*2)];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BC at the inlet and outlet
  // ---------------------------------------------------------------------------------
  // IO_BC_HERE
  for (int i =0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtInOutCond"))
    {
      double TermIO  = 0.0;
      // Get the in/out terminal condition
      string TerminalType = *(ele->Nodes()[i]->GetCondition("ArtInOutCond")->Get<string>("terminaltype"));
      if(TerminalType=="inlet")
        TermIO = -1.0;
      else if (TerminalType=="outlet")
        TermIO = 1.0;
      else
        dserror("Something is severely wrong! In/Out terminal condition should be either \"outlet\" or \"inlet\"");

      RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
      RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

      if (bcval==null||dbctog==null)
        dserror("Cannot get state vectors 'bcval' and 'dbctog'");


      th_(0,0)    = t1;
      th_(1,0)    = t2;
      young_(0,0) = E1;
      young_(1,0) = E2;
      const double beta = sqrt(PI)*young_(i)*th_(i)/(1.0-nue*nue);
      double Wf, Wb;

      // -----------------------------------------------------------------------------
      // fill the required parameters to solve the inlet BC
      // -----------------------------------------------------------------------------
      ParameterList Cparams;
      Cparams.set<int>   ("in out flag",int(TermIO));
      Cparams.set<double>("total time", params.get<double>("total time"));
      Cparams.set<double>("artery beta",beta);
      Cparams.set<double>("artery area",area0_(i));
      Cparams.set<double>("blood density",dens);
      int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
      if (local_id< 0 )
      {
        dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
        exit(1);
      }
                            
      if (TermIO == -1)
      {
         Cparams.set<double>("backward characteristic wave speed",(*Wbnp)[local_id]);
      }
      else
      {
         Cparams.set<double>("forward characteristic wave speed",(*Wfnp)[local_id]);
      }
      Cparams.set<double>("external pressure",pext_(i));

      // -----------------------------------------------------------------------------
      // Solve any possible prescribed boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtPrescribedCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtPrescribedCond");
        Cparams.set<string>("Condition Name","ArtPrescribedCond");
        ART::UTILS::SolvePrescribedTerminalBC(rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible 3-D/reduced-D coupled boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
      {
        RCP<ParameterList> CoupledTo3DParams = params.get<RCP<ParameterList> >("coupling with 3D fluid params");
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");
        Cparams.set<RCP<ParameterList > >("coupling with 3D fluid params",CoupledTo3DParams);
        Cparams.set<string>("Condition Name","Art_redD_3D_CouplingCond");

        ART::UTILS::SolvePrescribedTerminalBC(rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible reflection boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtRfCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtRfCond");
        ART::UTILS::SolveReflectiveTerminal(rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // Solve any possible windkessel boundary condition
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtWkCond"))
      {
        const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("ArtWkCond");
        Cparams.set<double>("time step size",dt);
        Cparams.set<double>("external pressure",pext_(i));
        Cparams.set<double>("terminal volumetric flow rate",qn_(i));
        Cparams.set<double>("terminal cross-sectional area",an_(i));
        ART::UTILS::SolveExplWindkesselBC(rcp(&discretization,false), condition, Cparams);
      }

      // -----------------------------------------------------------------------------
      // break the for loopIf the boundary condition is a junction, since it will be solved later
      // -----------------------------------------------------------------------------
      if(ele->Nodes()[i]->GetCondition("ArtJunctionCond")==NULL)
      {
        Wf = Cparams.get<double>("forward characteristic wave speed");
        Wb = Cparams.get<double>("backward characteristic wave speed");

        // -----------------------------------------------------------------------------
        // Modify the global forward and backward characteristics speeds vector
        // -----------------------------------------------------------------------------
        int myrank  = discretization.Comm().MyPID();
        if(myrank == ele->Nodes()[i]->Owner())
        {
          int    gid  = ele->Nodes()[i]->Id();
          if (TermIO == -1.0)
          {
            double val1 = Wf;
            Wfnp->ReplaceGlobalValues(1,&val1,&gid);
          }
          else
          {
            double val2 = Wb;
            Wbnp->ReplaceGlobalValues(1,&val2,&gid);
            int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
            Wf = (*Wfnp)[local_id];
          }
        }
        
        // calculating A at node i
        int    gid; 
        double val; 
        double cross_area;

        cross_area = pow(2.0*dens*area0_(i)/beta,2)*pow((Wf - Wb)/8.0,4);
        cout<<"CROSS_SECTIONAL AREA IS: "<<cross_area<<endl;
        
        gid = lm[2*i];
        val = cross_area;
        bcval->ReplaceGlobalValues(1,&val,&gid);

        gid = lm[2*i];
        val = 1;
        dbctog->ReplaceGlobalValues(1,&val,&gid);

        // calculating Q at node i
        gid = lm[2*i+1];
        val = (cross_area)*(Wf + Wb)/2.0;
        bcval->ReplaceGlobalValues(1,&val,&gid);

        gid = lm[2*i+1];
        val = 1;
        dbctog->ReplaceGlobalValues(1,&val,&gid);
      }
    } // End of node i has a condition
  } //End of for loop

  // ---------------------------------------------------------------------------------
  // Solve the any available junction boundary conditions
  // ---------------------------------------------------------------------------------   
  for (int i =0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("ArtJunctionCond"))
    {
      RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > junc_nodal_vals;
      junc_nodal_vals = params.get<RCP<map<const int, RCP<ART::UTILS::JunctionNodeParams> > > >("Junctions Parameters");

      RefCountPtr<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
      RefCountPtr<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

      // -------------------------------------------------------------------------------
      // Update the Dirichlet BC vector
      // -------------------------------------------------------------------------------
      int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
      int    gid;
      double val;
      // set A at node i
      gid = lm[2*i  ];
      val = (*junc_nodal_vals)[local_id]->A_;
      bcval->ReplaceGlobalValues(1,&val,&gid);

      val = 1;
      dbctog->ReplaceGlobalValues(1,&val,&gid);

      // set Q at node i
      gid = lm[2*i+1];
      val = (*junc_nodal_vals)[local_id]->Q_;
      bcval->ReplaceGlobalValues(1,&val,&gid);

      val = 1;
      dbctog->ReplaceGlobalValues(1,&val,&gid);
    }
  }
#endif 
}

#endif
#endif
