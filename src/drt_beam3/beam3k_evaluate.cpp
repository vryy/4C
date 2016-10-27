/*----------------------------------------------------------------------*/
/*!
\file beam3k_evaluate.cpp

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

\maintainer Christoph Meier
*/
/*-----------------------------------------------------------------------------------------------------------*/


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_structure_new/str_timint_basedatasdyn.H"


#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "beam3k.H"


/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 01/16|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3k::Evaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        std::vector<int>& lm,
                                        Epetra_SerialDenseMatrix& elemat1, //stiffness matrix
                                        Epetra_SerialDenseMatrix& elemat2, //mass matrix
                                        Epetra_SerialDenseVector& elevec1, //internal forces
                                        Epetra_SerialDenseVector& elevec2, //inertia forces
                                        Epetra_SerialDenseVector& elevec3)
{
  SetParamsInterfacePtr(params);

  // Set statmech params interface pointer
  if (IsParamsInterface())
    SetStatMechParamsInterfacePtr();

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
   act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if     (action == "calc_none")         dserror("No action supplied");
    else if (action=="calc_struct_linstiff")                               act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")                               act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce")                          act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")                           act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")                           act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass")                          act = ELEMENTS::struct_calc_nlnstifflmass; //with lumped mass matrix
    else if (action=="calc_struct_stress")                                 act = ELEMENTS::struct_calc_stress;
    else if (action=="calc_struct_eleload")                                act = ELEMENTS::struct_calc_eleload;
    else if (action=="calc_struct_fsiload")                                act = ELEMENTS::struct_calc_fsiload;
    else if (action=="calc_struct_update_istep")                           act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")                            act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_ptcstiff")                               act = ELEMENTS::struct_calc_ptcstiff;
    else if (action=="calc_struct_energy")                                 act = ELEMENTS::struct_calc_energy;
    else     dserror("Unknown type of action for Beam3k");
  }

  switch(act)
  {
    case ELEMENTS::struct_calc_ptcstiff:
    {
      dserror("no ptc implemented for Beam3k element");
    }
    break;

    case ELEMENTS::struct_calc_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    case ELEMENTS::struct_calc_internalinertiaforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");

      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);     // ToDo unused?

      // get element velocity vector
      std::vector<double> myvel(lm.size());

      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      if (act == ELEMENTS::struct_calc_nlnstiffmass)
      {
        if(weakkirchhoff_)
          CalculateInternalForcesWK(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
        else
          CalculateInternalForcesSK(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
      }
      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        dserror("The action calc_struct_nlnstifflmass is not implemented yet!");
      }
      else if (act == ELEMENTS::struct_calc_nlnstiff)
      {
        if(weakkirchhoff_)
          CalculateInternalForcesWK(params,mydisp,&elemat1,NULL,&elevec1,NULL);
        else
          CalculateInternalForcesSK(params,mydisp,&elemat1,NULL,&elevec1,NULL);
      }
      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        if(weakkirchhoff_)
          CalculateInternalForcesWK(params,mydisp,NULL,NULL,&elevec1,NULL);
        else
          CalculateInternalForcesSK(params,mydisp,NULL,NULL,&elevec1,NULL);
      }
      else if (act == ELEMENTS::struct_calc_internalinertiaforce)
      {
        if(weakkirchhoff_)
          CalculateInternalForcesWK(params,mydisp,NULL,NULL,&elevec1,&elevec2);
        else
          CalculateInternalForcesSK(params,mydisp,NULL,NULL,&elevec1,&elevec2);
      }

      //ATTENTION: In order to perform a brief finite difference check of the nonlinear stiffness matrix the code block
      //"FD-CHECK" from the end of this file has to be copied to this place:
      //***************************************************************************************************************
                                                 //Insert code block here!
      //***************************************************************************************************************

    }
    break;

    case ELEMENTS::struct_calc_brownianforce:
    case ELEMENTS::struct_calc_brownianstiff:
    {
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get element velocity
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      if (act == ELEMENTS::struct_calc_brownianforce)
        CalcBrownianForcesAndStiff<2,2,3>(params,myvel,mydisp,NULL,&elevec1);
      else if (act == ELEMENTS::struct_calc_brownianstiff)
        CalcBrownianForcesAndStiff<2,2,3>(params,myvel,mydisp,&elemat1,&elevec1);
      else
        dserror("You shouldn't be here.");

      break;
    }

    case ELEMENTS::struct_calc_energy:
    {
      elevec1(0)=Eint_;
      //elevec1(1)=Ekin_;
    }
    break;

    case ELEMENTS::struct_calc_stress:
    {
      dserror("No stress output implemented for beam3k elements");
    }
    break;

    case ELEMENTS::struct_calc_update_istep:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam at the end of the time step has to be stored*/
      Qconvmass_ = Qnewmass_;
      wconvmass_ = wnewmass_;
      aconvmass_ = anewmass_;
      amodconvmass_ = amodnewmass_;
      rttconvmass_ = rttnewmass_;
      rttmodconvmass_ = rttmodnewmass_;
      rtconvmass_ = rtnewmass_;
      rconvmass_ = rnewmass_;

      Qrefconv_ = Qrefnew_;
    }
    break;

    case ELEMENTS::struct_calc_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitable timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/

      Qnewmass_=Qconvmass_;
      wnewmass_=wconvmass_;
      anewmass_=aconvmass_;
      amodnewmass_=amodconvmass_;
      rttnewmass_=rttconvmass_;
      rttmodnewmass_=rttmodconvmass_;
      rtnewmass_=rtconvmass_;
      rnewmass_=rconvmass_;

      Qrefnew_ = Qrefconv_;
    }

    case ELEMENTS::struct_calc_recover:
    {
      // do nothing here
      break;
    }

    default:
      std::cout << "\ncalled element with action type " << ActionType2String(act);
      dserror("This action type is not implemented for Beam3k");
     break;

  }//switch(act)

  return 0;

}  //DRT::ELEMENTS::Beam3k::Evaluate

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                          meier 01/16|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::Lumpmass(Epetra_SerialDenseMatrix* emass)
{
  dserror("Lumped mass matrix not implemented yet!!!");
}

/*------------------------------------------------------------------------------------------------------------*
 | internal forces: weak Kirchhoff constraint (private)                                            meier 01/16|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::CalculateInternalForcesWK( Teuchos::ParameterList& params,
                                                        std::vector<double>&      disp,
                                                        Epetra_SerialDenseMatrix* stiffmatrix,
                                                        Epetra_SerialDenseMatrix* massmatrix,
                                                        Epetra_SerialDenseVector* force,
                                                        Epetra_SerialDenseVector* inertia_force)
{
  TEUCHOS_FUNC_TIME_MONITOR("Beam3k::CalculateInternalForcesWK");
  if(BEAM3K_COLLOCATION_POINTS!=2 and BEAM3K_COLLOCATION_POINTS!=3 and BEAM3K_COLLOCATION_POINTS!=4)
    dserror("Only the values 2,3 and 4 are valid for BEAM3K_COLLOCATION_POINTS!!!");

  //number of nodes fixed for these element
  const int nnode = 2;

  if(disp.size()!=(6*nnode+BEAM3K_COLLOCATION_POINTS))
    dserror("Number of BEAM3K_COLLOCATION_POINTS does not match number of nodes defined in the input file!!!");

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int_aux(true);

  //vector for current nodal DoFs in total Lagrangian style, i.e. displacement + initial values:
  //rotvec_==true:  disp_totlag=[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3]
  //rotvec_==false: disp_totlag=[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3]
  //! The Number of collocation points can take on the values 2, 3 and 4. 3 and 4 are interior nodes.
  //This leads e.g. in the case rotvec_==true to the following ordering:
  //if BEAM3K_COLLOCATION_POINTS = 2: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2]
  //if BEAM3K_COLLOCATION_POINTS = 3: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3]
  //if BEAM3K_COLLOCATION_POINTS = 4: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3 \alpha_4]
  std::vector<FAD> disp_totlag(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);

  //vector containing locally assembled nodal positions and tangents required for centerline:
  //r(s)=N(s)*disp_totlag_centerline, with disp_totlag_centerline=[\v{d}_1, \v{t}_1, 0, \v{d}_2, \v{t}_2, 0, 0]
  std::vector<FAD> disp_totlag_centerline(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);

  //CP values of strains and their variations needed for interpolation
  std::vector<FAD> epsilon_cp(BEAM3K_COLLOCATION_POINTS); // axial tension epsilon=|r_s|-1 at collocation points
  std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> > v_epsilon_cp(BEAM3K_COLLOCATION_POINTS);
  std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> > v_thetaperp_cp(BEAM3K_COLLOCATION_POINTS);
  std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> > v_thetapar_cp(BEAM3K_COLLOCATION_POINTS);

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  //interpolated values of strains and their variations evaluated at Gauss points
  FAD epsilon;
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> v_epsilon(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetaperp_s(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetapar_s(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_theta_s(true); //=v_thetaperp_s+v_thetapar_s

  //interpolated spin vector variation required for inertia forces: v_theta=v_thetaperp+v_thetapar
  std::vector<LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,3> > v_theta(gausspoints.nquad,LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,3>(true));

  //Further material and spatial strains and forces to be evaluated at the Gauss points
  LINALG::TMatrix<FAD,3,1> K; //material curvature
  LINALG::TMatrix<FAD,3,1> Omega; //material deformation measure Omega:=K-K0
  LINALG::TMatrix<FAD,3,1> m; //spatial moment stress resultant
  LINALG::TMatrix<FAD,3,1> M; //material moment stress resultant
  FAD f_par; //material=spatial axial force component

  //Triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_ref_cp(BEAM3K_COLLOCATION_POINTS);  //reference triads (SR-system) at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(BEAM3K_COLLOCATION_POINTS);  //material triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,1> > theta_cp(BEAM3K_COLLOCATION_POINTS);  //relative angle at collocation points

  //Interpolated material triad and angles evaluated at Gauss point
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat(gausspoints.nquad,LINALG::TMatrix<FAD,3,3>(true)); //vector of material triads at gps
  LINALG::TMatrix<FAD,3,1> theta;    //interpolated angle theta
  LINALG::TMatrix<FAD,3,1> theta_s;  //derivative of theta with respect to arc-length s

  //matrices holding the assembled shape functions and s-derivatives
  LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS> N_s;
  LINALG::TMatrix<FAD,1,6*nnode+BEAM3K_COLLOCATION_POINTS> L;
  LINALG::TMatrix<FAD,1,6*nnode+BEAM3K_COLLOCATION_POINTS> L_s;

  //Matrices for individual shape functions and xi-derivatives
  LINALG::TMatrix<FAD,1,2*nnode> N_i_xi;
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i;
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i_xi;

  //Matrices for individual s-derivatives
  LINALG::TMatrix<FAD,1,2*nnode> N_i_s;
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i_s;

  //Additional kinematic quantities
  LINALG::TMatrix<FAD,3,1> r_s; //Matrix to store r'
  FAD abs_r_s; // ||r'||

  //MISC
  double xi=0.0;  //parameter coordinated
  int ind=0; //position index where CP quantities have to be stored

  // unshift node positions, i.e. manipulate element displacement vector
  // as if there where no periodic boundary conditions
  if(StatMechParamsInterfacePtr() != Teuchos::null)
    UnShiftNodePosition(disp,nnode);

  //Set current positions and orientations at all nodes:
  UpdateDispTotlag(disp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  //Store nodal tangents in class variable
  for(int i=0;i<3;i++)
  {
    T_[0](i)=disp_totlag_centerline[3+i].val();
    T_[1](i)=disp_totlag_centerline[10+i].val();
  }

  //********begin: evaluate quantities at collocation points********************************
  for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    //calculate xi of cp
    //node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;

    //get value of interpolating function of theta (lagrange polynomials) at xi
    L_i.Clear();
    N_i_xi.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

    //Determine storage position for the node node
    ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi,jacobi_cp_[ind],N_s);
    r_s.Clear();
    //Calculation of r' at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];
      }
    }
    //calculate epsilon at collocation point
    abs_r_s=FADUTILS::Norm<FAD>(r_s);
    epsilon_cp[ind]=abs_r_s-1.0;

    AssembleShapefunctionsL(L_i,L);

    v_epsilon_cp[ind].Clear();
    v_epsilon_cp[ind].MultiplyTN(N_s,r_s);
    v_epsilon_cp[ind].Scale(1.0/abs_r_s);

    v_thetaperp_cp[ind].Clear();
    LINALG::TMatrix<FAD,3,3> auxmatrix(true);
    LARGEROTATIONS::computespin<FAD>(auxmatrix,r_s);
    v_thetaperp_cp[ind].MultiplyTN(N_s,auxmatrix);
    v_thetaperp_cp[ind].Scale(-1.0/(abs_r_s*abs_r_s));

    v_thetapar_cp[ind].Clear();
    for (int i=0;i<6*nnode+BEAM3K_COLLOCATION_POINTS;i++)
      for (int j=0;j<3;j++)
        (v_thetapar_cp[ind])(i,j)=L(i)*r_s(j)/abs_r_s;
  } //for (int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)

  //calculate angle at cp (this has to be done in a SEPARATE loop as follows)
  for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    theta_cp[node].Clear();
    triadtoangleright(theta_cp[node],triad_mat_cp[REFERENCE_NODE],triad_mat_cp[node]);
  }
  //********end: evaluate quantities at collocation points********************************

  //Clear energy in the beginning
  Eint_=0.0;

  //******begin: gauss integration for internal force vector and stiffness matrix*********
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //Evaluate shape functions
    L_i.Clear();
    L_i_xi.Clear();
    L_i_s.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_1D_deriv1(L_i_xi,xi,Shape());
    L_i_s.Update(1.0/jacobi_[numgp],L_i_xi,0.0);


    //Calculate collocation point interpolations ("v"-vectors and epsilon)
    v_epsilon.Clear();
    v_thetaperp_s.Clear();
    v_thetapar_s.Clear();
    v_theta_s.Clear();
    epsilon=0.0;
    theta.Clear();
    theta_s.Clear();
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      v_epsilon.Update(L_i(node),v_epsilon_cp[node],1.0);
      v_thetaperp_s.Update(L_i_s(node),v_thetaperp_cp[node],1.0);
      v_thetapar_s.Update(L_i_s(node),v_thetapar_cp[node],1.0);
      epsilon+=L_i(node)*epsilon_cp[node];
      theta.Update(L_i(node),theta_cp[node],1.0);
      theta_s.Update(L_i_s(node),theta_cp[node],1.0);
    }
    v_theta_s.Update(1.0,v_thetaperp_s,1.0);
    v_theta_s.Update(1.0,v_thetapar_s,1.0);
    //"v"-vector which is required for inertia forces is already calculated here
    if (massmatrix != NULL or inertia_force != NULL)
    {
      v_theta[numgp].Clear();
      for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
      {
        v_theta[numgp].Update(L_i(node),v_thetaperp_cp[node],1.0);
        v_theta[numgp].Update(L_i(node),v_thetapar_cp[node],1.0);
      }
    }

    //compute material strain K
    K.Clear();
    Omega.Clear();
    computestrain(theta,theta_s,K);

    for(int i=0; i<3; i++)
    {
      Omega(i)=K(i)-K0_[numgp](i);
    }

    //compute material stress resultants
    M.Clear();
    f_par=0.0;
    straintostress(Omega,epsilon,M,f_par);

    //Calculate internal energy and store it in class variable
    Eint_ += 0.5*epsilon.val()*f_par.val()*wgt*jacobi_[numgp];
    for(int i=0; i<3; i++)
    {
      Eint_ += 0.5*Omega(i).val()*M(i).val()*wgt*jacobi_[numgp];
    }

    //compute material triad at gp
    triad_mat[numgp].Clear();
    angletotriad(theta,triad_mat_cp[REFERENCE_NODE],triad_mat[numgp]);

    //pushforward of stress resultants
    m.Clear();
    m.Multiply(triad_mat[numgp],M);

    //residual contribution from moments
    f_int_aux.Clear();
    f_int_aux.Multiply(v_theta_s,m);
    f_int_aux.Scale(wgt*jacobi_[numgp]);
    f_int.Update(1.0,f_int_aux,1.0);

    //residual contribution from axial force
    f_int_aux.Clear();
    f_int_aux.Update(1.0,v_epsilon,0.0);
    f_int_aux.Scale(wgt*jacobi_[numgp]*f_par);
    f_int.Update(1.0,f_int_aux,1.0);
  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  //******end: gauss integration for internal force vector and stiffness matrix*********

  if(rotvec_==true)
  {
    ApplyRotVecTrafo(disp_totlag_centerline,f_int);
  }

  //Update stiffness matrix and force vector
  if(stiffmatrix!=NULL)
  {
    //Calculating stiffness matrix with FAD
    for(int i = 0; i < 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
    {
      for(int j = 0; j < 6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }
    if(rotvec_==true)
    {
      TransformStiffMatrixMultipl(stiffmatrix,disp_totlag);
    }
  }
  if(force!=NULL)
  {
    for (int i=0; i< 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
    {
      (*force)(i)=f_int(i).val();
    }
  }

  //calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
  //finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
  //needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
  //steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
  //consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
  //a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
  //our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
  //according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)

  if ( (massmatrix != NULL or inertia_force != NULL) and !statmechprob_)
  {
    double dt = 1000.0;
    double beta = -1.0;
    double alpha_f = -1.0;
    double alpha_m = -1.0;

    if (this->IsParamsInterface())
    {
      dt = ParamsInterface().GetDeltaTime();
      beta = ParamsInterface().GetBeamParamsInterfacePtr()->GetBeta();
      alpha_f = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlphaf();
      alpha_m = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlpham();
    }
    else
    {
      beta = params.get<double>("rot_beta",1000);
      alpha_f = params.get<double>("rot_alphaf",1000);
      alpha_m = params.get<double>("rot_alpham",1000);
      dt = params.get<double>("delta time",1000);
    }


    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_inert(true);
    CalculateInertiaForces(params,triad_mat,disp_totlag_centerline,v_theta,f_inert);


    if(rotvec_==true)
    {
      ApplyRotVecTrafo(disp_totlag_centerline,f_inert);
    }

    //Update mass matrix and inertia force vector
    if(massmatrix!=NULL)
    {
      //Calculating stiffness matrix with FAD
      for(int i = 0; i < 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
        for(int j = 0; j < 6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
          (*massmatrix)(i,j)=f_inert(i).dx(j);

      if(rotvec_==true)
        TransformStiffMatrixMultipl(massmatrix,disp_totlag);


      // In Lie group GenAlpha algorithm, the mass matrix is multiplied with factor (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)) later.
      // so we apply inverse factor here because the correct prefactors for linearization of displacement/velocity/acceleration dependent terms have been applied automatically by FAD
      massmatrix->Scale(beta*dt*dt*(1.0-alpha_f)/(1.0-alpha_m));
    }

    if(inertia_force!=NULL)
    {
      for (int i=0; i< 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
        (*inertia_force)(i)=f_inert(i).val();
    }

  }

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | internal forces: strong Kirchhoff constraint (private)                                          meier 02/16|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::CalculateInternalForcesSK( Teuchos::ParameterList& params,
                                                      std::vector<double>&      disp,
                                                      Epetra_SerialDenseMatrix* stiffmatrix,
                                                      Epetra_SerialDenseMatrix* massmatrix,
                                                      Epetra_SerialDenseVector* force,
                                                      Epetra_SerialDenseVector* inertia_force)
{
  TEUCHOS_FUNC_TIME_MONITOR("Beam3k::CalculateInternalForcesSK");
  if(BEAM3K_COLLOCATION_POINTS!=2 and BEAM3K_COLLOCATION_POINTS!=3 and BEAM3K_COLLOCATION_POINTS!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for BEAM3K_COLLOCATION_POINTS!!!");
  }

  //number of nodes fixed for these element
  const int nnode = 2;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int_aux(true);

  //vector for current nodal DoFs in total Lagrangian style, i.e. displacement + initial values:
  //rotvec_==true:  disp_totlag=[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3]
  //rotvec_==false: disp_totlag=[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3]
  //! The Number of collocation points can take on the values 2, 3 and 4. 3 and 4 are interior nodes.
  //This leads e.g. in the case rotvec_==true to the following ordering:
  //if BEAM3K_COLLOCATION_POINTS = 2: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2]
  //if BEAM3K_COLLOCATION_POINTS = 3: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3]
  //if BEAM3K_COLLOCATION_POINTS = 4: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3 \alpha_4]
  std::vector<FAD> disp_totlag(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);

  //vector containing locally assembled nodal positions and tangents required for centerline:
  //r(s)=N(s)*disp_totlag_centerline, with disp_totlag_centerline=[\v{d}_1, \v{t}_1, 0, \v{d}_2, \v{t}_2, 0, 0]
  std::vector<FAD> disp_totlag_centerline(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);

  //CP values of strains and their variations needed for interpolation
  std::vector<FAD> epsilon_cp(BEAM3K_COLLOCATION_POINTS); // axial tension epsilon=|r_s|-1 at collocation points
  std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> > v_epsilon_cp(BEAM3K_COLLOCATION_POINTS);

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  //interpolated values of strains and their variations evaluated at Gauss points
  FAD epsilon;
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> v_epsilon(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetaperp_s(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetapartheta_s(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_theta_s(true); //=v_thetaperp_s+v_thetapartheta_s(+v_thetapard_s)

  //interpolated spin vector variation required for inertia forces: v_theta=v_thetaperp+v_thetapartheta(+v_thetapard)
  std::vector<LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,3> > v_theta(gausspoints.nquad,LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,3>(true));

  //Further material and spatial strains and forces to be evaluated at the Gauss points
  LINALG::TMatrix<FAD,3,1> K(true); //material curvature
  LINALG::TMatrix<FAD,3,1> Omega(true); //material deformation measure Omega:=K-K0
  LINALG::TMatrix<FAD,3,1> m(true); //spatial moment stress resultant
  LINALG::TMatrix<FAD,3,1> M(true); //material moment stress resultant
  FAD f_par=0.0; //material=spatial axial force component

  //Triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_ref_cp(BEAM3K_COLLOCATION_POINTS,LINALG::TMatrix<FAD,3,3>(true));  //reference triads (SR-system) at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(BEAM3K_COLLOCATION_POINTS,LINALG::TMatrix<FAD,3,3>(true));  //material triads at collocation points
  std::vector<FAD > phi_cp(BEAM3K_COLLOCATION_POINTS,0.0);  //relative angle at collocation points

  //Interpolated material triad and angles evaluated at Gauss point
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat(gausspoints.nquad,LINALG::TMatrix<FAD,3,3>(true)); //vector of material triads at gps
  FAD phi=0.0;    //interpolated relative angle phi
  FAD phi_s=0.0;  //derivative of interpolated relative angle phi with respect to arc-length s

  //matrices holding the assembled shape functions and s-derivatives
  LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS> N_s(true);
  LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS> N_ss(true);
  LINALG::TMatrix<FAD,1,6*nnode+BEAM3K_COLLOCATION_POINTS> L(true);
  LINALG::TMatrix<FAD,1,6*nnode+BEAM3K_COLLOCATION_POINTS> L_s(true);

  //Matrices for individual shape functions and xi-derivatives
  LINALG::TMatrix<FAD,1,2*nnode> N_i_xi(true);
  LINALG::TMatrix<FAD,1,2*nnode> N_i_xixi(true);
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i(true);
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i_xi(true);

  //Matrices for individual s-derivatives
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i_s(true);

  //Additional kinematic quantities
  LINALG::TMatrix<FAD,3,1> r_s(true); //Matrix to store r'
  LINALG::TMatrix<FAD,3,1> r_ss(true); //Matrix to store r''
  LINALG::TMatrix<FAD,3,1> g1(true); //g1:=r'/||r'||
  LINALG::TMatrix<FAD,3,1> g1_s(true); //g1'
  LINALG::TMatrix<FAD,3,1> ttilde(true); //\tilde{t}:=g1/||r'||=r'/||r'||^2
  LINALG::TMatrix<FAD,3,1> ttilde_s(true);//\tilde{t}'
  LINALG::TMatrix<FAD,3,1> kappacl(true); //centerline (cl) curvature vector
  FAD abs_r_s=0.0; // ||r'||
  FAD rsTrss=0.0; // r'^Tr''
  LINALG::TMatrix<FAD,3,3> auxmatrix1(true); //auxilliary matrix
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> auxmatrix2(true); //auxilliary matrix

  #ifdef CONSISTENTSPINSK
    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetapard_s(true);
    std::vector<LINALG::TMatrix<FAD,3,1> > g1_cp(BEAM3K_COLLOCATION_POINTS, LINALG::TMatrix<FAD,3,1>(true));
    std::vector<LINALG::TMatrix<FAD,3,1> > ttilde_cp(BEAM3K_COLLOCATION_POINTS, LINALG::TMatrix<FAD,3,1>(true));
    std::vector<LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS> > N_s_cp(BEAM3K_COLLOCATION_POINTS, LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS>(true));
    std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> > v1_cp(BEAM3K_COLLOCATION_POINTS, LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1>(true));
    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> v1(true);
    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> v1_s(true);
  #endif

  //MISC
  double xi=0.0;  //parameter coordinated
  int ind=0; //position index where CP quantities have to be stored

  // unshift node positions, i.e. manipulate element displacement vector
  // as if there where no periodic boundary conditions
  if(StatMechParamsInterfacePtr() != Teuchos::null)
    UnShiftNodePosition(disp,nnode);

  //Set current positions and orientations at all nodes:
  UpdateDispTotlag(disp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  //Store nodal tangents in class variable
  for(int i=0;i<3;i++)
  {
    T_[0](i)=disp_totlag_centerline[3+i].val();
    T_[1](i)=disp_totlag_centerline[10+i].val();
  }


  //********begin: evaluate quantities at collocation points********************************
  for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    //calculate xi of cp
    //node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;

    //get value of interpolating function of theta (lagrange polynomials) at xi
    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

    //Determine storage position for the node node
    ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi,jacobi_cp_[ind],N_s);
    r_s.Clear();
    //Calculation of r' at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];
      }
    }
    //calculate epsilon at collocation point
    abs_r_s=FADUTILS::Norm<FAD>(r_s);
    epsilon_cp[ind]=abs_r_s-1.0;

    v_epsilon_cp[ind].Clear();
    v_epsilon_cp[ind].MultiplyTN(N_s,r_s);
    v_epsilon_cp[ind].Scale(1.0/abs_r_s);

    #ifdef CONSISTENTSPINSK
      N_s_cp[ind].Update(1.0,N_s,0.0);
      g1_cp[ind].Update(1.0/abs_r_s,r_s,0.0);
      ttilde_cp[ind].Update(1.0/(abs_r_s*abs_r_s),r_s,0.0);
    #endif

  } //for (int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)

  //calculate angle at cp (this has to be done in a SEPARATE loop as follows)
  for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    LINALG::TMatrix<FAD,3,3> Lambdabarref(true);
    LINALG::TMatrix<FAD,3,1> tangentref(true);
    LINALG::TMatrix<FAD,3,1> phivec(true);
    for(int i=0;i<3;i++)
    {
      tangentref(i)=triad_mat_cp[node](i,0);
    }
    CalculateSRTriads<FAD>(tangentref,triad_mat_cp[REFERENCE_NODE],Lambdabarref);
    triadtoangleleft(phivec,Lambdabarref,triad_mat_cp[node]);
    phi_cp[node]=0.0;
    for(int i=0;i<3;i++)
    {
      phi_cp[node]+=tangentref(i)*phivec(i);
    }

    #ifdef CONSISTENTSPINSK
      LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> auxmatrix3(true);
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s_cp[node], g1_cp[REFERENCE_NODE], ttilde_cp[node], auxmatrix3);
      v1_cp[node].Update(1.0,auxmatrix3,0.0);
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s_cp[REFERENCE_NODE], g1_cp[node], ttilde_cp[REFERENCE_NODE], auxmatrix3);
      v1_cp[node].Update(-1.0,auxmatrix3,1.0);
      LINALG::TMatrix<FAD,1,1> auxscalar(true);
      auxscalar.MultiplyTN(g1_cp[node],g1_cp[REFERENCE_NODE]);
      v1_cp[node].Scale(1.0/(1.0+auxscalar(0,0)));
    #endif
  }
  //********end: evaluate quantities at collocation points********************************

  //Clear energy in the beginning
  Eint_=0.0;

  //******begin: gauss integration for internal force vector and stiffness matrix*********
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //Evaluate and assemble shape functions
    L_i.Clear();
    L_i_xi.Clear();
    L_i_s.Clear();
    L.Clear();
    L_s.Clear();
    N_i_xi.Clear();
    N_i_xixi.Clear();
    N_s.Clear();
    N_ss.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_1D_deriv1(L_i_xi,xi,Shape());
    L_i_s.Update(1.0/jacobi_[numgp],L_i_xi,0.0);
    AssembleShapefunctionsL(L_i,L);
    //The assemble routine is identical for L and L_s
    AssembleShapefunctionsL(L_i_s,L_s);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    AssembleShapefunctionsNs(N_i_xi,jacobi_[numgp],N_s);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi,xi,length_,line2);
    AssembleShapefunctionsNss(N_i_xi,N_i_xixi,jacobi_[numgp],jacobi2_[numgp],N_ss);

    //Calculate collocation piont interpolations
    v_epsilon.Clear();
    epsilon=0.0;
    phi=0.0;
    phi_s=0.0;
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      //calculate interpolated axial tension and variation
      v_epsilon.Update(L_i(node),v_epsilon_cp[node],1.0);
      epsilon+=L_i(node)*epsilon_cp[node];
      //calculate interpolated relative angle
      phi+=L_i(node)*phi_cp[node];
      phi_s+=L_i_s(node)*phi_cp[node];
    }

    //Calculation of r' and r'' at xi
    r_s.Clear();
    r_ss.Clear();
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];
        r_ss(i)+=N_ss(i,j)*disp_totlag_centerline[j];
      }
    }

    //*****************************************************************************************************************************
    //************************Begin: Determine "v"-vectors representing the discrete strain variations*****************************
    //*****************************************************************************************************************************
    //Auxilliary quantities
    abs_r_s=0.0;
    rsTrss=0.0;
    abs_r_s=FADUTILS::Norm<FAD>(r_s);
    for (int i=0;i<3;i++)
    {
      rsTrss+=r_s(i)*r_ss(i);
    }
    g1.Clear();
    g1_s.Clear();
    g1.Update(1.0/abs_r_s,r_s,0.0);
    g1_s.Update(1.0/abs_r_s,r_ss,0.0);
    g1_s.Update(-rsTrss/(abs_r_s*abs_r_s*abs_r_s),r_s,1.0);
    ttilde.Clear();
    ttilde_s.Clear();
    ttilde.Update(1.0/(abs_r_s*abs_r_s),r_s,0.0);
    ttilde_s.Update(1.0/(abs_r_s*abs_r_s),r_ss,0.0);
    ttilde_s.Update(-2*rsTrss/(abs_r_s*abs_r_s*abs_r_s*abs_r_s),r_s,1.0);

    //************** I) Compute v_theta_s=v_thetaperp_s+v_thetapartheta_s(+v_thetapard_s) *****************************************
    // I a) Compute v_thetapartheta_s
    v_thetapartheta_s.Clear();
    for (int row=0;row<6*nnode+BEAM3K_COLLOCATION_POINTS;row++)
    {
      for (int column=0;column<3;column++)
      {
        v_thetapartheta_s(row,column)+=L_s(0,row)*g1(column, 0) + L(0,row)*g1_s(column, 0);
      }
    }

    // I b) Compute v_thetaperp_s
    v_thetaperp_s.Clear();
    auxmatrix1.Clear();
    auxmatrix2.Clear();
    LARGEROTATIONS::computespin(auxmatrix1,ttilde);
    auxmatrix2.MultiplyTN(N_ss,auxmatrix1);
    v_thetaperp_s.Update(-1.0,auxmatrix2,0.0);
    auxmatrix1.Clear();
    auxmatrix2.Clear();
    LARGEROTATIONS::computespin(auxmatrix1,ttilde_s);
    auxmatrix2.MultiplyTN(N_s,auxmatrix1);
    v_thetaperp_s.Update(-1.0,auxmatrix2,1.0);

    // I c) Calculate sum v_theta_s=v_thetaperp_s+v_thetapartheta_s
    v_theta_s.Clear();
    v_theta_s.Update(1.0,v_thetaperp_s,1.0);
    v_theta_s.Update(1.0,v_thetapartheta_s,1.0);

    //************** II) Compute v_theta=v_thetaperp_+v_thetapartheta_(+v_thetapard_)  which is required for inertia forces ********
    if (massmatrix != NULL or inertia_force != NULL)
    {
      // II a) v_thetapartheta contribution
      v_theta[numgp].Clear();
      for (int row=0;row<6*nnode+BEAM3K_COLLOCATION_POINTS;row++)
      {
        for (int column=0;column<3;column++)
        {
          v_theta[numgp](row,column)+=L(0,row)*g1(column, 0);
        }
      }

      // II b) Compute v_thetaperp contribution
      auxmatrix1.Clear();
      auxmatrix2.Clear();
      LARGEROTATIONS::computespin(auxmatrix1,ttilde);
      auxmatrix2.MultiplyTN(N_s,auxmatrix1);
      v_theta[numgp].Update(-1.0,auxmatrix2,1.0);
    }

    // Compute contributions stemming from CONSISTENTSPINSK
    #ifdef CONSISTENTSPINSK
      //************** to I) Compute v_thetapard_s of v_theta_s=v_thetaperp_s+v_thetapartheta_s(+v_thetapard_s) *******************
      LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> auxmatrix3(true);
      LINALG::TMatrix<FAD,1,1> auxscalar1(true);

      //Calculate v1:
      v1.Clear();
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s, g1_cp[REFERENCE_NODE], ttilde, auxmatrix3);
      v1.Update(1.0,auxmatrix3,0.0);
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s_cp[REFERENCE_NODE], g1, ttilde_cp[REFERENCE_NODE], auxmatrix3);
      v1.Update(-1.0,auxmatrix3,1.0);
      auxscalar1.Clear();
      auxscalar1.MultiplyTN(g1,g1_cp[REFERENCE_NODE]);
      v1.Scale(1.0/(1.0+auxscalar1(0,0)));

      //Calculate v1_s:
      v1_s.Clear();
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s, g1_cp[REFERENCE_NODE], ttilde_s, auxmatrix3);
      v1_s.Update(1.0,auxmatrix3,0.0);
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_ss, g1_cp[REFERENCE_NODE], ttilde, auxmatrix3);
      v1_s.Update(1.0,auxmatrix3,1.0);
      auxmatrix3.Clear();
      ComputeTripleProduct<6*nnode+BEAM3K_COLLOCATION_POINTS>(N_s_cp[REFERENCE_NODE], g1_s, ttilde_cp[REFERENCE_NODE], auxmatrix3);
      v1_s.Update(-1.0,auxmatrix3,1.0);
      auxscalar1.Clear();
      auxscalar1.MultiplyTN(g1_s,g1_cp[REFERENCE_NODE]);
      v1_s.Update(-auxscalar1(0,0),v1,1.0);
      auxscalar1.Clear();
      auxscalar1.MultiplyTN(g1,g1_cp[REFERENCE_NODE]);
      v1_s.Scale(1.0/(1.0+auxscalar1(0,0)));

      //Calculate vec1 and vec2
      LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> vec1(true);
      LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> vec2(true);
      for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
      {
        vec1.Update(L_i_s(node),v1_cp[node],1.0);
        vec2.Update(L_i(node),v1_cp[node],1.0);
      }
      vec1.Update(-1.0,v1_s,1.0);
      vec2.Update(-1.0,v1,1.0);

      //Compute v_thetapard_s
      v_thetapard_s.Clear();
      for (int i=0; i<6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
      {
        for (int j=0; j<3; j++)
        {
          v_thetapard_s(i,j)+=vec1(i)*g1(j)+vec2(i)*g1_s(j);
        }
      }
      // I d) Add v_thetapard_s contribution according to v_theta_s+=v_thetapard_s
      v_theta_s.Update(1.0,v_thetapard_s,1.0);

      // to II)  Compute v_thetapard_ of v_theta=v_thetaperp_+v_thetapartheta_(+v_thetapard_)  which is required for inertia forces*******
      // II c) v_thetapard_ contribution
      if (massmatrix != NULL or inertia_force != NULL)
      {

        LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> vec3(true);
        for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
        {
          vec3.Update(L_i(node),v1_cp[node],1.0);
        }
        vec3.Update(-1.0,v1,1.0);
        for (int i=0; i<6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
        {
          for (int j=0; j<3; j++)
          {
            v_theta[numgp](i,j)+=vec3(i)*g1(j);
          }
        }
      }
    #endif
    //***************************************************************************************************************************
    //************************End: Determine "v"-vectors representing the discrete strain variations*****************************
    //***************************************************************************************************************************

    //Compute material triad and centerline curvature at Gauss point
    triad_mat[numgp].Clear();
    ComputeTriadSK(phi,r_s,triad_mat_cp[REFERENCE_NODE],triad_mat[numgp]);
    kappacl.Clear();
    Calculate_clcurvature(r_s, r_ss, kappacl);

    //compute material strain K at Gauss point
    K.Clear();
    Omega.Clear();
    computestrainSK(phi_s,kappacl,triad_mat_cp[REFERENCE_NODE],triad_mat[numgp],K);
    for(int i=0; i<3; i++)
    {
      Omega(i)=K(i)-K0_[numgp](i);
    }

    //compute material stress resultants at Gauss point
    M.Clear();
    f_par=0.0;
    straintostress(Omega,epsilon,M,f_par);

    //Calculate internal energy and store it in class variable
    Eint_ += 0.5*epsilon.val()*f_par.val()*wgt*jacobi_[numgp];
    for(int i=0; i<3; i++)
    {
      Eint_ += 0.5*Omega(i).val()*M(i).val()*wgt*jacobi_[numgp];
    }

    //pushforward of stress resultants
    m.Clear();
    m.Multiply(triad_mat[numgp],M);

    //residual contribution from moments
    f_int_aux.Clear();
    f_int_aux.Multiply(v_theta_s,m);
    f_int_aux.Scale(wgt*jacobi_[numgp]);
    f_int.Update(1.0,f_int_aux,1.0);

    //residual contribution from axial forces
    f_int_aux.Clear();
    f_int_aux.Update(1.0,v_epsilon,0.0);
    f_int_aux.Scale(wgt*jacobi_[numgp]*f_par);
    f_int.Update(1.0,f_int_aux,1.0);
  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  //******end: gauss integration for internal force vector and stiffness matrix*********

  if(rotvec_==true)
  {
    ApplyRotVecTrafo(disp_totlag_centerline,f_int);
  }

  //Update stiffness matrix and force vector
  if(stiffmatrix!=NULL)
  {
    //Calculating stiffness matrix with FAD
    for(int i = 0; i < 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
    {
      for(int j = 0; j < 6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }
    if(rotvec_==true)
    {
        TransformStiffMatrixMultipl(stiffmatrix,disp_totlag);
    }

  }
  if(force!=NULL)
  {
    for (int i=0; i< 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
    {
      (*force)(i)=f_int(i).val();
    }
  }

  //calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
  //finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
  //needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
  //steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
  //consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
  //a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
  //our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
  //according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)

  if ( (massmatrix != NULL or inertia_force != NULL) and !statmechprob_)
  {
    double dt = 1000.0;
    double beta = -1.0;
    double alpha_f = -1.0;
    double alpha_m = -1.0;

    if (this->IsParamsInterface())
    {
      dt = ParamsInterface().GetDeltaTime();
      beta = ParamsInterface().GetBeamParamsInterfacePtr()->GetBeta();
      alpha_f = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlphaf();
      alpha_m = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlpham();
    }
    else
    {
      beta = params.get<double>("rot_beta",1000);
      alpha_f = params.get<double>("rot_alphaf",1000);
      alpha_m = params.get<double>("rot_alpham",1000);
      dt = params.get<double>("delta time",1000);
    }

    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_inert(true);
    CalculateInertiaForces(params,triad_mat,disp_totlag_centerline,v_theta,f_inert);

    if(rotvec_==true)
    {
      ApplyRotVecTrafo(disp_totlag_centerline,f_inert);
    }

    //Update mass matrix and inertia force vector
    if(massmatrix!=NULL)
    {
      //Calculating stiffness matrix with FAD
      for(int i = 0; i < 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
        for(int j = 0; j < 6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
          (*massmatrix)(i,j)=f_inert(i).dx(j);

      if(rotvec_==true)
          TransformStiffMatrixMultipl(massmatrix,disp_totlag);

      // In Lie group GenAlpha algorithm, the mass matrix is multiplied with factor (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)) later.
      // so we apply inverse factor here because the correct prefactors for linearization of displacement/velocity/acceleration dependent terms have been applied automatically by FAD
      massmatrix->Scale(beta*dt*dt*(1.0-alpha_f)/(1.0-alpha_m));
    }
    if(inertia_force!=NULL)
    {
      for (int i=0; i< 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
        (*inertia_force)(i)=f_inert(i).val();
    }

  }

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | Calculate inertia forces and moments                                                           meier 02/16|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::CalculateInertiaForces(Teuchos::ParameterList& params,
                                                    std::vector<LINALG::TMatrix<FAD,3,3> >& triad_mat,
                                                    std::vector<FAD>& disp_totlag_centerline,
                                                    std::vector<LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,3> >& v_theta,
                                                    LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,1>& f_inert)
{
  double dt = 1000.0;
  double beta = -1.0;
  double gamma = -1.0;
  double alpha_f = -1.0;
  double alpha_m = -1.0;

  if (this->IsParamsInterface())
  {
    dt = ParamsInterface().GetDeltaTime();
    beta = ParamsInterface().GetBeamParamsInterfacePtr()->GetBeta();
    gamma = ParamsInterface().GetBeamParamsInterfacePtr()->GetGamma();
    alpha_f = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlphaf();
    alpha_m = ParamsInterface().GetBeamParamsInterfacePtr()->GetAlpham();
  }
  else
  {
    beta = params.get<double>("rot_beta",1000);
    gamma = params.get<double>("rot_gamma",1000);
    alpha_f = params.get<double>("rot_alphaf",1000);
    alpha_m = params.get<double>("rot_alpham",1000);
    dt = params.get<double>("delta time",1000);
  }

  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double rho = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      rho = actmat->Density();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
    break;
  }

  //tensor of mass moments of inertia and cross-section value. These values are used in order to artificially scale
  //the the translational and rotational inertia terms with given input parameters if necessary:
  LINALG::TMatrix<FAD,3,3> Jp(true);
  Jp(0,0)=inertscalerot1_*(Iyy_+Izz_);
  Jp(1,1)=inertscalerot2_*Iyy_;
  Jp(2,2)=inertscalerot2_*Izz_;
  Jp.Scale(rho);

  const int nnode=2;

  double scaledcrosssec=inertscaletrans_*crosssec_;

  LINALG::TMatrix<FAD,3,6*nnode+BEAM3K_COLLOCATION_POINTS> N(true);
  LINALG::TMatrix<FAD,1,2*nnode> N_i(true);
  LINALG::TMatrix<FAD,3,1> rnewmass(true); //Matrix to store r
  LINALG::TMatrix<FAD,3,3> triad_mat_old(true);

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_inert_aux(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetaperp(true);
  LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetapar(true);

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  Ekin_=0.0;

  for (int numgp=0; numgp<gausspoints.nquad; numgp++)//loop through Gauss points
  {
    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //compute quaterion of material triad at gp
    LINALG::TMatrix<FAD,4,1> Qnewmass(true);
    LARGEROTATIONS::triadtoquaternion(triad_mat[numgp],Qnewmass);

    N_i.Clear();
    N.Clear();
    DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
    AssembleShapefunctionsN(N_i,N);
    rnewmass.Clear();
    //Calculation of r at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
      {
        rnewmass(i)+=N(i,j)*disp_totlag_centerline[j];
      }
    }

    triad_mat_old.Clear();
    LINALG::TMatrix<FAD,4,1>  Qconv(true);
    for(int i=0;i<4;i++)
      Qconv(i)=(Qconvmass_[numgp])(i);

    LARGEROTATIONS::quaterniontotriad(Qconv,triad_mat_old);

    LINALG::TMatrix<FAD,3,3>  deltatriad(true);
    deltatriad.MultiplyNT(triad_mat[numgp],triad_mat_old);
    LINALG::TMatrix<FAD,4,1>  deltaQ(true);
    LARGEROTATIONS::triadtoquaternion(deltatriad,deltaQ);
    LINALG::TMatrix<FAD,3,1> deltatheta(true);
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //compute material counterparts of spatial vectors
    LINALG::TMatrix<FAD,3,1> deltaTHETA(true);
    LINALG::TMatrix<FAD,3,1> Wconvmass(true);
    LINALG::TMatrix<FAD,3,1> Aconvmass(true);
    LINALG::TMatrix<FAD,3,1> Amodconvmass(true);

    deltaTHETA.MultiplyTN(triad_mat[numgp],deltatheta);

    LINALG::TMatrix<FAD,3,1> auxvector(true);
    for(int i=0;i<3;i++)
      auxvector(i)=(wconvmass_[numgp])(i);
    Wconvmass.MultiplyTN(triad_mat_old,auxvector);

    for(int i=0;i<3;i++)
      auxvector(i)=(aconvmass_[numgp])(i);
    Aconvmass.MultiplyTN(triad_mat_old,auxvector);

    for(int i=0;i<3;i++)
      auxvector(i)=(amodconvmass_[numgp])(i);
    Amodconvmass.MultiplyTN(triad_mat_old,auxvector);

    LINALG::TMatrix<FAD,3,1> deltar(true);
    for (int i=0;i<3;i++)
    {
      deltar(i)=rnewmass(i)-rconvmass_[numgp](i);
    }

    LINALG::TMatrix<FAD,3,1> Anewmass(true);
    LINALG::TMatrix<FAD,3,1> Wnewmass(true);
    LINALG::TMatrix<FAD,3,1> Amodnewmass(true);
    LINALG::TMatrix<FAD,3,1> rttnewmass(true);
    LINALG::TMatrix<FAD,3,1> rtnewmass(true);
    LINALG::TMatrix<FAD,3,1> rttmodnewmass(true);

    //update angular velocities and accelerations according to Newmark time integration scheme in
    //material description (see Jelenic, 1999, p. 146, equations (2.8) and (2.9)).
    //The corresponding equations are adapted according to
    //the gen-alpha Lie group time integration scheme proposed in [Arnold, Brüls (2007)], [Brüls, Cardona, 2010]
    //and [Brüls, Cardona, Arnold (2012)]. In the predictor step of the time integration the following
    //formulas automatically deliver a constant displacement (deltatheta=0), consistent velocity and consistent acceleration
    //predictor. This fact has to be reflected in a consistent manner by the choice of the predictor in the input file:
    for (int i=0;i<3;i++)
    {
      Anewmass(i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltaTHETA(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*Wconvmass(i)      \
                    -alpha_f/(1.0-alpha_f)*Aconvmass(i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*Amodconvmass(i);

      Wnewmass(i)=gamma/(beta*dt)*deltaTHETA(i)+(1-gamma/beta)*Wconvmass(i)+dt*(1-gamma/(2*beta))*Amodconvmass(i);


      Amodnewmass(i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*Anewmass(i) + alpha_f*Aconvmass(i) - alpha_m*Amodconvmass(i));
    }

    for (int i=0;i<3;i++)
    {
      rttnewmass(i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltar(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*rtconvmass_[numgp](i)      \
                    -alpha_f/(1.0-alpha_f)*rttconvmass_[numgp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*rttmodconvmass_[numgp](i);

      rtnewmass(i)=gamma/(beta*dt)*deltar(i)+(1-gamma/beta)*rtconvmass_[numgp](i)+dt*(1-gamma/(2*beta))*rttmodconvmass_[numgp](i);

      rttmodnewmass(i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*rttnewmass(i) + alpha_f*rttconvmass_[numgp](i) - alpha_m*rttmodconvmass_[numgp](i) );
    }

    //spin matrix of the material angular velocity, i.e. S(W)
    LINALG::TMatrix<FAD,3,3> SWnewmass(true);
    LARGEROTATIONS::computespin(SWnewmass,Wnewmass);
    LINALG::TMatrix<FAD,3,1> Jp_Wnewmass(true);
    LINALG::TMatrix<FAD,3,1> auxvector1(true);
    LINALG::TMatrix<FAD,3,1> Pi_t(true);
    Jp_Wnewmass.Multiply(Jp,Wnewmass);
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
        auxvector1(i)+=SWnewmass(i,j)*Jp_Wnewmass(j)+Jp(i,j)*Anewmass(j);

    Pi_t.Multiply(triad_mat[numgp],auxvector1);
    LINALG::TMatrix<FAD,3,1> L_t(true);
    L_t.Update(rho*scaledcrosssec,rttnewmass,1.0);

    f_inert_aux.Clear();
    f_inert_aux.Multiply(v_theta[numgp],Pi_t);
    f_inert_aux.Scale(wgt*jacobi_[numgp]);
    f_inert.Update(1.0,f_inert_aux,1.0);

    f_inert_aux.Clear();
    f_inert_aux.MultiplyTN(N,L_t);
    f_inert_aux.Scale(wgt*jacobi_[numgp]);
    f_inert.Update(1.0,f_inert_aux,1.0);

    //Calculation of kinetic energy
    LINALG::TMatrix<FAD,1,1> ekinrot(true);
    LINALG::TMatrix<FAD,1,1> ekintrans(true);
    ekinrot.MultiplyTN(Wnewmass,Jp_Wnewmass);
    ekintrans.MultiplyTN(rtnewmass,rtnewmass);
    Ekin_+=0.5*(ekinrot(0,0).val() + rho*scaledcrosssec*ekintrans(0,0).val())*wgt*jacobi_[numgp];

    //**********begin: update class variables needed for storage**************
    LINALG::TMatrix<FAD,3,1> wnewmass(true);
    LINALG::TMatrix<FAD,3,1> anewmass(true);
    LINALG::TMatrix<FAD,3,1> amodnewmass(true);
    wnewmass.Multiply(triad_mat[numgp],Wnewmass);
    anewmass.Multiply(triad_mat[numgp],Anewmass);
    amodnewmass.Multiply(triad_mat[numgp],Amodnewmass);

    for(int i=0;i<4;i++)
    {
      (Qnewmass_[numgp])(i)=(Qnewmass(i)).val();
    }

    for(int i=0;i<3;i++)
    {
      (wnewmass_[numgp])(i)=(wnewmass(i)).val();
      (anewmass_[numgp])(i)=(anewmass(i)).val();
      (amodnewmass_[numgp])(i)=(amodnewmass(i)).val();

      (rnewmass_[numgp])(i)=(rnewmass(i)).val();
      (rtnewmass_[numgp])(i)=(rtnewmass(i)).val();
      (rttnewmass_[numgp])(i)=(rttnewmass(i)).val();
      (rttmodnewmass_[numgp])(i)=(rttmodnewmass(i)).val();
    }
    //**********end: update class variables needed for storage**************
  }//for (int gp=0; gp<gausspoints.nquad; gp++)
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition (public)                                  meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3k::EvaluateNeumann(Teuchos::ParameterList& params,
                                               DRT::Discretization& discretization,
                                               DRT::Condition& condition,
                                               std::vector<int>& lm,
                                               Epetra_SerialDenseVector& elevec1,
                                               Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  //As long as only endpoint forces and moments as well as distributed forces (i.e. no distributed moments)
  //are considered, the method EvaluateNeumann is identical for the WK and the SK case.

  const int twistdofs = BEAM3K_COLLOCATION_POINTS;
  if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for BEAM3K_COLLOCATION_POINTS!!!");
  }
  //dimensions of freedom per node
  const int nnode=2;

  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

  std::vector<FAD> disp_totlag(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);
  std::vector<FAD> disp_totlag_centerline(6*nnode+BEAM3K_COLLOCATION_POINTS, 0.0);
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(BEAM3K_COLLOCATION_POINTS);

  UpdateDispTotlag(mydisp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  // find out whether we will use a time curve
  bool usetime = true;
  double time = -1.0;
  if (this->IsParamsInterface())
    time = this->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  // amplitude of load curve at current time called
  std::vector<double> curvefac(6,1.0);

  for (int i=0; i<6; ++i)
  {
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[i];

    if (curvenum>=0 && usetime)
      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // get values and switches from the condition:
  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition

  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");

  //find out which node is correct
  const std::vector< int > * nodeids = condition.Nodes();


  //if a point neumann condition needs to be linearized
  if(condition.Type() == DRT::Condition::PointNeumannEB)
  {
    //external force vector
    LINALG::TMatrix<FAD,6*nnode+twistdofs,1> f_ext(true);
    //r' at node
    LINALG::TMatrix<FAD,3,1> r_s(true);
    //|r'| at node
    FAD abs_r_s=0.0;
    //S(r') at node
    LINALG::TMatrix<FAD,3,3> Srs(true);

    //auxiliary quantities
    LINALG::TMatrix<FAD,3,1> auxvector(true);
    LINALG::TMatrix<FAD,1,1> auxscalar(true);

    //matrix for moment at node
    LINALG::TMatrix<FAD,3,1> moment(true);


    //find out local node number --> this is done since the first element of a neumann point condition is used for this function
    //in this case we do not know whether it is the left or the right node. In addition to that, xi is assigned
    //in order to determine the index of the base vectors for the smallest rotation system
    int node = -1;

    if((*nodeids)[0] == Nodes()[0]->Id())
    {
      node = 0;
    }
    else if((*nodeids)[0] == Nodes()[1]->Id())
    {
      node = 1;
    }

    if (node == -1)
      dserror("\nNode could not be found on nodemap!\n");

    //IMPORTANT: fext is multiplied by (-1) in BACI, consequently we need no minus sign here
    if (rotvec_==false)
    {
      for(int i = 0; i < 3 ; i++)
      {
        f_ext(node*7+i)+=(*onoff)[i]*(*val)[i]*curvefac[i];
        moment(i)=(*onoff)[i+3]*(*val)[i+3]*curvefac[i+3];
        r_s(i)=disp_totlag_centerline[node*7+3+i];
        abs_r_s+=r_s(i)*r_s(i);
      }
      abs_r_s=sqrt(abs_r_s);

      LARGEROTATIONS::computespin(Srs,r_s);
      auxvector.Multiply(Srs,moment);
      auxvector.Scale(-1.0/(abs_r_s*abs_r_s));
      auxscalar.MultiplyTN(r_s,moment);
      auxscalar.Scale(1.0/abs_r_s);

      for(int j=0;j<3;j++)
      {
        f_ext(node*7+3+j)+=auxvector(j);
      }
      f_ext(node*7+6)+=auxscalar(0,0);
    }
    else
    {
      for(int i = 0; i < 3 ; i++)
      {
        f_ext(node*7+i)+=(*onoff)[i]*(*val)[i]*curvefac[i];
        f_ext(node*7+3+i)+=(*onoff)[i+3]*(*val)[i+3]*curvefac[i+3];
      }
    }
    //elemat1 is zero for rotvec_==true -> It does not have to be transformed in the case BEAM3K_MULTIPLICATIVEUPDATES
    //IMPORTANT: in contrast to f_ext, elemat1 it is directly added to the stiffness matrix, therefore there has to be a sign change!
    if(elemat1!=NULL)
    {
      //Calculating stiffness matrix with FAD
      for(int i = 0; i < 6*nnode+twistdofs; i++)
      {
        for(int j = 0; j < 6*nnode+twistdofs; j++)
        {
          (*elemat1)(i,j)-=f_ext(i).dx(j);
        }
      }
    }
    for (int i=0; i< 6*nnode+twistdofs; i++)
    {
      elevec1(i)=f_ext(i).val();
    }
  }//if a point neumann condition needs to be linearized
  else if(condition.Type() == DRT::Condition::LineNeumann)//if a line neumann condition needs to be linearized
  {
    //external force vector
    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_ext(true);
    LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,1> f_ext_aux(true);
    std::vector< LINALG::Matrix<3,3> > Gref(2);
    for (int node=0;node<2;node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node],Gref[node]);
    }

    // gaussian points
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);
    LINALG::TMatrix<FAD,1,4> N_i;
    LINALG::TMatrix<FAD,3,6*nnode + BEAM3K_COLLOCATION_POINTS> N;

    //funct is related to the 6 "funct" fields after the val field of the Neumann condition
    //in the input file; funct gives the number of the function defined in the section FUNCT
    const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

    //integration loops
    for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
    {
      //integration points in parameter space and weights
      const double xi = gausspoints.qxg[numgp][0];
      const double wgt = gausspoints.qwgt[numgp];

      //Clear matrix for shape functions
      N_i.Clear();
      N.Clear();
      DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
      AssembleShapefunctionsN(N_i,N);

      //position vector at the gauss point at reference configuration needed for function evaluation
      std::vector<double> X_ref(3,0.0);
      //calculate coordinates of corresponding Guass point in reference configuration
      for (int node=0;node<2;node++)
      {
        for (int dof=0;dof<3;dof++)
        {
          X_ref[dof]+=Nodes()[node]->X()[dof]*(N_i(2*node)).val()+(Gref[node])(dof,0)*(N_i(2*node + 1)).val();
        }
      }

      int functnum = -1;

      //Check if also moment line Neumann conditions are implemented accidentally and throw error
      for (int dof=3; dof<6; ++dof)
      {
        if (functions) functnum = (*functions)[dof];
        else functnum = -1;

        if (functnum>0)
        {
          dserror("Line Neumann conditions for distributed moments are not implemented for beam3k so far! Only the function flag 1, 2 and 3 can be set!");
        }
      }

      double functionfac = 1.0;
      LINALG::TMatrix<FAD,3,1> force(true);

      //sum up load components
      for (int dof=0; dof<3; ++dof)
      {
        if (functions) functnum = (*functions)[dof];
        else functnum = -1;

        if (functnum>0)
        {
          // evaluate function at the position of the current node       --> dof here correct?
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof, &X_ref[0], time, NULL);
        }
        else functionfac = 1.0;

        force(dof)= (*onoff)[dof]*(*val)[dof]*curvefac[dof]*functionfac;
      }

      f_ext_aux.Clear();
      f_ext_aux.MultiplyTN(N,force);
      f_ext_aux.Scale(wgt*jacobi_[numgp]);
      f_ext.Update(1.0,f_ext_aux,1.0);
    } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

    if(rotvec_==true)
    {
      ApplyRotVecTrafo(disp_totlag_centerline,f_ext);
    }

    //Update stiffness matrix and force vector
    //IMPORTANT: in contrast to f_ext, elemat1 it is directly added to the stiffness matrix, therefore there has to be a sign change!
    if(elemat1!=NULL)
    {
      //Calculating stiffness matrix with FAD
      for(int i = 0; i < 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
      {
        for(int j = 0; j < 6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
        {
          (*elemat1)(i,j)=-f_ext(i).dx(j);
        }
      }
      if(rotvec_==true)
      {
          TransformStiffMatrixMultipl(elemat1,disp_totlag);
      }
    }
    for (int i=0; i< 6*nnode+BEAM3K_COLLOCATION_POINTS; i++)
    {
      elevec1(i)=f_ext(i).val();
    }
  }//if a line neumann condition needs to be linearized

  return 0;
}  //DRT::ELEMENTS::Beam3k::EvaluateNeumann


template<unsigned int nnode, unsigned int vpernode, unsigned int ndim>
inline void DRT::ELEMENTS::Beam3k::CalcBrownianForcesAndStiff(Teuchos::ParameterList& params,
                                              std::vector<double>&      vel,  //!< element velocity vector
                                              std::vector<double>&      disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  if (weakkirchhoff_==false)
    dserror("calculation of viscous damping moments not implemented for"
        "WK=0 ('strong' Kirchhoff) yet. Use BEAM3WK elements (set WK=1)!");

  // unshift node positions, i.e. manipulate element displacement vector
  // as if there where no periodic boundary conditions
  if(StatMechParamsInterfacePtr() != Teuchos::null)
    UnShiftNodePosition(disp,nnode);

  //internal force vector
  LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int(true);

  // update current total position state of element

  //vector for current nodal DoFs in total Lagrangian style, i.e. displacement + initial values:
  //rotvec_==true:  disp_totlag=[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3]
  //rotvec_==false: disp_totlag=[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3]
  //! The Number of collocation points can take on the values 2, 3 and 4. 3 and 4 are interior nodes.
  //This leads e.g. in the case rotvec_==true to the following ordering:
  //if BEAM3K_COLLOCATION_POINTS = 2: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2]
  //if BEAM3K_COLLOCATION_POINTS = 3: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3]
  //if BEAM3K_COLLOCATION_POINTS = 4: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3 \alpha_4]
  std::vector<FAD> disp_totlag(nnode*vpernode*ndim+BEAM3K_COLLOCATION_POINTS, 0.0);

  // vector containing locally assembled nodal positions and tangents required for centerline:
  // r(s)=N(s)*disp_totlag_centerline, with disp_totlag_centerline=[\v{d}_1, \v{t}_1, 0, \v{d}_2, \v{t}_2, 0, 0]
  std::vector<FAD> disp_totlag_centerline(nnode*vpernode*ndim+BEAM3K_COLLOCATION_POINTS, 0.0);

  // material triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(BEAM3K_COLLOCATION_POINTS);

  // Set current positions and tangents and triads at all nodes:
  UpdateDispTotlag(disp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  LINALG::TMatrix<FAD,nnode*vpernode*ndim,1> disp_totlag_centerlineDOFs_only(true);
  ExtractCenterlineDofValues<nnode,vpernode,FAD>(disp_totlag_centerline,disp_totlag_centerlineDOFs_only);

  // export current velocity state of element to fixed size matrix
  LINALG::Matrix<nnode*vpernode*ndim,1> vel_centerline(true);

  // update current values of centerline (i.e. translational) velocity
  ExtractCenterlineDofValues<nnode,vpernode,double>(vel,vel_centerline);

  // Evaluation of force vectors and stiffness matrices

  // add stiffness and forces due to translational damping effects
  EvaluateTranslationalDamping<nnode,vpernode,ndim>(params,vel_centerline,disp_totlag_centerlineDOFs_only,stiffmatrix,f_int);

  // add stochastic forces and (if required) resulting stiffness
  EvaluateStochasticForces<nnode,vpernode,ndim,3>(params,disp_totlag_centerlineDOFs_only,stiffmatrix,f_int);

  // add stiffness and forces (i.e. moments) due to rotational damping effects
  EvaluateRotationalDamping<nnode,vpernode,ndim>(params,disp_totlag_centerline,triad_mat_cp,stiffmatrix,f_int);


  if (rotvec_==true)
  {
    ApplyRotVecTrafo(disp_totlag_centerline,f_int);
  }

  // Update stiffness matrix and force vector
  if (stiffmatrix!=NULL)
  {
    // Calculating stiffness matrix with FAD
    for (unsigned int i=0; i<ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS; i++)
      for (unsigned int j=0; j<ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS; j++)
        (*stiffmatrix)(i,j)=f_int(i).dx(j);

    if(rotvec_==true)
      TransformStiffMatrixMultipl(stiffmatrix,disp_totlag);

  }

  if (force!=NULL)
  {
    for (unsigned int i=0; i<ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS; i++)
      (*force)(i)=f_int(i).val();

  }

}


template<unsigned int nnode, unsigned int vpernode, int ndim>
void DRT::ELEMENTS::Beam3k::EvaluateTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                          const LINALG::TMatrix<double,ndim*vpernode*nnode,1>& vel,
                                                          const LINALG::TMatrix<FAD,ndim*vpernode*nnode,1>& disp_totlag,
                                                          Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                          LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,1>& f_int)//!< element internal force vector
{
  // get time step size
  const double dt = ParamsInterface().GetDeltaTime();

  const unsigned int dofpernode = ndim*vpernode+1;

  // get damping coefficients for translational and rotational degrees of freedom (the latter is unused in this element)
  LINALG::Matrix<ndim,1> gamma(true);
  GetDampingCoefficients(gamma);

  // velocity and gradient of background velocity field
  LINALG::TMatrix<FAD,ndim,1> velbackground(true);
  LINALG::TMatrix<FAD,ndim,ndim> velbackgroundgrad(true);

  // evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::TMatrix<FAD,ndim,1> evaluationpoint(true);
  // tangent vector (derivative of beam centerline curve r with respect to arc-length parameter s)
  LINALG::TMatrix<FAD,ndim,1> r_s(true);
  // velocity of beam centerline point relative to background fluid velocity
  LINALG::TMatrix<FAD,ndim,1> vel_rel(true);

  // viscous force vector per unit length at current GP
  LINALG::TMatrix<FAD,ndim,1> f_visc(true);
  // damping matrix
  LINALG::TMatrix<FAD,ndim,ndim> damp_mat(true);


  // get Gauss points and weights
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  // matrix to store individual Hermite shape functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode*vpernode> N_i(true);
  LINALG::Matrix<1,nnode*vpernode> N_i_xi(true);


  for(int gp=0; gp<gausspoints.nquad; gp++)
  {
    EvaluateShapeFunctionsAndDerivsAtXi<nnode,vpernode>(gausspoints.qxg[gp][0],N_i,N_i_xi,this->Shape());

    // compute position vector r of point in physical space corresponding to Gauss point
    Calc_r<nnode,vpernode,FAD>(disp_totlag, N_i, evaluationpoint);

    // compute tangent vector t_{\par}=r' at current Gauss point
    Calc_r_s<nnode,vpernode,FAD>(disp_totlag, N_i_xi, jacobi_[gp], r_s);

    // compute velocity and gradient of background flow field at point r
    GetBackgroundVelocity<ndim,FAD>(params,evaluationpoint,velbackground,velbackgroundgrad);

    // compute velocity vector at this Gauss point via same interpolation as for centerline position vector
    /* Todo: as long as we use FAD, we need to track the dependency of velocity vector on primary variables, i.e. we
             must calculate it like this using the class variable rconvmass_ for position vector at GP of last converged state
             this again is tricky in case of periodoc boundary conditions, because the position of this element (nodes)
             might have been shifted outside; we therefore manually adapt the calculated velocity and compare it with the
             velocity calculated in time integrator and handed in from outside for safety reasons */
    LINALG::TMatrix<double,ndim,1> vel_rel_test;
    CalcInterpolation<nnode,vpernode,3,double>(vel, N_i, vel_rel_test);

    // ************ the following is obsolete as soon as we get rid of FAD linearization ************
    LINALG::Matrix<3,1> diff(true);

    LINALG::TMatrix<FAD,3,1> delta_r_ost(true);
    std::vector<double> periodlength = *(StatMechParamsInterface().GetPeriodLength());
    for (unsigned int idim=0; idim<ndim; idim++)
    {
      // difference in position of this GP as compared to last time step
      delta_r_ost(idim) = evaluationpoint(idim) - rconvmass_[gp](idim);

      // manually adapt this difference in case of shifting due to PBC
      if (periodlength[idim] > 0.0 and delta_r_ost(idim) >= 0.5*periodlength[idim])
        delta_r_ost(idim) -= periodlength[idim];
      else if (periodlength[idim] > 0.0 and delta_r_ost(idim) <= -0.5*periodlength[idim])
        delta_r_ost(idim) += periodlength[idim];

      // velocity according to Backward Euler scheme
      vel_rel(idim) = delta_r_ost(idim)/dt;

      diff(idim)=vel_rel(idim).val()-vel_rel_test(idim);

      // set class variable, such that rconvmass_ is available in next time step
      rnewmass_[gp](idim) = evaluationpoint(idim).val();
    }

    // safety check
    if (diff.NormInf() > 1e-10)
    {
      std::cout << "\nrnewmass = "; evaluationpoint.Print(std::cout);
      std::cout << "\nrconvmass_ = "; rconvmass_[gp].Print(std::cout);
      std::cout << "\nvel_rel = "; vel_rel.Print(std::cout);
      std::cout << "\nvel_rel_test = "; vel_rel_test.Print(std::cout);
      dserror("velocity vector computed locally in beam3k differs from OneStepTheta velocity vector (see values above)!");
    }
    // ************ end *****************************************************************************

    vel_rel -= velbackground;

    // loop over lines and columns of damping matrix
    for (unsigned int idim=0; idim<ndim; idim++)
      for (unsigned int jdim=0; jdim<ndim; jdim++)
        damp_mat(idim,jdim) = (idim==jdim)*gamma(1) + (gamma(0) - gamma(1))*r_s(idim)*r_s(jdim);

    // compute viscous force vector per unit length at current GP
    f_visc.Multiply(damp_mat,vel_rel);


    // loop over all nodes used for centerline interpolation
    for (unsigned int inode=0; inode<nnode; inode++)
      // loop over dimensions
      for (unsigned int idim=0; idim<ndim; idim++)
      {
        f_int(inode*dofpernode+idim) += N_i(vpernode*inode)*jacobi_[gp]*gausspoints.qwgt[gp]*f_visc(idim);
        f_int(inode*dofpernode+3+idim) += N_i(vpernode*inode+1)*jacobi_[gp]*gausspoints.qwgt[gp]*f_visc(idim);
      }

  }

}

template<unsigned int nnode, unsigned int vpernode, unsigned int ndim, unsigned int randompergauss>
void DRT::ELEMENTS::Beam3k::EvaluateStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const LINALG::TMatrix<FAD,ndim*vpernode*nnode,1>& disp_totlag,
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,1>& f_int)//!< element internal force vector
{
  const unsigned int dofpernode = ndim*vpernode+1;

  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  GetDampingCoefficients(gamma);

  /* get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5 */
  Teuchos::RCP<Epetra_MultiVector> randomforces = StatMechParamsInterface().GetRandomForces();

  // tangent vector (derivative of beam centerline curve r with respect to arc-length parameter s)
  LINALG::TMatrix<FAD,ndim,1> r_s(true);

  // my random number vector at current GP
  LINALG::Matrix<ndim,1> randnumvec(true);

  // stochastic force vector per unit length at current GP
  LINALG::TMatrix<FAD,ndim,1> f_stoch(true);


  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  //matrix to store hermite shape functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode*vpernode> N_i;
  LINALG::Matrix<1,nnode*vpernode> N_i_xi;

  for(int gp=0; gp<gausspoints.nquad; gp++)
  {
    EvaluateShapeFunctionsAndDerivsAtXi<nnode,vpernode>(gausspoints.qxg[gp][0],N_i,N_i_xi,this->Shape());

    // compute tangent vector t_{\par}=r' at current Gauss point
    Calc_r_s<nnode,vpernode,FAD>(disp_totlag, N_i_xi, jacobi_[gp], r_s);

    // extract random numbers from global vector
    for (unsigned int idim=0; idim<ndim; idim++)
      randnumvec(idim) = (*randomforces)[gp*randompergauss+idim][LID()];

    // compute stochastic force vector per unit length at current GP
    f_stoch.Clear();
    for (unsigned int idim=0; idim<ndim; idim++)
      for (unsigned int jdim=0; jdim<ndim; jdim++)
        f_stoch(idim) += (std::sqrt(gamma(1))*(idim==jdim) + (std::sqrt(gamma(0)) - std::sqrt(gamma(1)))*r_s(idim)*r_s(jdim))*randnumvec(jdim);

    // loop over all nodes used for centerline interpolation
    for (unsigned int inode=0; inode<nnode; inode++)
      // loop over dimensions
      for (unsigned int idim=0; idim<ndim; idim++)
      {
        f_int(inode*dofpernode+idim) -= N_i(vpernode*inode)*f_stoch(idim)*std::sqrt(jacobi_[gp]*gausspoints.qwgt[gp]);
        f_int(inode*dofpernode+3+idim) -= N_i(vpernode*inode+1)*f_stoch(idim)*std::sqrt(jacobi_[gp]*gausspoints.qwgt[gp]);
      }

  }

}

template<unsigned int nnode, unsigned int vpernode, unsigned int ndim>
void DRT::ELEMENTS::Beam3k::EvaluateRotationalDamping(Teuchos::ParameterList&          params,  //!<parameter list
                                              const std::vector<FAD>& disp_totlag_centerline,
                                              const std::vector<LINALG::TMatrix<FAD,3,3> >& triad_mat_cp,
                                              Epetra_SerialDenseMatrix*                stiffmatrix,  //!< element stiffness matrix
                                              LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,1>& f_int)  //!< element internal force vector
{
  // get time step size
  const double dt = ParamsInterface().GetDeltaTime();

  // get damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  GetDampingCoefficients(gamma);

  // get Gauss points and weights for evaluation of viscous damping contributions
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3k);

  LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,1> f_int_aux(true);

  // CP values of strains and their variations needed for interpolation
  std::vector<LINALG::TMatrix<FAD,6*nnode+BEAM3K_COLLOCATION_POINTS,3> > v_thetapar_cp(BEAM3K_COLLOCATION_POINTS);

  // interpolated values of strains and their variations evaluated at Gauss points
  LINALG::TMatrix<FAD,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS,3> v_thetapar(true);

  std::vector<LINALG::TMatrix<FAD,3,1> > theta_cp(BEAM3K_COLLOCATION_POINTS);  //relative angle at collocation points

  // Interpolated material triad and angle evaluated at Gauss point
  LINALG::TMatrix<FAD,3,3> triad_mat(true);
  LINALG::TMatrix<FAD,3,1> theta(true);


  // matrices holding the assembled shape functions and s-derivatives
  LINALG::TMatrix<FAD,3,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS> N_s;
  LINALG::TMatrix<FAD,1,ndim*vpernode*nnode+BEAM3K_COLLOCATION_POINTS> L;

  // Matrices for individual shape functions and xi-derivatives
  LINALG::TMatrix<FAD,1,vpernode*nnode> N_i_xi;
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i;

  // Additional kinematic quantities
  LINALG::TMatrix<FAD,3,1> r_s; //Matrix to store r'
  FAD abs_r_s; // ||r'||


  //********begin: evaluate quantities at collocation points********************************
  for (unsigned int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    // calculate xi of cp
    // node=0->xi=-1  node=1->xi=0  node=2->xi=1
    const double xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;

    // get value of interpolating function of theta (Lagrange polynomials) at xi
    L_i.Clear();
    N_i_xi.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

    // Determine storage position for the node node
    const int ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi,jacobi_cp_[ind],N_s);
    r_s.Clear();
    // Calculation of r' at xi
    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<6*nnode+BEAM3K_COLLOCATION_POINTS; j++)
        r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];

    // calculate epsilon at collocation point
    abs_r_s=FADUTILS::Norm<FAD>(r_s);

    AssembleShapefunctionsL(L_i,L);

    v_thetapar_cp[ind].Clear();
    for (unsigned int i=0;i<6*nnode+BEAM3K_COLLOCATION_POINTS;i++)
      for (unsigned int j=0;j<3;j++)
        (v_thetapar_cp[ind])(i,j)=L(i)*r_s(j)/abs_r_s;
  }

  // calculate angle at cp (this has to be done in a SEPARATE loop as follows)
  for (unsigned int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
  {
    theta_cp[node].Clear();
    triadtoangleright(theta_cp[node],triad_mat_cp[REFERENCE_NODE],triad_mat_cp[node]);
  }


  for (int gp=0; gp<gausspoints.nquad; gp++)//loop through Gauss points
  {
    // Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[gp][0];
    const double wgt = gausspoints.qwgt[gp];

    // Evaluate shape functions
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,this->Shape());

    theta.Clear();
    v_thetapar.Clear();
    for (unsigned int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      theta.Update(L_i(node),theta_cp[node],1.0);
      v_thetapar.Update(L_i(node),v_thetapar_cp[node],1.0);
    }

    // compute material triad at gp
    triad_mat.Clear();
    angletotriad(theta,triad_mat_cp[REFERENCE_NODE],triad_mat);


    // compute quaterion of material triad at gp
    LINALG::TMatrix<FAD,4,1> Qnewmass(true);
    LARGEROTATIONS::triadtoquaternion(triad_mat,Qnewmass);

    for (unsigned int i=0; i<4; i++)
      (Qnewmass_[gp])(i)=(Qnewmass(i)).val();


    LINALG::TMatrix<FAD,4,1> Qconv(true);
    for (unsigned int i=0; i<4; i++)
      Qconv(i)=(Qconvmass_[gp])(i);


    LINALG::TMatrix<FAD,4,1> deltaQ(true);
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconv),Qnewmass,deltaQ);

    LINALG::TMatrix<FAD,3,1> deltatheta(true);
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);


    // angular velocity at this Gauss point according to backward Euler scheme
    LINALG::TMatrix<FAD,3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);

    // compute matrix Lambda*[gamma(2) 0 0 \\ 0 0 0 \\ 0 0 0]*Lambda^t = gamma(2) * g_1 \otimes g_1
    // where g_1 is first base vector, i.e. first column of Lambda
    LINALG::TMatrix<FAD,3,3> g1g1gamma;
    for (unsigned int k=0; k<3; k++)
      for (unsigned int j=0; j<3; j++)
        g1g1gamma(k,j) = triad_mat(k,0)*triad_mat(j,0)*gamma(2);

    // compute vector gamma(2) * g_1 \otimes g_1 * \omega
    LINALG::TMatrix<FAD,3,1> g1g1gammaomega;
    g1g1gammaomega.Multiply(g1g1gamma,omega);

    // residual contribution from viscous damping moment
    f_int_aux.Clear();
    f_int_aux.Multiply(v_thetapar,g1g1gammaomega);
    f_int_aux.Scale(wgt*jacobi_[gp]);
    f_int.Update(1.0,f_int_aux,1.0);
  }

}

/*----------------------------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3k::HowManyRandomNumbersINeed() const
{
  // get Gauss rule for evaluation of stochastic force contributions
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::mygaussrulebeam3k;
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  /* at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   * random numbers for the translational degrees of freedom */
  return (3*gausspoints.nquad);
}


/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble C shape function                                                                     meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::AssembleShapefunctionsL(LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS>& L_i,
                                                    LINALG::TMatrix<FAD,1,2*6+BEAM3K_COLLOCATION_POINTS>& L)
{

  #if defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==2)

  int assembly_L[2*6+BEAM3K_COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2};

  #elif defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==3)

  int assembly_L[2*6+BEAM3K_COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3};

  #elif defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==4)

  int assembly_L[2*6+BEAM3K_COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3,4};

  #else
  dserror("BEAM3K_COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for BEAM3K_COLLOCATION_POINTS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+BEAM3K_COLLOCATION_POINTS; i++)
  {
    if(assembly_L[i]==0)
    {
      L(i)=0.0;
    }
    else
    {
      L(i)=L_i(assembly_L[i]-1);
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble the N_s shape functions                                                              meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::AssembleShapefunctionsNss( LINALG::TMatrix<FAD,1,4>& N_i_xi,
                                                        LINALG::TMatrix<FAD,1,4>& N_i_xixi,
                                                        FAD jacobi,
                                                        FAD jacobi2,
                                                        LINALG::TMatrix<FAD,3,2*6+BEAM3K_COLLOCATION_POINTS>& N_ss)
{
  LINALG::TMatrix<FAD,1,4> N_i_ss(true);
  N_i_ss.Update(1.0/(jacobi*jacobi),N_i_xixi,1.0);
  N_i_ss.Update(-jacobi2/(jacobi*jacobi*jacobi*jacobi),N_i_xi,1.0);

  AssembleShapefunctionsN(N_i_ss,N_ss);

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble the N_s shape functions                                                              meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::AssembleShapefunctionsNs( LINALG::TMatrix<FAD,1,4>& N_i_xi,
                                                       FAD jacobi,
                                                       LINALG::TMatrix<FAD,3,2*6+BEAM3K_COLLOCATION_POINTS>& N_s)
{

  LINALG::TMatrix<FAD,1,4> N_i_s(true);

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1.0/jacobi);

  AssembleShapefunctionsN(N_i_s,N_s);

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble the N shape functions                                                                meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::AssembleShapefunctionsN( LINALG::TMatrix<FAD,1,4>& N_i,
                                                     LINALG::TMatrix<FAD,3,2*6+BEAM3K_COLLOCATION_POINTS>& N)
{

#if defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==2)

int assembly_N[3][2*6+BEAM3K_COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

#elif defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==3)

int assembly_N[3][2*6+BEAM3K_COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

#elif defined(BEAM3K_COLLOCATION_POINTS) && (BEAM3K_COLLOCATION_POINTS==4)

int assembly_N[3][2*6+BEAM3K_COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                       {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                       {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

#else
dserror("BEAM3K_COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for BEAM3K_COLLOCATION_POINTS!!!");
#endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+BEAM3K_COLLOCATION_POINTS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N(j,i)=0;
      }
      else
      {
        N(j,i)=N_i(assembly_N[j][i]-1);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate position vectors from displacement +  initial position                              meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::UpdateDispTotlag(const std::vector<double>& disp,
                                            std::vector<FAD>& disp_totlag)
{
  for (int dof=0;dof<2*6+BEAM3K_COLLOCATION_POINTS;dof++)
    disp_totlag[dof]=disp[dof];

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
          disp_totlag[7*node+ndof]+=Nodes()[node]->X()[ndof];
        }
        else if (ndof<6)
        {
          disp_totlag[7*node+ndof]+=(Tref()[node])(ndof-3);
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
          disp_totlag[7*node+ndof]+=Nodes()[node]->X()[ndof];
        }
        else if (ndof<6)
        {
          //Nothing to do here, rotations are treated below
        }
        else
        {
          //here we have to add the initial length of the tangents at the boundary nodes, i.e. ||r'_i(t=0)||=1:
          disp_totlag[7*node+ndof]+=1.0;
        }
      }//for (int ndof=0;ndof<7;ndof++)//loop over dofs per node

        LINALG::Matrix<3,1> thetanew(true);
        LINALG::Matrix<4,1> deltaQ(true);
        for(int i=0; i<3; i++)
        {
          dispthetanew_[node](i) = disp[7*node+3+i];
        }

        LARGEROTATIONS::angletoquaternion(dispthetanew_[node],deltaQ);

        LINALG::Matrix<4,1> Q0;
        LARGEROTATIONS::angletoquaternion(theta0_[node],Q0);
        LARGEROTATIONS::quaternionproduct(Q0,deltaQ,Qnew_[node]);

        //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
        Qnew_[node].Scale(1.0/Qnew_[node].Norm2());

        //Calculate the new nodal angle thetanew \in ]-PI,PI] -> Here, thetanew \in ]-PI,PI] by quaterniontoangle()
        LARGEROTATIONS::quaterniontoangle(Qnew_[node],thetanew);

        //Finally set rotation values in disp_totlag
        for(int i=0; i<3; i++)
        {
          disp_totlag[7*node+3+i]=thetanew(i);
        }
    }//for (int node=0;node<2;node++)//loop over boundary nodes
  }

  //Next, we have to set variables for FAD
  for (int dof=0;dof<2*6+BEAM3K_COLLOCATION_POINTS;dof++)
  {
    disp_totlag[dof].diff(dof,2*6+BEAM3K_COLLOCATION_POINTS);
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Set positions vectors and tangents at boundary nodes and triads at all CPs                    meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::SetNodalVariables( std::vector<FAD>& disp_totlag,
                                                std::vector<FAD>& disp_totlag_centerline,
                                                std::vector<LINALG::TMatrix<FAD,3,3> >& triad_mat_cp)
{
  //Set positions at boundary nodes
  for (int i=0;i<3;i++)
  {
    disp_totlag_centerline[i]=disp_totlag[i];
    disp_totlag_centerline[7+i]=disp_totlag[7+i];
  }

  LINALG::TMatrix<FAD,3,1> tangent(true);
  LINALG::TMatrix<FAD,3,3> triad_ref(true);
  LINALG::Matrix<3,3> triad_aux(true);
  LINALG::TMatrix<FAD,3,3> triad_aux2(true);
  FAD alpha = 0.0;
  //next, set triads and tangents at boundary nodes
  if(rotvec_==false)
  {
    for (int i=0;i<3;i++)
    {
      disp_totlag_centerline[3+i]=disp_totlag[3+i];
      disp_totlag_centerline[10+i]=disp_totlag[10+i];
    }

    for (int node=0;node<2;node++)
    {
      for (int i=0;i<3;i++)
        tangent(i)=disp_totlag[7*node+3+i];

      alpha=disp_totlag[7*node+6];

      //calculate new sr triads
      triad_ref.Clear();
      triad_aux.Clear();
      triad_aux2.Clear();
      LARGEROTATIONS::quaterniontotriad(Qrefconv_[node],triad_aux);
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
          triad_aux2(i,j)=triad_aux(i,j);

      CalculateSRTriads<FAD>(tangent,triad_aux2,triad_ref);

      //Store nodal reference triad
      LINALG::TMatrix<FAD,4,1> Qref(true);
      LARGEROTATIONS::triadtoquaternion(triad_ref,Qref);
      for(int i=0;i<4;i++)
        Qrefnew_[node](i)=Qref(i).val();

      //calculate material triad
      triad_mat_cp[node].Clear();
      RotateTriad(triad_ref,alpha,triad_mat_cp[node]);
    }
  }
  else
  {
    LINALG::TMatrix<FAD,3,1> theta(true);
    LINALG::TMatrix<FAD,3,3> unity(true);
    for (int node=0;node<2;node++)
    {
      for (int i=0;i<3;i++)
      {
        theta(i)=disp_totlag[7*node+3+i];
        unity(i,i)=1.0;
      }
      angletotriad(theta,unity,triad_mat_cp[node]);
      for (int i=0;i<3;i++)
      {
        disp_totlag_centerline[7*node+3+i]=(triad_mat_cp[node])(i,0)*disp_totlag[7*node+6];
      }
    }
  }

  LINALG::TMatrix<FAD,1,4> N_i_xi(true);
  LINALG::TMatrix<FAD,3,6*2+BEAM3K_COLLOCATION_POINTS> N_s(true);
  LINALG::TMatrix<FAD,1,BEAM3K_COLLOCATION_POINTS> L_i;
  double xi=0.0;
  int ind=0;
  //********begin: evaluate quantities at collocation points********************************
  for(int node=1; node<BEAM3K_COLLOCATION_POINTS-1; node++)
  {
    //calculate xi of cp
    //node=0->xi=-1  node=1->xi=0 node=2->xi=1
    xi=(double)node/(BEAM3K_COLLOCATION_POINTS-1)*2-1.0;
    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());

    //Determine storage position for the node node
    ind=LARGEROTATIONS::NumberingTrafo(node+1, BEAM3K_COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctionsNs(N_i_xi,jacobi_cp_[ind],N_s);
    tangent.Clear();
    //Calculation of r' at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*2+BEAM3K_COLLOCATION_POINTS; j++)
      {
        tangent(i)+=N_s(i,j)*disp_totlag_centerline[j];
      }
    }

    alpha=disp_totlag[7*2+ind-2];

    //calculate new sr triads
    triad_ref.Clear();
    triad_aux.Clear();
    triad_aux2.Clear();
    LARGEROTATIONS::quaterniontotriad(Qrefconv_[ind],triad_aux);
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
        triad_aux2(i,j)=triad_aux(i,j);

    CalculateSRTriads<FAD>(tangent,triad_aux2,triad_ref);

    //Store nodal reference triad
    LINALG::TMatrix<FAD,4,1> Qref(true);
    LARGEROTATIONS::triadtoquaternion(triad_ref,Qref);
    for(int i=0;i<4;i++)
      Qrefnew_[ind](i)=Qref(i).val();

    //calculate material triad
    triad_mat_cp[ind].Clear();
    RotateTriad(triad_ref,alpha,triad_mat_cp[ind]);

  }//for (int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)

  //Store the triads as quaternions also in the case rotvec_==false
  if(rotvec_==false)
  {
    for(int node=0; node<BEAM3K_COLLOCATION_POINTS; node++)
    {
      LINALG::TMatrix<FAD,4,1> Q(true);
      LARGEROTATIONS::triadtoquaternion(triad_mat_cp[node],Q);
      for(int i=0;i<4;i++)
       Qnew_[node](i)=Q(i).val();
    }
  }

  return;
}


/*-----------------------------------------------------------------------------------------------------------*
 |  Pre-multiply trafo matrix if rotvec_==true: \tilde{\vec{f}_int}=\mat{T}^T*\vec{f}_int         meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::ApplyRotVecTrafo( std::vector<FAD>& disp_totlag_centerline,
                                               LINALG::TMatrix<FAD,6*2+BEAM3K_COLLOCATION_POINTS,1>& f_int)
{
  //Trafo matrices:
  LINALG::TMatrix<FAD,4,4> T(true);
  LINALG::TMatrix<FAD,3,1> g_1(true);
  LINALG::TMatrix<FAD,3,3> auxmatrix(true);
  FAD t=0.0;
  LINALG::TMatrix<FAD,4,1> f_aux1(true);
  LINALG::TMatrix<FAD,4,1> f_aux2(true);
  for (int node=0;node<2;node++)
  {
    g_1.Clear();
    t=0.0;
    auxmatrix.Clear();

    for(int i=0;i<3;i++)
    {
      g_1(i)=disp_totlag_centerline[7*node+3+i];
    }

    t=FADUTILS::Norm<FAD>(g_1);
    g_1.Scale(1.0/t);
    LARGEROTATIONS::computespin(auxmatrix,g_1);
    auxmatrix.Scale(-1.0*t);
    T.Clear();

    for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
      {
        T(i,j)=auxmatrix(i,j);
      }
      T(i,3)=g_1(i);
      T(3,i)=g_1(i);
    }

    f_aux1.Clear();
    f_aux2.Clear();
    for(int i=0;i<4;i++)
    {
      f_aux1(i)=f_int(7*node+3+i);
    }

    f_aux2.MultiplyTN(T,f_aux1);

    for(int i=0;i<4;i++)
    {
      f_int(7*node+3+i)=f_aux2(i);
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Transform stiffness matrix in order to solve for multiplicative rotation vector increments    meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3k::TransformStiffMatrixMultipl(Epetra_SerialDenseMatrix* stiffmatrix,
                                                        std::vector<FAD>& disp_totlag)
{
  // we need to transform the stiffmatrix because its entries are derivatives with respect to additive rotational increments
  // we want a stiffmatrix containing derivatives with respect to multiplicative rotational increments
  // therefore apply a trafo matrix to all those 3x3 blocks in stiffmatrix which correspond to derivation with respect to rotational DOFs
  // the trafo matrix is simply the T-Matrix (see Jelenic1999, (2.4)): \Delta_{mult} \vec \theta_{inode} = \mat T(\vec \theta_{inode} * \Delta_{addit} \vec \theta_{inode}
  LINALG::Matrix<2*6+BEAM3K_COLLOCATION_POINTS,3> tempmat(true);
  LINALG::Matrix<2*6+BEAM3K_COLLOCATION_POINTS,3> newstiffmat(true);
  LINALG::Matrix<3,3> Tmat(true);
  std::vector<LINALG::Matrix<3,1> > theta(2,LINALG::Matrix<3,1>(true));

  //Loop over the two boundary nodes
  for(int node=0;node<2;node++)
  {
    for(int i=0;i<3;i++)
      theta[node](i)=disp_totlag[7*node+3+i].val();
  }

  //Loop over the two boundary nodes
  for(int node=0;node<2;node++)
  {
    Tmat.Clear();
    Tmat = LARGEROTATIONS::Tmatrix(theta[node]);

    tempmat.Clear();
    for (int i=0; i<2*6+BEAM3K_COLLOCATION_POINTS; ++i)
      for (int j=0; j<3; ++j)
        tempmat(i,j) = (*stiffmatrix)(i,7*node+3+j);

    newstiffmat.Clear();
    newstiffmat.MultiplyNN(tempmat,Tmat);

    for (int i=0; i<2*6+BEAM3K_COLLOCATION_POINTS; ++i)
      for (int j=0; j<3; ++j)
        (*stiffmatrix)(i,7*node+3+j)=newstiffmat(i,j);
  }

  return;
}

void DRT::ELEMENTS::Beam3k::straintostress(LINALG::TMatrix<FAD,3,1>& Omega, FAD epsilon, LINALG::TMatrix<FAD,3,1>& M, FAD& f_par)
{
  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = actmat->ShearMod();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
    break;
  }

  M.Clear();
  f_par=0.0;

  f_par=ym*crosssec_*epsilon;
  M(0)=sm*Irr_*Omega(0);
  M(1)=ym*Iyy_*Omega(1);
  M(2)=ym*Izz_*Omega(2);

}

//      //*******************************Begin: FD-CHECK************************************************************
//      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
//      //correctly or not by means of a numerically approximated stiffness matrix. Uncomment this code block and copy
//      //it to the marked place in the method DRT::ELEMENTS::Beam3k::Evaluate() on the top of this file!
//      if(Id() == 0) //limiting the following tests to certain element numbers
//      {
//
//        //variable to store numerically approximated stiffness matrix
//        Epetra_SerialDenseMatrix stiff_approx;
//        stiff_approx.Shape(6*2+BEAM3K_COLLOCATION_POINTS,6*2+BEAM3K_COLLOCATION_POINTS);
//
//
//        //relative error of numerically approximated stiffness matrix
//        Epetra_SerialDenseMatrix stiff_relerr;
//        stiff_relerr.Shape(6*2+BEAM3K_COLLOCATION_POINTS,6*2+BEAM3K_COLLOCATION_POINTS);
//
//        //characteristic length for numerical approximation of stiffness
//        double h_rel = 1e-7;
//
//        //flag indicating whether approximation leads to significant relative error
//        int outputflag = 0;
//
//        //calculating strains in new configuration
//        for(int i=0; i<6*2+BEAM3K_COLLOCATION_POINTS; i++) //for all dof
//        {
//          Epetra_SerialDenseVector force_aux;
//          force_aux.Size(6*2+BEAM3K_COLLOCATION_POINTS);
//
//          //create new displacement and velocity vectors in order to store artificially modified displacements
//          std::vector<double> vel_aux(myvel);
//          std::vector<double> disp_aux(mydisp);
//
//          //modifying displacement artificially (for numerical derivative of internal forces):
//          disp_aux[i] += h_rel;
//          vel_aux[i] += h_rel / params.get<double>("delta time",0.01);
//
//          if(weakkirchhoff_)
//            CalculateInternalForcesWK(params,disp_aux,NULL,NULL,&force_aux,NULL,false);
//          else
//            CalculateInternalForcesSK(params,disp_aux,NULL,NULL,&force_aux,NULL,false);
//
//          //computing derivative d(fint)/du numerically by finite difference
//          for(int u = 0 ; u < 6*2+BEAM3K_COLLOCATION_POINTS ; u++ )
//          {
//            stiff_approx(u,i)= ( force_aux[u] - elevec1(u) )/ h_rel ;
//          }
//        } //for(int i=0; i<3; i++) //for all dof
//
//
//        for(int line=0; line<6*2+BEAM3K_COLLOCATION_POINTS; line++)
//        {
//          for(int col=0; col<6*2+BEAM3K_COLLOCATION_POINTS; col++)
//          {
//            if (fabs(elemat1(line,col)) > 1.0e-10)
//              stiff_relerr(line,col)= fabs( ( elemat1(line,col) - stiff_approx(line,col) )/ elemat1(line,col) );
//            else if (fabs(stiff_approx(line,col)) < 1.0e-5)
//              stiff_relerr(line,col)=0.0;
//            else
//              stiff_relerr(line,col)=1000.0;
//
//            //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
//            if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col))) //isnan = is not a number
//              stiff_relerr(line,col) = 0;
//
//            //if ( stiff_relerr(line,col) > 0)
//              outputflag = 1;
//          } //for(int col=0; col<3*nnode; col++)
//
//        } //for(int line=0; line<3*nnode; line++)
//
//        if(outputflag ==1)
//        {
//
//          std::cout<<"\n\n acutally calculated stiffness matrix\n";
//          for(int line=0; line<6*2+BEAM3K_COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+BEAM3K_COLLOCATION_POINTS; col++)
//            {
//              if(isnan(elemat1(line,col)))
//                std::cout<<"     nan   ";
//              else if(elemat1(line,col) == 0)
//                std::cout<<"     0     ";
//              else if(elemat1(line,col) >= 0)
//                std::cout<<"  "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
//              else
//                std::cout<<" "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
//            }
//            std::cout<<"\n";
//          }
//
//          std::cout<<"\n\n approximated stiffness matrix\n";
//          for(int line=0; line<6*2+BEAM3K_COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+BEAM3K_COLLOCATION_POINTS; col++)
//            {
//              if(isnan(stiff_approx(line,col)))
//                std::cout<<"     nan   ";
//              else if(stiff_approx(line,col) == 0)
//                std::cout<<"     0     ";
//              else if(stiff_approx(line,col) >= 0)
//                std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
//              else
//                std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
//            }
//            std::cout<<"\n";
//          }
//
//          std::cout<<"\n\n rel error stiffness matrix\n";
//          for(int line=0; line<6*2+BEAM3K_COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+BEAM3K_COLLOCATION_POINTS; col++)
//            {
//              if(isnan(stiff_relerr(line,col)))
//                std::cout<<"     nan   ";
//              else if(stiff_relerr(line,col) == 0)
//                std::cout<<"     0     ";
//              else if(stiff_relerr(line,col) >= 0)
//                std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
//              else
//                std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
//            }
//            std::cout<<"\n";
//          }
//        }
//
//      } //end of section in which numerical approximation for stiffness matrix is computed
//      //*******************************End: FD-CHECK************************************************************
