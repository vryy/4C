/*!----------------------------------------------------------------------
\file beam3wk_evaluate.cpp

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/


#include "beam3wk.H"
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

#include <Teuchos_TimeMonitor.hpp>


/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 10/12|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3wk::Evaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        std::vector<int>& lm,
                                        Epetra_SerialDenseMatrix& elemat1, //stiffness matrix
                                        Epetra_SerialDenseMatrix& elemat2, //mass matrix
                                        Epetra_SerialDenseVector& elevec1, //internal forces
                                        Epetra_SerialDenseVector& elevec2, //inertia forces
                                        Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam3wk::ActionType act = Beam3wk::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");

  if     (action == "calc_none")         dserror("No action supplied");
  else if (action=="calc_struct_linstiff")     act = Beam3wk::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")     act = Beam3wk::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3wk::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")   act = Beam3wk::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")   act = Beam3wk::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3wk::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress")     act = Beam3wk::calc_struct_stress;
  else if (action=="calc_struct_eleload")     act = Beam3wk::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")     act = Beam3wk::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam3wk::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   act = Beam3wk::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")    act = Beam3wk::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy")     act = Beam3wk::calc_struct_energy;
  else     dserror("Unknown type of action for Beam3wk");

  std::string test = params.get<std::string>("action","calc_none");

  switch(act)
  {
    case Beam3wk::calc_struct_ptcstiff:
    {
      dserror("no ptc implemented for Beam3wk element");
    }
    break;

    case Beam3wk::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3wk::calc_struct_nlnstiffmass:
    case Beam3wk::calc_struct_nlnstifflmass:
    case Beam3wk::calc_struct_nlnstiff:
    case Beam3wk::calc_struct_internalforce:
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
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      if (act == Beam3wk::calc_struct_nlnstiffmass)
      {
        CalculateInternalForces(params,mydisp,&elemat1,&elemat2,&elevec1,&elevec2,false);
      }
      else if (act == Beam3wk::calc_struct_nlnstifflmass)
      {
        dserror("The action calc_struct_nlnstifflmass is not implemented yet!");
      }
      else if (act == Beam3wk::calc_struct_nlnstiff)
      {
        CalculateInternalForces(params,mydisp,&elemat1,NULL,&elevec1,NULL,false);
      }
      else if (act == Beam3wk::calc_struct_internalforce)
      {
        CalculateInternalForces(params,mydisp,NULL,NULL,&elevec1,NULL,false);
      }

//      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
//      //correctly or not by means of a numerically approximated stiffness matrix
//      //The code block will work for all higher order elements.
//      if(Id() == 0) //limiting the following tests to certain element numbers
//      {
//
//        //variable to store numerically approximated stiffness matrix
//        Epetra_SerialDenseMatrix stiff_approx;
//        stiff_approx.Shape(6*2+COLLOCATION_POINTS,6*2+COLLOCATION_POINTS);
//
//
//        //relative error of numerically approximated stiffness matrix
//        Epetra_SerialDenseMatrix stiff_relerr;
//        stiff_relerr.Shape(6*2+COLLOCATION_POINTS,6*2+COLLOCATION_POINTS);
//
//        //characteristic length for numerical approximation of stiffness
//        double h_rel = 1e-7;
//
//        //flag indicating whether approximation leads to significant relative error
//        int outputflag = 0;
//
//        //calculating strains in new configuration
//        for(int i=0; i<6*2+COLLOCATION_POINTS; i++) //for all dof
//        {
//          Epetra_SerialDenseVector force_aux;
//          force_aux.Size(6*2+COLLOCATION_POINTS);
//
//          //create new displacement and velocity vectors in order to store artificially modified displacements
//          std::vector<double> vel_aux(myvel);
//          std::vector<double> disp_aux(mydisp);
//
//          //modifying displacement artificially (for numerical derivative of internal forces):
//          disp_aux[i] += h_rel;
//          vel_aux[i] += h_rel / params.get<double>("delta time",0.01);
//
//          CalculateInternalForces(params,disp_aux,NULL,NULL,&force_aux,NULL,false);
//
//          //computing derivative d(fint)/du numerically by finite difference
//          for(int u = 0 ; u < 6*2+COLLOCATION_POINTS ; u++ )
//          {
//            stiff_approx(u,i)= ( force_aux[u] - elevec1(u) )/ h_rel ;
//          }
//        } //for(int i=0; i<3; i++) //for all dof
//
//
//        for(int line=0; line<6*2+COLLOCATION_POINTS; line++)
//        {
//          for(int col=0; col<6*2+COLLOCATION_POINTS; col++)
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
//          for(int line=0; line<6*2+COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+COLLOCATION_POINTS; col++)
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
//          for(int line=0; line<6*2+COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+COLLOCATION_POINTS; col++)
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
//          for(int line=0; line<6*2+COLLOCATION_POINTS; line++)
//          {
//            for(int col=0; col<6*2+COLLOCATION_POINTS; col++)
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
    }
    break;

    case calc_struct_energy:
    {
      dserror("No energy output implemented for beam3wk elements");
    }
    break;

    case calc_struct_stress:
    {
      dserror("No stress output implemented for beam3wk elements");
    }
    break;

    case calc_struct_update_istep:
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

      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get element velocities
      std::vector<double> myvel(lm.size());
      std::vector<double> myacc(lm.size());

      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        Teuchos::RCP<const Epetra_Vector> acc  = discretization.GetState("acceleration");
        if (acc==Teuchos::null) dserror("Cannot get state vectors 'acceleration'");
        DRT::UTILS::ExtractMyValues(*acc,myacc,lm);
      }
      UpdateTriads(params,myvel,mydisp);

    }
    break;

    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitable timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
       dserror("The action calc_struct_reset_istep is not implemented yet!");
    }
    break;

    default:
      dserror("Unknown type of action for Beam3wk %d", act);
     break;

  }//switch(act)

  return 0;

}  //DRT::ELEMENTS::Beam3wk::Evaluate

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                           meier 05/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  dserror("Lumped mass matrix not implemented yet!!!");
}

/*------------------------------------------------------------------------------------------------------------*
 | stiffness matrix            (private)                                                           meier 09/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::CalculateInternalForces( Teuchos::ParameterList& params,
                                                      std::vector<double>&      disp,
                                                      Epetra_SerialDenseMatrix* stiffmatrix,
                                                      Epetra_SerialDenseMatrix* massmatrix,
                                                      Epetra_SerialDenseVector* force,
                                                      Epetra_SerialDenseVector* inertia_force,
                                                      bool update)
{
  TEUCHOS_FUNC_TIME_MONITOR("Beam3wk::CalculateInternalForces");
  if(COLLOCATION_POINTS!=2 and COLLOCATION_POINTS!=3 and COLLOCATION_POINTS!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
  }

  //number of nodes fixed for these element
  const int nnode = 2;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_int(true);
  LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_int_aux(true);

  //vector for current nodal DoFs in total Lagrangian style, i.e. displacement + initial values:
  //rotvec_==true:  disp_totlag=[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3]
  //rotvec_==false: disp_totlag=[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3]
  //! The Number of collocation points can take on the values 2, 3 and 4. 3 and 4 are interior nodes.
  //This leads e.g. in the case rotvec_==true to the following ordering:
  //if COLLOCATION_POINTS = 2: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2]
  //if COLLOCATION_POINTS = 3: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3]
  //if COLLOCATION_POINTS = 4: [ \v{d}_1 \v{theta}_1 t_1 \v{d}_2 \v{theta}_1 t_2 \alpha_3 \alpha_4]
  std::vector<FAD> disp_totlag(6*nnode+COLLOCATION_POINTS, 0.0);

  //vector containing locally assembled nodal positions and tangents required for centerline:
  //r(s)=N(s)*disp_totlag_centerline, with disp_totlag_centerline=[\v{d}_1, \v{t}_1, 0, \v{d}_2, \v{t}_2, 0, 0]
  std::vector<FAD> disp_totlag_centerline(6*nnode+COLLOCATION_POINTS, 0.0);

  //CP values of strains and their variations needed for interpolation
  std::vector<FAD> epsilon_cp(COLLOCATION_POINTS); // axial tension epsilon=|r_s|-1 at collocation points
  std::vector<LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> > v_epsilon_cp(COLLOCATION_POINTS);
  std::vector<LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> > v_thetaperp_cp(COLLOCATION_POINTS);
  std::vector<LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> > v_thetapar_cp(COLLOCATION_POINTS);

  //interpolated values of strains and their variations evaluated at Gauss points
  FAD epsilon;
  LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> v_epsilon(true);
  LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> v_thetaperp_s(true);
  LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> v_thetapar_s(true);

  //Further material and spatial strains and forces to be evaluated at the Gauss points
  LINALG::TMatrix<FAD,3,1> kappa; //material curvature
  LINALG::TMatrix<FAD,3,1> m; //spatial moment stress resultant
  LINALG::TMatrix<FAD,3,1> M; //material moment stress resultant
  FAD f_par; //material=spatial axial force component

  //Triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_ref_cp(COLLOCATION_POINTS);  //reference triads (SR-system) at collocation points
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(COLLOCATION_POINTS);  //material triads at collocation points
  std::vector<LINALG::TMatrix<FAD,3,1> > theta_cp(COLLOCATION_POINTS);  //angle at collocation points

  //Interpolated material triad and angles evaluated at Gauss point
  LINALG::TMatrix<FAD,3,3> triad_mat; //material triad at gp
  LINALG::TMatrix<FAD,3,1> theta;    //interpolated angle theta
  LINALG::TMatrix<FAD,3,1> theta_s;  //derivative of theta with respect to arc-length s

  //matrices holding the assembled shape functions and s-derivatives
  LINALG::TMatrix<FAD,3,6*nnode+COLLOCATION_POINTS> N_s;
  LINALG::TMatrix<FAD,1,6*nnode+COLLOCATION_POINTS> L;
  LINALG::TMatrix<FAD,1,6*nnode+COLLOCATION_POINTS> L_s;

  //Matrices for individual shape functions and xi-derivatives
  LINALG::TMatrix<FAD,1,2*nnode> N_i_xi;
  LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> L_i;
  LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> L_i_xi;

  //Matrices for individual s-derivatives
  LINALG::TMatrix<FAD,1,2*nnode> N_i_s;
  LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> L_i_s;

  //Additional kinematic quantities
  LINALG::TMatrix<FAD,3,1> r_s; //Matrix to store r'
  FAD abs_r_s; // ||r'||

  //MISC
  double xi=0.0;  //parameter coordinated
  int ind=0; //position index where CP quantities have to be stored

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3wk);

  //Set current positions and orientations at all nodes:
  UpdateDispTotlag(disp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  //********begin: evaluate quantities at collocation points********************************
  for(int node=0; node<COLLOCATION_POINTS; node++)
  {
    //calculate xi of cp
    //node=0->xi=-1  node=1->xi=0  node=2->xi=1
    xi=(double)node/(COLLOCATION_POINTS-1)*2-1.0;

    //get value of interpolating function of theta (lagrange polynomials) at xi
    L_i.Clear();
    N_i_xi.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

    //Determine storage position for the node node
    ind=LARGEROTATIONS::NumberingTrafo(node+1, COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctions(N_i_xi,jacobi_cp_[ind],N_s);
    r_s.Clear();
    //Calculation of r' at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+COLLOCATION_POINTS; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];
      }
    }
    //calculate epsilon at collocation point
    abs_r_s=Norm<FAD>(r_s);
    epsilon_cp[ind]=abs_r_s-1.0;

    AssembleShapefunctions(L_i,L);

    v_epsilon_cp[ind].Clear();
    v_epsilon_cp[ind].MultiplyTN(N_s,r_s);
    v_epsilon_cp[ind].Scale(1.0/abs_r_s);

    v_thetaperp_cp[ind].Clear();
    LINALG::TMatrix<FAD,3,3> auxmatrix(true);
    LARGEROTATIONS::computespin<FAD>(auxmatrix,r_s);
    v_thetaperp_cp[ind].MultiplyTN(N_s,auxmatrix);
    v_thetaperp_cp[ind].Scale(-1.0/(abs_r_s*abs_r_s));

    v_thetapar_cp[ind].Clear();
    for (int i=0;i<6*nnode+COLLOCATION_POINTS;i++)
      for (int j=0;j<3;j++)
        (v_thetapar_cp[ind])(i,j)=L(i)*r_s(j)/abs_r_s;
  } //for (int node=0; node<COLLOCATION_POINTS; node++)

  //calculate angle at cp (this has to be done in a SEPARATE loop as follows)
  for(int node=0; node<COLLOCATION_POINTS; node++)
  {
    theta_cp[node].Clear();
    triadtoangle(theta_cp[node],triad_mat_cp[REFERENCE_NODE],triad_mat_cp[node]);
  }
  //********end: evaluate quantities at collocation points********************************

  //******begin: gauss integration for internal force vector and stiffness matrix*********
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    L_i.Clear();
    L_i_xi.Clear();
    L_i_s.Clear();

    DRT::UTILS::shape_function_1D(L_i,xi,Shape());
    DRT::UTILS::shape_function_1D_deriv1(L_i_xi,xi,Shape());
    L_i_s.Update(1.0/jacobi_[numgp],L_i_xi,0.0);

    //Calculate collocation piont interpolations
    v_epsilon.Clear();
    v_thetaperp_s.Clear();
    v_thetapar_s.Clear();
    epsilon=0.0;
    theta.Clear();
    theta_s.Clear();

    for(int node=0; node<COLLOCATION_POINTS; node++)
    {
      v_epsilon.Update(L_i(node),v_epsilon_cp[node],1.0);
      v_thetaperp_s.Update(L_i_s(node),v_thetaperp_cp[node],1.0);
      v_thetapar_s.Update(L_i_s(node),v_thetapar_cp[node],1.0);
      epsilon+=L_i(node)*epsilon_cp[node];
      theta.Update(L_i(node),theta_cp[node],1.0);
      theta_s.Update(L_i_s(node),theta_cp[node],1.0);
    }

    //compute material strain resultant kappa
    kappa.Clear();
    computestrain(theta,theta_s,kappa);

    for(int i=0; i<3; i++)
    {
      kappa(i)-=kappa0_[numgp](i);
    }

    //compute material stress resultants
    M.Clear();
    f_par=0.0;
    straintostress(kappa,epsilon,M,f_par);

    //compute material triad at gp
    triad_mat.Clear();
    angletotriad(theta,triad_mat_cp[REFERENCE_NODE],triad_mat);

    //pushforward of stress resultants
    m.Clear();
    m.Multiply(triad_mat,M);

    f_int_aux.Clear();
    f_int_aux.Multiply(v_thetaperp_s,m);
    f_int_aux.Scale(wgt*jacobi_[numgp]);
    f_int.Update(1.0,f_int_aux,1.0);

    f_int_aux.Clear();
    f_int_aux.Multiply(v_thetapar_s,m);
    f_int_aux.Scale(wgt*jacobi_[numgp]);
    f_int.Update(1.0,f_int_aux,1.0);

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
    for(int i = 0; i < 6*nnode+COLLOCATION_POINTS; i++)
    {
      for(int j = 0; j < 6*nnode+COLLOCATION_POINTS; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }
  }
  if(force!=NULL)
  {
    for (int i=0; i< 6*nnode+COLLOCATION_POINTS; i++)
    {
      (*force)(i)=f_int(i).val();
    }
  }

  //****************Update of the old reference triads at the nodes***********************************************
  if (update and rotvec_==false)
  {
      UpdateRefTriads(disp_totlag_centerline);
  }
  //****************end: Update of the old reference triads at the nodes********************************************

  //calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
  //finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
  //needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
  //steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
  //consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
  //a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
  //our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
  //according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)

  double beta = params.get<double>("rot_beta",1000);

  if (massmatrix != NULL and inertia_force != NULL)
  {
    if(beta < 999)
    {
      double gamma = params.get<double>("rot_gamma",1000);
      double alpha_f = params.get<double>("rot_alphaf",1000);
      double alpha_m = params.get<double>("rot_alpham",1000);
      double dt = params.get<double>("delta time",1000);

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

      double scaledcrosssec=inertscaletrans_*crosssec_;

      LINALG::TMatrix<FAD,3,6*nnode+COLLOCATION_POINTS> N(true);
      LINALG::TMatrix<FAD,1,2*nnode> N_i(true);
      LINALG::TMatrix<FAD,3,1> rnewmass(true); //Matrix to store r
      LINALG::TMatrix<FAD,3,3> triad_mat_old(true);

      //internal force vector
      LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_inert(true);
      LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_inert_aux(true);
      LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> v_thetaperp(true);
      LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,3> v_thetapar(true);

      for (int numgp=0; numgp<gausspoints.nquad; numgp++)//loop through Gauss points
      {
        //Get location and weight of GP in parameter space
        const double xi = gausspoints.qxg[numgp][0];
        const double wgt = gausspoints.qwgt[numgp];

        L_i.Clear();
        DRT::UTILS::shape_function_1D(L_i,xi,Shape());
        v_thetaperp.Clear();
        v_thetapar.Clear();
        theta.Clear();

        for(int node=0; node<COLLOCATION_POINTS; node++)
        {
          v_thetaperp.Update(L_i(node),v_thetaperp_cp[node],1.0);
          v_thetapar.Update(L_i(node),v_thetapar_cp[node],1.0);
          theta.Update(L_i(node),theta_cp[node],1.0);
        }

        //compute material triad at gp
        triad_mat.Clear();
        angletotriad(theta,triad_mat_cp[REFERENCE_NODE],triad_mat);
        LINALG::TMatrix<FAD,4,1> Qnewmass(true);
        triadtoquaternion(triad_mat,Qnewmass);

        N_i.Clear();
        N.Clear();
        DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
        AssembleShapefunctionsN(N_i,N);
        rnewmass.Clear();
        //Calculation of r at xi
        for (int i=0; i<3; i++)
        {
          for (int j=0; j<6*nnode+COLLOCATION_POINTS; j++)
          {
            rnewmass(i)+=N(i,j)*disp_totlag_centerline[j];
          }
        }

        triad_mat_old.Clear();
        LINALG::TMatrix<FAD,4,1>  Qconv(true);
        for(int i=0;i<4;i++)
          Qconv(i)=(Qconvmass_[numgp])(i);

        quaterniontotriad(Qconv,triad_mat_old);

        LINALG::TMatrix<FAD,3,3>  deltatriad(true);
        deltatriad.MultiplyNT(triad_mat,triad_mat_old);
        LINALG::TMatrix<FAD,4,1>  deltaQ(true);
        triadtoquaternion(deltatriad,deltaQ);
        LINALG::TMatrix<FAD,3,1> deltatheta(true);
        quaterniontoangle(deltaQ,deltatheta);

        //compute material counterparts of spatial vectors
        LINALG::TMatrix<FAD,3,1> deltaTHETA(true);
        LINALG::TMatrix<FAD,3,1> Wconvmass(true);
        LINALG::TMatrix<FAD,3,1> Aconvmass(true);
        LINALG::TMatrix<FAD,3,1> Amodconvmass(true);

        deltaTHETA.MultiplyTN(triad_mat,deltatheta);

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

        Pi_t.Multiply(triad_mat,auxvector1);
        LINALG::TMatrix<FAD,3,1> L_t(true);
        L_t.Update(rho*scaledcrosssec,rttnewmass,1.0);

        f_inert_aux.Clear();
        f_inert_aux.Multiply(v_thetaperp,Pi_t);
        f_inert_aux.Scale(wgt*jacobi_[numgp]);
        f_inert.Update(1.0,f_inert_aux,1.0);

        f_inert_aux.Clear();
        f_inert_aux.Multiply(v_thetapar,Pi_t);
        f_inert_aux.Scale(wgt*jacobi_[numgp]);
        f_inert.Update(1.0,f_inert_aux,1.0);

        f_inert_aux.Clear();
        f_inert_aux.MultiplyTN(N,L_t);
        f_inert_aux.Scale(wgt*jacobi_[numgp]);
        f_inert.Update(1.0,f_inert_aux,1.0);

        //**********begin: update class variables needed for storage**************
        LINALG::TMatrix<FAD,3,1> wnewmass(true);
        LINALG::TMatrix<FAD,3,1> anewmass(true);
        LINALG::TMatrix<FAD,3,1> amodnewmass(true);
        wnewmass.Multiply(triad_mat,Wnewmass);
        anewmass.Multiply(triad_mat,Anewmass);
        amodnewmass.Multiply(triad_mat,Amodnewmass);

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

      if(rotvec_==true)
      {
        ApplyRotVecTrafo(disp_totlag_centerline,f_inert);
      }

      //Update mass matrix and inertia force vector
      if(massmatrix!=NULL)
      {
        //Calculating stiffness matrix with FAD
        for(int i = 0; i < 6*nnode+COLLOCATION_POINTS; i++)
        {
          for(int j = 0; j < 6*nnode+COLLOCATION_POINTS; j++)
          {
            (*massmatrix)(i,j)=f_inert(i).dx(j);
          }
        }
      }
      for (int i=0; i< 6*nnode+COLLOCATION_POINTS; i++)
      {
        (*inertia_force)(i)=f_inert(i).val();
      }
    }
    else
    {
      //This is a dummy mass matrix which is necessary for statmech simulations
      for (int i=0; i<6*nnode; i++)
        (*massmatrix)(i,i) = 1;
    }//if(beta < 999)
  }//if (massmatrix != NULL)

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition (public)                                  meier 09/12|
 *-----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3wk::EvaluateNeumann(Teuchos::ParameterList& params,
                                               DRT::Discretization& discretization,
                                               DRT::Condition& condition,
                                               std::vector<int>& lm,
                                               Epetra_SerialDenseVector& elevec1,
                                               Epetra_SerialDenseMatrix* elemat1)
{

  const int twistdofs = COLLOCATION_POINTS;
  if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
  }
  //dimensions of freedom per node
  const int nnode=2;

  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

  std::vector<FAD> disp_totlag(6*nnode+COLLOCATION_POINTS, 0.0);
  std::vector<FAD> disp_totlag_centerline(6*nnode+COLLOCATION_POINTS, 0.0);
  std::vector<LINALG::TMatrix<FAD,3,3> > triad_mat_cp(COLLOCATION_POINTS);

  UpdateDispTotlag(mydisp, disp_totlag);
  SetNodalVariables(disp_totlag,disp_totlag_centerline,triad_mat_cp);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
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
    LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_ext(true);
    LINALG::TMatrix<FAD,6*nnode+COLLOCATION_POINTS,1> f_ext_aux(true);
    std::vector< LINALG::Matrix<3,3> > Gref(2);
    for (int node=0;node<2;node++)
    {
      Gref[node].Clear();
      LARGEROTATIONS::angletotriad(theta0_[node],Gref[node]);
    }

    // gaussian points
    DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussrulebeam3wk);
    LINALG::TMatrix<FAD,1,4> N_i;
    LINALG::TMatrix<FAD,3,6*nnode + COLLOCATION_POINTS> N;

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
          dserror("Line Neumann conditions for distributed moments are not implemented for beam3wk so far! Only the function flag 1, 2 and 3 can be set!");
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
      for(int i = 0; i < 6*nnode+COLLOCATION_POINTS; i++)
      {
        for(int j = 0; j < 6*nnode+COLLOCATION_POINTS; j++)
        {
          (*elemat1)(i,j)=-f_ext(i).dx(j);
        }
      }
    }
    if(elevec1!=NULL)
    {
      for (int i=0; i< 6*nnode+COLLOCATION_POINTS; i++)
      {
        elevec1(i)=f_ext(i).val();
      }
    }
  }//if a line neumann condition needs to be linearized

//  std::cout << "elemat1: " << std::endl;
//  if(elemat1!=NULL)
//  {
//    //Calculating stiffness matrix with FAD
//    for(int i = 0; i < 6*nnode+COLLOCATION_POINTS; i++)
//    {
//      for(int j = 0; j < 6*nnode+COLLOCATION_POINTS; j++)
//      {
//        std::cout << (*elemat1)(i,j) << "  ";
//      }
//      std::cout << std::endl;
//    }
//  }
//  std::cout << "fext: " << std::endl;
//  {
//    for (int i=0; i< 6*nnode+COLLOCATION_POINTS; i++)
//    {
//      std::cout << elevec1(i) << std::endl;
//    }
//  }

  return 0;
}  //DRT::ELEMENTS::Beam3wk::EvaluateNeumann

/*-----------------------------------------------------------------------------------------------------------*
 |  Rotate a triad around the first base vector (tangent)                                         meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::RotateTriad(LINALG::TMatrix<FAD,3,3> triad, FAD alpha, LINALG::TMatrix<FAD,3,3>& triad_rot)
{

  for (int i=0; i<3; i++)
  {
  triad_rot(i,0)=triad(i,0);
  triad_rot(i,1)=triad(i,1)*cos(alpha)+triad(i,2)*sin(alpha);
  triad_rot(i,2)=triad(i,2)*cos(alpha)-triad(i,1)*sin(alpha);
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble all shape functions                                                                   meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i_xi,
                                                     LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> C_i,
                                                     LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> C_i_xi,
                                                     FAD jacobi_local,
                                                     LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N_s,
                                                     LINALG::TMatrix<FAD,1,2*6+COLLOCATION_POINTS>& C,
                                                     LINALG::TMatrix<FAD,1,2*6+COLLOCATION_POINTS>& C_s)
{

  //Matrices for N_i,s and N_i,ss
  LINALG::TMatrix<FAD,1,4> N_i_s;
  LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> C_i_s;

  N_i_s.Clear();
  C_i_s.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1/jacobi_local);
  C_i_s=C_i_xi;
  C_i_s.Scale(1/jacobi_local);

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C


  #if defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==2)

  int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2};

  #elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==3)

  int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3};

  #elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==4)

  int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                         {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                         {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3,4};

  #else
  dserror("COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+COLLOCATION_POINTS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N_s(j,i)=0;
      }
      else
      {
        N_s(j,i)=N_i_s(assembly_N[j][i]-1);
      }
    }
    if(assembly_C[i]==0)
    {
      C(i)=0.0;
      C_s(i)=0.0;
    }
    else
    {
        C(i)=C_i(assembly_C[i]-1);
        C_s(i)=C_i_s(assembly_C[i]-1);
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble C shape function                                                                meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::AssembleShapefunctions(LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> C_i,
                                                    LINALG::TMatrix<FAD,1,2*6+COLLOCATION_POINTS>& C)
{



  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C

  #if defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==2)

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2};

  #elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==3)

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3};

  #elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==4)

  int assembly_C[2*6+COLLOCATION_POINTS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3,4};

  #else
  dserror("COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+COLLOCATION_POINTS; i++)
  {
    if(assembly_C[i]==0)
    {
      C(i)=0.0;
    }
    else
    {
        C(i)=C_i(assembly_C[i]-1);
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble the N_s shape functions                                                              meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i_xi,
                                                     FAD jacobi_local,
                                                     LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N_s)
{

  LINALG::TMatrix<FAD,1,4> N_i_s;

  N_i_s.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;

  N_i_s.Scale(1/jacobi_local);

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C

#if defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==2)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

#elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==3)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

#elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==4)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                       {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                       {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

#else
dserror("COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
#endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+COLLOCATION_POINTS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N_s(j,i)=0;
      }
      else
      {
        N_s(j,i)=N_i_s(assembly_N[j][i]-1);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble the N shape functions                                                              meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::AssembleShapefunctionsN( LINALG::TMatrix<FAD,1,4> N_i,
                                                     LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N)
{

#if defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==2)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

#elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==3)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

#elif defined(COLLOCATION_POINTS) && (COLLOCATION_POINTS==4)

int assembly_N[3][2*6+COLLOCATION_POINTS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                       {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                       {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

#else
dserror("COLLOCATION_POINTS has to be defined. Only the values 2,3 and 4 are valid for COLLOCATION_POINTS!!!");
#endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+COLLOCATION_POINTS; i++)
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
 |  Calculate position vectors from displacement +  initial position                               meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::UpdateDispTotlag(std::vector<double> disp, std::vector<FAD>& disp_totlag)
{
  for (int dof=0;dof<2*6+COLLOCATION_POINTS;dof++)
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
          disp_totlag[7*node+ndof]+=(GetNodalRefTangent(node))(ndof-3);
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
      for (int ndof=0;ndof<7;ndof++)//loop over dofs per node
      {
        if(ndof<3)
        {
          disp_totlag[7*node+ndof]+=Nodes()[node]->X()[ndof];
        }
        else if (ndof<6)
        {
          disp_totlag[7*node+ndof]+=(theta0_[node])(ndof-3);
        }
        else
        {
          //here we have to add the initial length of the tangents at the boundary nodes, i.e. ||r'_i(t=0)||=1:
          disp_totlag[7*node+ndof]+=1;
        }
      }
    }
  }

  //Next, we have to set variables for FAD
  for (int dof=0;dof<2*6+COLLOCATION_POINTS;dof++)
    disp_totlag[dof].diff(dof,2*6+COLLOCATION_POINTS);

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Set positions vectors and tangents at boundary nodes and triads at all CPs                    meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::SetNodalVariables( std::vector<FAD>& disp_totlag,
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
      LARGEROTATIONS::quaterniontotriad(qrefconv_[node],triad_aux);
      CalculateSRTriads<FAD>(tangent,triad_aux,triad_ref);

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
  LINALG::TMatrix<FAD,3,6*2+COLLOCATION_POINTS> N_s(true);
  LINALG::TMatrix<FAD,1,COLLOCATION_POINTS> L_i;
  double xi=0.0;
  int ind=0;
  //********begin: evaluate quantities at collocation points********************************
  for(int node=1; node<COLLOCATION_POINTS-1; node++)
  {
    //calculate xi of cp
    //node=0->xi=-1  node=1->xi=0 node=2->xi=1
    xi=(double)node/(COLLOCATION_POINTS-1)*2-1.0;
    N_i_xi.Clear();
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    L_i.Clear();
    DRT::UTILS::shape_function_1D(L_i,xi,Shape());

    //Determine storage position for the node node
    ind=LARGEROTATIONS::NumberingTrafo(node+1, COLLOCATION_POINTS);

    N_s.Clear();
    AssembleShapefunctions(N_i_xi,jacobi_cp_[ind],N_s);
    tangent.Clear();
    //Calculation of r' at xi
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*2+COLLOCATION_POINTS; j++)
      {
        tangent(i)+=N_s(i,j)*disp_totlag_centerline[j];
      }
    }

    alpha=disp_totlag[7*2+ind-2];

    //calculate new sr triads
    triad_ref.Clear();
    triad_aux.Clear();
    LARGEROTATIONS::quaterniontotriad(qrefconv_[ind],triad_aux);
    CalculateSRTriads<FAD>(tangent,triad_aux,triad_ref);

    //calculate material triad
    triad_mat_cp[ind].Clear();
    RotateTriad(triad_ref,alpha,triad_mat_cp[ind]);

  } //for (int node=0; node<COLLOCATION_POINTS; node++)

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Pre-multiply trafo matrix if rotvec_==true: \tilde{\vec{f}_int}=\mat{T}^T*\vec{f}_int         meier 01/16|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::ApplyRotVecTrafo( std::vector<FAD>& disp_totlag_centerline,
                                               LINALG::TMatrix<FAD,6*2+COLLOCATION_POINTS,1>& f_int)
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

    t=Norm<FAD>(g_1);;
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
 |  Calculate tangent at an arbitrary point                                                       meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::Calculate_r_s(std::vector<FAD> disp_totlag_centerline, double xi, double jacobi_local, LINALG::TMatrix<FAD,3,1>& r_s)
{

  FAD abs_r_s=0.0;
  LINALG::TMatrix<FAD,1,4> N_i_xi;
  LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS> N_s;
  N_i_xi.Clear();
  N_s.Clear();

  //Get hermite derivatives N'xi and N''xi
  DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

  AssembleShapefunctions(N_i_xi, jacobi_local, N_s);

  //Calculation of r' at xi
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<2*6+COLLOCATION_POINTS; j++)
    {
      r_s(i)+=N_s(i,j)*disp_totlag_centerline[j];
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Update the reference triad at the element nodes at the end of time step                       meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::UpdateRefTriads(std::vector<FAD> disp_totlag_centerline)
{
  double xi=0.0;
  int ind=0;
  LINALG::Matrix<3,1> r_xi(true);
  LINALG::Matrix<3,1> t(true);
  LINALG::Matrix<3,3> triad(true);
  LINALG::Matrix<3,3> triad_aux(true);

  LINALG::Matrix<1,COLLOCATION_POINTS> L_i(true);
  LINALG::Matrix<1,4> N_i_xi(true);

  for(int node=0; node<COLLOCATION_POINTS; node++)
  {
      //calculate xi of cp
      xi=(double)node/(COLLOCATION_POINTS-1)*2-1;

      //get value of interpolating function of theta (lagrange polynomials) at xi
      L_i.Clear();
      DRT::UTILS::shape_function_1D(L_i,xi,Shape());

      //Determine storage position for the node node
      ind=LARGEROTATIONS::NumberingTrafo(node+1, COLLOCATION_POINTS);

      N_i_xi.Clear();
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

      //initialize
      r_xi.Clear();
      t.Clear();
      for (int i=0; i<3; i++)
      {
        r_xi(i)+=disp_totlag_centerline[i].val()*N_i_xi(0)+disp_totlag_centerline[7+i].val()*N_i_xi(2)+disp_totlag_centerline[3+i].val()*N_i_xi(1)+disp_totlag_centerline[10+i].val()*N_i_xi(3);
      }
      t=r_xi;
      t.Scale(1.0/r_xi.Norm2());

      //calculate new sr triads
      triad.Clear();
      triad_aux.Clear();
      LARGEROTATIONS::quaterniontotriad(qrefconv_[ind],triad_aux);
      CalculateSRTriads<double>(t,triad_aux,triad);

      //extract quaternion
      LARGEROTATIONS::triadtoquaternion(triad,qrefconv_[ind]);
  }

  return;
}

void DRT::ELEMENTS::Beam3wk::computestrain(LINALG::Matrix<3,1>& theta, LINALG::Matrix<3,1>& theta_deriv, LINALG::Matrix<3,1>& kappa)
{
  LINALG::Matrix<3,3> Tinv;
  kappa.Clear();
  Tinv.Clear();

  Tinv=LARGEROTATIONS::Tinvmatrix(theta);
  kappa.MultiplyTN(Tinv,theta_deriv);
}

void DRT::ELEMENTS::Beam3wk::computestrain(LINALG::TMatrix<FAD,3,1>& theta, LINALG::TMatrix<FAD,3,1>& theta_deriv, LINALG::TMatrix<FAD,3,1>& kappa)
{
  LINALG::TMatrix<FAD,3,3> Tinv;
  kappa.Clear();
  Tinv.Clear();

  Tinv=Tinvmatrix(theta);
  kappa.MultiplyTN(Tinv,theta_deriv);
}

void DRT::ELEMENTS::Beam3wk::straintostress(LINALG::TMatrix<FAD,3,1>& kappa, FAD epsilon, LINALG::TMatrix<FAD,3,1>& M, FAD& f_par)
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
  M(0)=sm*Irr_*kappa(0);
  M(1)=ym*Iyy_*kappa(1);
  M(2)=ym*Izz_*kappa(2);

}
/*----------------------------------------------------------------------------------------------------------*
 |  Update of the nodal triads at the end of time step (public)                                  meier 05/13|
 *----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::UpdateTriads( Teuchos::ParameterList& params,
                                              std::vector<double>& vel,
                                              std::vector<double>& disp)
{
  CalculateInternalForces(params,disp,NULL,NULL,NULL,NULL,true);

  return;
}
///*--------------------------------------------------------------------------------------------------------*
// |  Cross product of two Tmatrices (public)                                                    meier 05/13|
// *--------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3wk::CrossProduct(LINALG::TMatrix<FAD,3,1> T1,LINALG::TMatrix<FAD,3,1> T2,LINALG::TMatrix<FAD,3,1>& T3)
{
  T3(0,0) = (T1(1,0)*T2(2,0)) - (T1(2,0)*T2(1,0));
  T3(1,0) = (T1(2,0)*T2(0,0)) - (T1(0,0)*T2(2,0));
  T3(2,0) = (T1(0,0)*T2(1,0)) - (T1(1,0)*T2(0,0));
}

/*----------------------------------------------------------------------------------------------------------*
 | Get position vector at xi for given nodal displacements                                        popp 02/16|
 *----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3wk::GetPos(double& xi, LINALG::Matrix<12,1>& disp_totlag) const
{
  LINALG::Matrix<3,1> r(true);
  LINALG::Matrix<4,1> N_i(true);

  DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

  for (int n=0;n<4;n++)
  {
    for (int i=0;i<3;i++)
    {
      r(i)+=N_i(n)*disp_totlag(3*n+i);
    }
  }

  return (r);
}
