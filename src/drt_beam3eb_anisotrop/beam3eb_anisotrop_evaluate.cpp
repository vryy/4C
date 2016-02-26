/*!----------------------------------------------------------------------
\file beam3eb_anisotrop_evaluate.cpp

\brief three dimensional nonlinear rod based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3eb_anisotrop.H"
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
#include "../drt_beamcontact/beam3contact_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 meier 10/12|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebanisotrop::Evaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        std::vector<int>& lm,
                                        Epetra_SerialDenseMatrix& elemat1, //stiffness matrix
                                        Epetra_SerialDenseMatrix& elemat2, //mass matrix
                                        Epetra_SerialDenseVector& elevec1, //internal forces
                                        Epetra_SerialDenseVector& elevec2, //inertia forces
                                        Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam3ebanisotrop::ActionType act = Beam3ebanisotrop::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");

  if     (action == "calc_none")         dserror("No action supplied");
  else if (action=="calc_struct_linstiff")     act = Beam3ebanisotrop::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")     act = Beam3ebanisotrop::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3ebanisotrop::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")   act = Beam3ebanisotrop::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")   act = Beam3ebanisotrop::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3ebanisotrop::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress")     act = Beam3ebanisotrop::calc_struct_stress;
  else if (action=="calc_struct_eleload")     act = Beam3ebanisotrop::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")     act = Beam3ebanisotrop::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam3ebanisotrop::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   act = Beam3ebanisotrop::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")    act = Beam3ebanisotrop::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy")     act = Beam3ebanisotrop::calc_struct_energy;
  else     dserror("Unknown type of action for Beam3ebanisotrop");

  std::string test = params.get<std::string>("action","calc_none");

  switch(act)
  {
    case Beam3ebanisotrop::calc_struct_ptcstiff:
    {
      dserror("no ptc implemented for Beam3ebanisotrop element");
    }
    break;

    case Beam3ebanisotrop::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3ebanisotrop::calc_struct_nlnstiffmass:
    case Beam3ebanisotrop::calc_struct_nlnstifflmass:
    case Beam3ebanisotrop::calc_struct_nlnstiff:
    case Beam3ebanisotrop::calc_struct_internalforce:
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

      // get element velocities
      std::vector<double> myvel(lm.size());

      // get element accelerations
      std::vector<double> myacc(lm.size());

      if(act == Beam3ebanisotrop::calc_struct_nlnstiffmass or act == Beam3ebanisotrop::calc_struct_nlnstifflmass)
      {
        #if defined (VP) or defined (SR1)
          dserror("Dynamic Calculation only implemented for NSRISR formulation so far!!!");
        #endif

        #if defined (MATERIALREF)
          dserror("Dynamic Calculation is not implemented for MATERIALREF so far!!!");
        #endif

        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        Teuchos::RCP<const Epetra_Vector> acc  = discretization.GetState("acceleration");
        if (acc==Teuchos::null) dserror("Cannot get state vectors 'acceleration'");
        DRT::UTILS::ExtractMyValues(*acc,myacc,lm);
      }
      if (act == Beam3ebanisotrop::calc_struct_nlnstiffmass)
      {
        CalculateInternalForces(params,myvel,mydisp,&elemat1,&elevec1,false);
        CalculateInertialForces(params,myacc,myvel,mydisp,&elemat2,&elevec2,false);
      }
      else if (act == Beam3ebanisotrop::calc_struct_nlnstifflmass)
      {
        CalculateInternalForces(params,myvel,mydisp,&elemat1,&elevec1,false);
        CalculateInertialForces(params,myacc,myvel,mydisp,&elemat2,&elevec2,false);
        lumpmass(&elemat2);
      }
      else if (act == Beam3ebanisotrop::calc_struct_nlnstiff)
      {
        CalculateInternalForces(params,myvel,mydisp,&elemat1,&elevec1,false);
      }
      else if (act == Beam3ebanisotrop::calc_struct_internalforce)
      {
        CalculateInternalForces(params,myvel,mydisp,NULL,&elevec1,false);
      }
    }
    break;

    case calc_struct_energy:
    {
      CalculateEnergy(elevec1);
    }
    break;

    case calc_struct_stress:
    {
      dserror("No stress output implemented for beam3eb_anisotrop elements");
    }
    break;

    case calc_struct_update_istep:
    {
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
      if(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP")!=INPAR::STR::dyna_statics)
      {
      UpdateNodalAngularVelocities(&elevec1,&elevec2,&elevec3);
      }

      timestepcount_++;
    }
    break;

    case calc_struct_reset_istep:
      //not necessary since no class variables are modified in predicting steps
    break;

    default:
      dserror("Unknown type of action for Beam3ebanisotrop %d", act);
     break;

  }//switch(act)

  return 0;

}  //DRT::ELEMENTS::Beam3ebanisotrop::Evaluate

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                           meier 05/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  dserror("Lumped mass matrix not implemented yet!!!");
}

/*------------------------------------------------------------------------------------------------------------*
 | stiffness matrix            (private)                                                           meier 09/12|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateInternalForces(Teuchos::ParameterList& params,
                      std::vector<double>&           vel,
                      std::vector<double>&           disp,
                      Epetra_SerialDenseMatrix* stiffmatrix,
                      Epetra_SerialDenseVector* force,
                      bool update)
{

  const int twistdofs = TWISTDOFS;
  if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  }

  #ifdef BEAM3EBANISOTROPORTHOPRESSURE
    const double time = params.get("total time",-1.0);
    double orthopressureload = 0.0;

    if(time > 1.0)
      orthopressureload = BEAM3EBANISOTROPORTHOPRESSURE * (time-1.0)/0.1;
    if(time > 1.1)
      orthopressureload = BEAM3EBANISOTROPORTHOPRESSURE;
  #endif

  //number of nodes fixed for these element
  const int nnode = 2;

  //number of values per node
  const int vnode = 2;  //value + first derivative

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,1> f_int(true);

  LINALG::Matrix<6*nnode+twistdofs,6*nnode+twistdofs> R_orthopressure_tot(true);
  LINALG::Matrix<6*nnode+twistdofs,1> Res_orthopressure_tot(true);

  //matrix for current nodal positions and nodal tangents
  std::vector<FAD> disp_totlag(6*nnode+twistdofs, 0.0);

  //quantities to store the total internal energy of the element
  FAD temp_energy=0.0;
  int_energy_=0.0;

  //values needed to calculate f_int and K_T

  FAD abs_r_s=0.0;  //|r'|
  FAD kappa_g3=0.0; // scalar product vec(kappa)*vec(g3)
  FAD kappa_g2=0.0; // scalar product vec(kappa)*vec(g2)
  FAD rxTrxx=0.0;    //scalar product: r'T r'' = r''T r'
  FAD epsilon=0.0; // axial tension epsilon=|r_s|-1

  //the reference tangent equals the material tangent: t=g1
  LINALG::TMatrix<FAD,3,1> tangent;
  LINALG::TMatrix<FAD,3,1> tangent_s; //t'=g1'

  //matrices holding the shape functions
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N;
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_s;
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_ss;
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_sss;
  LINALG::TMatrix<FAD,1,6*nnode+twistdofs> C;
  LINALG::TMatrix<FAD,1,6*nnode+twistdofs> C_s;

  //Matrices for N_i, N_i,xi and N_i,xixi. vnode*nnode due to hermite shapefunctions
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xi;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xixi;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xixixi;
  LINALG::TMatrix<FAD,1,twistdofs> C_i;
  LINALG::TMatrix<FAD,1,twistdofs> C_i_xi;

  //Matrices for N_i,s and N_i,ss
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_s;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_ss;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_sss;
  LINALG::TMatrix<FAD,1,twistdofs> C_i_s;

  //Matrices to store r',r''
  LINALG::TMatrix<FAD,3,1> r_s;
  LINALG::TMatrix<FAD,3,1> r_ss;
  LINALG::TMatrix<FAD,3,1> r_sss;

  //tempvector temp_vector needed to store 1/|r_s|^2 r_ss -2(r_s^T*r_ss)/|r_s|^4 r_s
  LINALG::TMatrix<FAD,3,1> temp_vector;

  //scalar values gamma and gamma'
  FAD gamma_;
  FAD gamma_s;

  //spinmatrices
  LINALG::TMatrix<FAD,3,3> Srx;
  LINALG::TMatrix<FAD,3,3> Srxx;
  LINALG::TMatrix<FAD,3,3> Stemp_vector;

  LINALG::TMatrix<FAD,6*nnode+twistdofs,3> M1;
  LINALG::TMatrix<FAD,6*nnode+twistdofs,3> M2;
  LINALG::TMatrix<FAD,6*nnode+twistdofs,3> M2_aux;

  //spatial force stress resultant n and moment stress resultant m
  LINALG::TMatrix<FAD,3,1> m;
  LINALG::TMatrix<FAD,3,1> n;


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

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebanisotrop);

  //Get DiscretizationType of beam element
  const DRT::Element::DiscretizationType distype = Shape();

  //! Calculates nodal positions, tangents and relative twist angles out of the corresponding displacements
  //if TWISTDOFS = 2: [ r1 t1 gamma1 r2 t2 gamma2]
  //if TWISTDOFS = 3: [ r1 t1 gamma1 r2 t2 gamma2 gamma3]
  //if TWISTDOFS = 4: [ r1 t1 gamma1 r2 t2 gamma2 gamma3 gamma4]
  UpdateDispTotlag(disp, disp_totlag);

  //begin: quantities which are needed for NSRISR calculation

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

  //Calculate the nodal reference triads as well as the difference angle between the reference triad triad_ref_nodes[1]
  // and the intermediate triad (normal_bar, binormal_bar) at the right element node

  FAD taubar0;
  taubar0=0.0;

  DetermineNodalTriads(disp_totlag, triads_mat_nodes, theta_nodes,true, triads_ref_nodes,taubar0);

  //end: quantities which are needed for NSRISR calculation

  //Calculation of tension at collocation points; this is needed for ANS approach
  LINALG::TMatrix<FAD,3,1> epsilonCP(true);
  LINALG::TMatrix<FAD,6*nnode+twistdofs,3> LINepsilonCP_T(true);
  LINALG::TMatrix<FAD,6*nnode+twistdofs,1> LINepsilonCP_T_GP(true);
  //Calculate first epsilon in the element midpoint (2nd CP)
  LINALG::TMatrix<FAD,3,1> r_s_midpoint;
  double xi=0.0;
  Calculate_r_s(disp_totlag, xi, r_s_midpoint);
  //Fill all three epsilon values at CPs in matrix epsilonCP
  CalculateEpsilonCP(disp_totlag, r_s_midpoint, epsilonCP, LINepsilonCP_T);

  //Calculate ANS interpolation of kappa
  #ifdef ANSKAPPA
    if(ANSKAPPA<2)
      dserror("ANSKAPPA<2 !!!");

    //Calculation of tension at collocation points; this is needed for ANS approach
    std::vector<LINALG::TMatrix<FAD,3,1> > kappaCP(ANSKAPPA, LINALG::TMatrix<FAD,3,1>(true));

    for(int i=0;i<ANSKAPPA;i++)
    {
      double xi= 0.0;

      if(i==0)
        xi=-1.0;
      else if(i==1)
        xi=1.0;
      else
        xi=-1.0+ (i-1)*2.0/(ANSKAPPA-1);

      //Calculate first epsilon in the element midpoint (2nd CP)
      LINALG::TMatrix<FAD,3,1> r_s_CP(true);
      LINALG::TMatrix<FAD,3,1> r_ss_CP(true);
      //std::cout << "xi: " << xi << std::endl;
      Calculate_r_s_and_r_ss(disp_totlag, xi, r_s_CP, r_ss_CP);
      kappaCP[i]=calculate_curvature(r_s_CP, r_ss_CP);
      //std::cout << "kappaCP[i](2): " << kappaCP[i](2) << std::endl;
    }
  #endif

//  if (this->Id()==0)
//  {
//    FILE *outFile_torsion;
//    outFile_torsion = fopen("beamtorsionalmoments.txt", "w");
//    fclose(outFile_torsion);
//
//    FILE *outFile_bending;
//    outFile_bending = fopen("beambendingmoments.txt", "w");
//    fclose(outFile_bending);
//
//    FILE *outFile_normalforce;
//    outFile_normalforce = fopen("beamnormalforce.txt", "w");
//    fclose(outFile_normalforce);
//  }

//**********************begin: gauss integration*********************************************************************

  double bending_energy=0.0;
  double tension_energy=0.0;

  //Loop through all GP and calculate their contribution to the internal forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {

  //Initialization and calculation of various values

    abs_r_s=0.0;  //|r'|
    kappa_g3=0.0; // scalar product vec(kappa)*vec(g3)
    kappa_g2=0.0; // scalar product vec(kappa)*vec(g2)
    rxTrxx=0.0;
    epsilon=0.0; // axial tension epsilon=|r_s|-1
    //gamma and gamma'
    gamma_=0.0;
    gamma_s=0.0;

    //initialization of matrices
    //shape function N_i
    N_i.Clear();
    //shape function derivatives in xi
    N_i_xi.Clear();
    N_i_xixi.Clear();
    N_i_xixixi.Clear();
    C_i.Clear();
    C_i_xi.Clear();

    //shape function derivatives in s
    N_i_s.Clear();
    N_i_ss.Clear();
    N_i_sss.Clear();
    C_i_s.Clear();

    //Assembled matrices of shape functions
    N.Clear();
    N_s.Clear();
    N_ss.Clear();
    N_sss.Clear();
    C.Clear();
    C_s.Clear();
    //r' and r''
    r_s.Clear();
    r_ss.Clear();
    r_sss.Clear();
    temp_vector.Clear();

    //spinmatrices
    Srx.Clear();
    Srxx.Clear();
    Stemp_vector.Clear();

    tangent.Clear();
    tangent_s.Clear();

    M1.Clear();
    M2.Clear();
    M2_aux.Clear();
    m.Clear();
    n.Clear();
    tau_bar_gp=0.0;
    LINepsilonCP_T_GP.Clear();

    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //specific-for------------Smallest Rotation and Cross Product
    //Get hermite derivatives N'xi and N''xi
    DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv3(N_i_xixixi,xi,length_,line2);
    DRT::UTILS::shape_function_1D(C_i,xi,distype);
    DRT::UTILS::shape_function_1D_deriv1(C_i_xi,xi,distype);

    AssembleShapefunctions(N_i, N_i_xi, N_i_xixi, N_i_xixixi, C_i, C_i_xi, jacobi_[numgp], jacobi2_[numgp], jacobi3_[numgp], N_s, N_ss, N_sss, C, C_s);
    AssembleShapefunctions(N_i,N);

    //Calculation of r' and r'', gamma and gamma' at the gp
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+twistdofs; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag[j];
        r_ss(i)+=N_ss(i,j)*disp_totlag[j];
        r_sss(i)+=N_sss(i,j)*disp_totlag[j];
      }
    }
    for (int i=0; i<6*nnode+twistdofs; i++)
    {
      gamma_+=C(i)*disp_totlag[i];
      gamma_s+=C_s(i)*disp_totlag[i];
    }

  #ifdef BEAM3EBANISOTROPORTHOPRESSURE
    LINALG::TMatrix<FAD,3,1> ortho_normal(true);
    LINALG::Matrix<6*nnode+twistdofs,6*nnode+twistdofs> R_orthopressure(true);
    LINALG::TMatrix<FAD,6*nnode+twistdofs,1> Res_orthopressure(true);

    ortho_normal(0)=r_s(1,0);
    ortho_normal(1)=-r_s(0,0);
    ortho_normal(2)=0.0;
    if (FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(ortho_normal))>1.0e-12)
      ortho_normal.Scale(1.0/(FADUTILS::VectorNorm<3>(ortho_normal)));

    Res_orthopressure.Clear();
    R_orthopressure.Clear();
    Res_orthopressure.MultiplyTN(N,ortho_normal);
    Res_orthopressure.Scale(-orthopressureload* wgt *jacobi_[numgp]);
    for (int i= 0; i<6*nnode+twistdofs; i++)
    {
      for (int j= 0; j<6*nnode+twistdofs; j++)
      {
        R_orthopressure(i,j)=Res_orthopressure(i).dx(j);
      }
    }

    for (int i= 0; i<6*nnode+twistdofs; i++)
    {
      for (int j= 0; j<6*nnode+twistdofs; j++)
      {
        R_orthopressure_tot(i,j)+=R_orthopressure(i,j);
      }
      Res_orthopressure_tot(i)+=Res_orthopressure(i).val();
    }
  #endif

    abs_r_s=Norm(r_s);

#ifndef ANS_BEAM3EB_ANISOTROP
  epsilon=abs_r_s - 1.0;
#else
    epsilon=epsilonCP(0)*(-0.5*xi*(1.0 - xi)) + epsilonCP(1)*(1 - xi*xi) + epsilonCP(2)*(0.5*xi*(1.0 + xi));
    //epsilon=epsilonCP(0)*(0.5*(1.0 - xi))  + epsilonCP(2)*(0.5*(1.0 + xi));
    for (int i=0;i<12+TWISTDOFS;i++)
    {
      LINepsilonCP_T_GP(i)=LINepsilonCP_T(i,0)*(-0.5*xi*(1.0 - xi))+LINepsilonCP_T(i,1)*(1 - xi*xi)+LINepsilonCP_T(i,2)*(0.5*xi*(1.0 + xi));
    }
#endif

    tangent=r_s;
    tangent.Scale(1.0/abs_r_s);
    rxTrxx=ScalarProduct(r_s,r_ss);

    for (int i=0; i<3; i++)
    {
      temp_vector(i)=r_ss(i)/pow(abs_r_s,2.0)-2*r_s(i)*rxTrxx/pow(abs_r_s,4.0);
    }

//      for (int i=0; i<3; i++)
//      {
//        temp_vector(i)=r_ss(i)/pow(epsilon+1,2.0)-2*r_s(i)*rxTrxx/pow(epsilon+1,4.0);
//      }

    LARGEROTATIONS::computespin<FAD>(Stemp_vector,temp_vector);
    for (int i=0;i<3;i++)
    {
      tangent_s(i)=r_ss(i)/abs_r_s-r_s(i)*rxTrxx/pow(abs_r_s,3.0);
    }

    //calculate Frenet-Serret triad: triad[0]=tangent triad[1]=normal_fs triad[2]=binormal_fs
    //attention: the vectors normal and binormal are set to zero for almost straight (|r'xr''|<1.0e-12) beams
    std::vector<LINALG::TMatrix<FAD,3,1> > triad_fs(3);
    triad_fs=calculate_fs_triad(r_s, r_ss);

    //calculate initial curvature
    LINALG::TMatrix<FAD,3,1> kappa_vec(true);

    #ifndef ANSKAPPA
      kappa_vec=calculate_curvature(r_s, r_ss);
    #else
      switch(ANSKAPPA)
      {
        case 2:
        {
          LINALG::Matrix<2,1> L(true);
          DRT::UTILS::shape_function_1D(L,xi,line2);

          for(int i=0;i<ANSKAPPA;i++)
          {
            for(int j=0;j<3;j++)
            {
              kappa_vec(j)+=kappaCP[i](j)*L(i);
            }
          }
        }
        break;

        case 3:
        {
          LINALG::Matrix<3,1> L(true);
          DRT::UTILS::shape_function_1D(L,xi,line3);

          for(int i=0;i<ANSKAPPA;i++)
          {
            for(int j=0;j<3;j++)
            {
              kappa_vec(j)+=kappaCP[i](j)*L(i);
            }
          }
        }
        break;

        case 4:
        {
          LINALG::Matrix<4,1> L(true);
          DRT::UTILS::shape_function_1D(L,xi,line4);

          for(int i=0;i<ANSKAPPA;i++)
          {
            for(int j=0;j<3;j++)
            {
              kappa_vec(j)+=kappaCP[i](j)*L(i);
            }
          }
        }
        break;

        case 5:
        {
          LINALG::Matrix<5,1> L(true);
          DRT::UTILS::shape_function_1D(L,xi,line5);

          for(int i=0;i<ANSKAPPA;i++)
          {
            for(int j=0;j<3;j++)
            {
              kappa_vec(j)+=kappaCP[i](j)*L(i);
            }
          }
        }
        break;

        default:
          dserror("Only the case ANSKAPPA=2, ANSKAPPA=3, ANSKAPPA=4 and ANSKAPPA=5 are implemented so far!");
        break;
      }
    #endif

    //Calculate the intermediate triad triad_bar_gp and the corresponding torsion tau_bar_gp at the gauss points
    CalculateIntermediateTriad(r_s, r_ss, numgp, triads_ref_nodes, triad_bar_gp, tau_bar_gp);

    //calculate the material triad triad_mat_gp at current gp and the mechanical torsion tau_gp out of the reference system
    //attention: for the initial configuration the angle gamma_0 is set to zero!
    CalculateMaterialTriad(xi, numgp, theta_nodes,triad_bar_gp, tau_bar_gp, gamma_, triad_mat_gp,tau_gp);

    //Needed spinmatrices
    LARGEROTATIONS::computespin<FAD>(Srx,r_s);
    LARGEROTATIONS::computespin<FAD>(Srxx,r_ss);

    kappa_g2=0.0;
    kappa_g3=0.0;
    for (int i=0;i<3;i++)
    {
      kappa_g2+=kappa_vec(i)*triad_mat_gp(0,i);
      kappa_g3+=kappa_vec(i)*triad_mat_gp(1,i);
    }

    //Calculation of stress resultants m and n
    for (int i=0;i<3;i++)
    {
      m(i,0)= sm*Irr_*(tau_gp + gamma_s - tau0_[numgp])*tangent(i) + ym*Iyy_*(kappa_g2 - kappa0g20_[numgp])*triad_mat_gp(0,i) + ym*Izz_*(kappa_g3-kappa0g30_[numgp])*triad_mat_gp(1,i);
      n(i,0)= ym*crosssec_*epsilon*tangent(i);
    }

    //for (int i=0;i<3;i++)
    //{
      //m(i)=ym*Izz_/pow(abs_r_s,2)*kappa_vec(i)+sm*Irr_*(tau_gp + gamma_s)*tangent(i);
    //}

//    double Mt = sm*Irr_*(tau_gp + gamma_s - tau0_[numgp]).val();
//    double Mb = pow(pow(ym*Iyy_*(kappa_g2 - kappa0g20_[numgp]).val(),2.0)+pow(ym*Izz_*(kappa_g3-kappa0g30_[numgp]).val(),2.0),0.5);
//    double N = (ym*crosssec_*epsilon).val();

//    FILE *outFile_torsion;
//    outFile_torsion = fopen("beamtorsionalmoments.txt", "a");
//    fprintf(outFile_torsion, "%.16e  ", Mt);
//    fprintf(outFile_torsion, "\n");
//    fclose(outFile_torsion);
//
//    FILE *outFile_bending;
//    outFile_bending = fopen("beambendingmoments.txt", "a");
//    fprintf(outFile_bending, "%.16e  ", Mb);
//    fprintf(outFile_bending, "\n");
//    fclose(outFile_bending);
//
//    FILE *outFile_normalforce;
//    outFile_normalforce = fopen("beamnormalforce.txt", "a");
//    fprintf(outFile_normalforce, "%.16e  ", N);
//    fprintf(outFile_normalforce, "\n");
//    fclose(outFile_normalforce);

    M2_aux.MultiplyTN(N_ss,Srx);
    M2.MultiplyTN(N_s,Stemp_vector);
    M2.Update(1.0/pow(abs_r_s,2.0),M2_aux,1.0);

    //Calculation of auxiliary matrices M1 and M2
    for (int row=0;row<6*nnode+twistdofs;row++)
    {
      for (int column=0;column<3;column++)
      {
        M1(row,column)+=C_s(0,row)*tangent(column, 0) + C(0,row)*tangent_s(column, 0);
      }
    }

#ifdef CONSISTENTANS
    for (int row=0; row<6*nnode+twistdofs; row++)
    {
      for (int column=0; column<3; column++)
      {
        f_int(row)+=wgt*jacobi_[numgp]*((M2(row,column)-M1(row,column))*m(column));
      }
      f_int(row)+=wgt*jacobi_[numgp]*(-ym*crosssec_*epsilon*LINepsilonCP_T_GP(row));
    }
#else
    //Calculation of the internal force vector (valid for all! triads)
    for (int row=0; row<6*nnode+twistdofs; row++)
    {
      for (int column=0; column<3; column++)
      {
        f_int(row)+=wgt*jacobi_[numgp]*(-N_s(column,row)*n(column)+(M2(row,column)-M1(row,column))*m(column));
      }
    }
#endif

    //Calculate internal energy
    temp_energy=0.0;
    temp_energy=ym*crosssec_*pow(epsilon,2.0);
    temp_energy+=sm*Irr_*pow(tau_gp + gamma_s - tau0_[numgp],2.0);
    temp_energy+=ym*Izz_*(pow(kappa_g2 - kappa0g20_[numgp],2.0)+pow(kappa_g3 - kappa0g30_[numgp],2.0));
    temp_energy=temp_energy*0.5*wgt*jacobi_[numgp];
    int_energy_+=temp_energy.val();
    //bending_energy+=(0.5*wgt*jacobi_[numgp]*(ym*Izz_*(pow(kappa_g2 - kappa0g20_[numgp],2.0)+pow(kappa_g3 - kappa0g30_[numgp],2.0)))).val();
    bending_energy+=0.5*wgt*jacobi_[numgp]*(pow(m(2),2.0)/(ym*Izz_)).val();
    tension_energy+=(0.5*wgt*jacobi_[numgp]*ym*crosssec_*pow(epsilon,2.0)).val();
    //End: Calculate internal energy

    //****************Update of the old reference triads at gauss points***********************************************
    if (update)
    {
      #if defined(SR1) or defined (VP)
        UpdateGPRefTriads(tangent, triad_mat_gp, triad_mat_gp, kappa_vec, tau_gp, gamma_s, numgp);
      #endif
    }
    //****************end: Update of the old reference triads at gauss points******************************************

  }//end of gauss integration

  //cout << std::setprecision(16) << "bending_energy: " << bending_energy << endl;
  //cout << "tension_energy: " << tension_energy << endl;

#ifdef OPTCURVE
  //**********************begin: optimal curve integration******************************************************************
      //Loop through the tension GP and calculate their contribution to the internal forcevector and stiffnessmatrix
      f_int.Clear();
      n.Clear();
      LINALG::TMatrix<FAD,3,1> r;
      LINALG::TMatrix<FAD,3,1> r_analyt;
      LINALG::TMatrix<FAD,3,1> rx_analyt;
      r_analyt.Clear();
      rx_analyt.Clear();

      FAD phi = 0.0;
      FAD R0= 100;
      double maxangle = PI/4;
      FAD ele_length=R0*maxangle/(NUMELE);
      int ele=this->Id()+1;
      FAD phi1=(ele-1)*ele_length/(R0);
      FAD phi2=(ele)*ele_length/(R0);

      for(int numgp=0; numgp < gausspoints.nquad; numgp++)
      {
        N_i.Clear();
        n.Clear();
        r.Clear();
        r_analyt.Clear();

        //Get location and weight of GP in parameter space
        const double xi = gausspoints.qxg[numgp][0];
        const double wgt = gausspoints.qwgt[numgp];

        //Get hermite derivatives N'xi and N''xi
        DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

        for (int i=0; i<3; i++)
        {
            r(i)+=N_i(0)*disp_totlag[i] + N_i(1)*disp_totlag[3 + i] + N_i(2)*disp_totlag[7 + i] + N_i(3)*disp_totlag[10 + i];
        }

        phi=(1-xi)/2*phi1 + (1+xi)/2*phi2;
        r_analyt(0)=R0*sin(phi);
        r_analyt(1)=R0*(1-cos(phi));
        r_analyt(2)=0.0;

        for (int i=0;i<3;i++)
        {
          n(i)=-r_analyt(i)+r(i);
        }

        for (int i=0;i<3;i++)
        {
          f_int(i)+= wgt*jacobi_[numgp]*N_i(0)*n(i);
        }

        for (int i=0;i<3;i++)
        {
          f_int(3+i)+= wgt*jacobi_[numgp]*N_i(1)*n(i);
        }

        for (int i=0;i<3;i++)
        {
          f_int(7+i)+= wgt*jacobi_[numgp]*N_i(2)*n(i);
        }

        for (int i=0;i<3;i++)
        {
          f_int(10+i)+= wgt*jacobi_[numgp]*N_i(3)*n(i);
        }
      }//end of gauss integration

      #ifdef OPTCURVEAX
      LINALG::TMatrix<FAD,3,1> fax1(true);
      FAD norm1=0.0;
      LINALG::TMatrix<FAD,3,1> fax2(true);
      FAD norm2=0.0;

      for (int i=0; i<3; i++)
      {
        fax1(i)=disp_totlag[3 + i];
        fax2(i)=disp_totlag[10 + i];
      }

      for (int i=0; i<3; i++)
      {
        fax1(i)=disp_totlag[3 + i];
        fax2(i)=disp_totlag[10 + i];
      }
      norm1=Norm(fax1);
      norm2=Norm(fax2);
      fax1.Scale((norm1-1)/norm1);
      fax2.Scale((norm2-1)/norm2);
      double axfac=1.0e12;

      for (int i=0; i<3; i++)
      {
        f_int(3 + i)+=axfac*fax1(i);
        f_int(10 + i)+=axfac*fax2(i);
      }
      #endif

  //**********************end: optimal curve integration******************************************************************
#endif

  //Update stiffness matrix and force vector
  if(stiffmatrix!=NULL)
  {
    //Calculating stiffness matrix with FAD
    for(int i = 0; i < 6*nnode+twistdofs; i++)
    {
      for(int j = 0; j < 6*nnode+twistdofs; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }
    #ifdef BEAM3EBANISOTROPORTHOPRESSURE
      for(int i = 0; i < 6*nnode+twistdofs; i++)
      {
        for(int j = 0; j < 6*nnode+twistdofs; j++)
        {
          (*stiffmatrix)(i,j)+=R_orthopressure_tot(i,j);
        }
      }
    #endif
    #ifdef OPTCURVE
          (*stiffmatrix)(6,6)=1.0;
          (*stiffmatrix)(13,13)=1.0;
          (*stiffmatrix)(14,14)=1.0;
    #endif
  }

  if(force!=NULL)
  {
    for (int i=0; i< 6*nnode+twistdofs; i++)
    {
      (*force)(i)=f_int(i).val();
    }
    #ifdef BEAM3EBANISOTROPORTHOPRESSURE
      for (int i=0; i< 6*nnode+twistdofs; i++)
      {
        (*force)(i)+=Res_orthopressure_tot(i);
      }
    #endif
  }

  //****************Update of the old reference triads at the nodes***********************************************
  if (update)
  {
    #ifdef NSRISR
      UpdateNodalRefTriads(tangent, triads_mat_nodes, triads_ref_nodes, theta_nodes);
    #endif

  }
  //****************end: Update of the old reference triads at the nodes********************************************

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | stiffness matrix            (private)                                                           meier 05/13|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateInertialForces(Teuchos::ParameterList& params,
                                                              std::vector<double>& acc,
                                                              std::vector<double>& vel,
                                                              std::vector<double>& disp,
                                                              Epetra_SerialDenseMatrix* massmatrix,
                                                              Epetra_SerialDenseVector* force_inertia,
                                                              bool update)
{
  const int twistdofs = TWISTDOFS;
  if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  }

  //number of nodes fixed for these element
  const int nnode = 2;

  //number of values per node
  const int vnode = 2;  //value + first derivative

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,1> f_inertia;

  //initialize
  f_inertia.Clear();

    //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,1> f_inertia_m;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,1> f_inertia_f;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> NTN;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> NsTNs;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> NTN_gp;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> NsTNs_gp;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> CTC;

  //internal force vector
  LINALG::TMatrix<FAD,6*nnode+twistdofs,6*nnode+twistdofs> CTC_gp;

  //initialize
  f_inertia.Clear();
  f_inertia_m.Clear();
  f_inertia_f.Clear();
  NTN.Clear();
  NsTNs.Clear();
  CTC.Clear();

  //matrix for current nodal positions and nodal tangents
  std::vector<FAD> disp_totlag(6*nnode+twistdofs, 0.0);
  std::vector<FAD> vel_totlag(6*nnode+twistdofs, 0.0);
  std::vector<FAD> acc_totlag(6*nnode+twistdofs, 0.0);

  //values needed to calculate f_int and K_T

  FAD abs_r_s=0.0;  //|r'|

  //the reference tangent equals the material tangent: t=g1
  LINALG::TMatrix<FAD,3,1> tangent;
  LINALG::TMatrix<FAD,3,1> tangent_t; //dot(t)=dot(g1)
  LINALG::TMatrix<FAD,3,1> tangent_tt_perp; // ddot(g1)-(ddot(g1)*g1)g1

  //angular velocities
  LINALG::TMatrix<FAD,3,1> w; //
  LINALG::TMatrix<FAD,3,1> w_t; //
  LINALG::TMatrix<FAD,3,1> w_perp; //w-(w*g1)g1
  LINALG::TMatrix<FAD,3,1> w_perp_t;
  FAD w_bar_gp;
  FAD w_bar_t_gp;
  std::vector<FAD> w_tilde_nodes;

  LINALG::TMatrix<FAD,3,1> tangent_ref;
  LINALG::TMatrix<FAD,3,1> tangent_ref_t;
  LINALG::TMatrix<FAD,3,1> w_ref;
  LINALG::TMatrix<FAD,2,3> w_ref_nodes;
  LINALG::TMatrix<FAD,2,1> w_bar_nodes;
  FAD theta_t=0.0;
  FAD theta_tt=0.0;
  FAD theta_t_gp=0.0;
  FAD theta_tt_gp=0.0;
  FAD gamma_t=0.0;
  FAD alpha_t=0.0;
  FAD gamma_tt=0.0;
  FAD alpha_tt=0.0;
  std::vector<LINALG::TMatrix<FAD,3,1> > w_ref_perp_t_nodes(2);
  std::vector<FAD> w_ref_parallel_t_nodes(2,0.0);


  //matrices holding the shape functions
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N;
  LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_s;
  LINALG::TMatrix<FAD,1,6*nnode+twistdofs> C;

  //Matrices for N_i, N_i,xi and N_i,xixi. vnode*nnode due to hermite shapefunctions
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i;
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xi;
  LINALG::TMatrix<FAD,1,twistdofs> C_i;

  //Matrices for N_i,s and N_i,ss
  LINALG::TMatrix<FAD,1,vnode*nnode> N_i_s;

  //Matrices to store r',r''
  LINALG::TMatrix<FAD,3,1> r_s;
  LINALG::TMatrix<FAD,3,1> r_s_t;
  LINALG::TMatrix<FAD,3,1> r_s_tt;
  LINALG::TMatrix<FAD,3,1> r_t;
  LINALG::TMatrix<FAD,3,1> r_tt;


  FAD kin_energy_trans=0.0;
  FAD kin_energy_rot=0.0;
  FAD kin_energy_alpha=0.0;

  //scalar values gamma and gamma'
  FAD gamma_;

  //spinmatrices
  LINALG::TMatrix<FAD,3,3> Stangent;

  LINALG::TMatrix<FAD,6*nnode+twistdofs,3> M_aux;

  //spatial force stress resultant n and moment stress resultant m
  LINALG::TMatrix<FAD,3,1> m_inertia;
  LINALG::TMatrix<FAD,3,1> n_inertia;

  //Get integration points for exact integration
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebanisotrop);

  //Get DiscretizationType of beam element
  const DRT::Element::DiscretizationType distype = Shape();

  //! Calculates nodal positions, tangents and relative twist angles out of the corresponding displacements
  //if TWISTDOFS = 2: [ r1 t1 gamma1 r2 t2 gamma2]
  //if TWISTDOFS = 3: [ r1 t1 gamma1 r2 t2 gamma2 gamma3]
  //if TWISTDOFS = 4: [ r1 t1 gamma1 r2 t2 gamma2 gamma3 gamma4]
  UpdateDispVelAccTotlag(disp, vel, acc, disp_totlag, vel_totlag, acc_totlag);

  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      density = actmat->Density();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
    break;
  }

  //begin: quantities which are needed for NSRISR calculation

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

  //Calculate the nodal reference triads as well as the difference angle between the reference triad triad_ref_nodes[1]
  // and the intermediate triad (normal_bar, binormal_bar) at the right element node
  FAD dummy;
  DetermineNodalTriads(disp_totlag, triads_mat_nodes, theta_nodes,true, triads_ref_nodes,dummy);

  DetermineWrefNodes(disp_totlag, vel_totlag, acc_totlag, w_ref_nodes, w_ref_perp_t_nodes, w_ref_parallel_t_nodes, theta_t, theta_tt);

  //end: quantities which are needed for NSRISR calculation

//**********************begin: gauss integration*********************************************************************

  //Loop through all GP and calculate their contribution to the internal forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {

    //Initialization and calculation of various values
    abs_r_s=0.0;
    gamma_=0.0;
    gamma_t=0.0;
    gamma_tt=0.0;
    alpha_t=0.0;
    alpha_tt=0.0;
    N_i.Clear();
    N_i_xi.Clear();
    C_i.Clear();
    N_i_s.Clear();
    N.Clear();
    N_s.Clear();
    C.Clear();
    r_s.Clear();
    r_s_t.Clear();
    r_s_tt.Clear();
    r_t.Clear();
    r_tt.Clear();
    Stangent.Clear();
    tangent.Clear();
    tangent_t.Clear();
    tangent_tt_perp.Clear();
    w.Clear();
    w_t.Clear();
    w_perp.Clear();
    w_perp_t.Clear();
    w_bar_gp=0.0;
    w_bar_t_gp=0.0;
    theta_t_gp=0.0;
    theta_tt_gp=0.0;
    M_aux.Clear();
    m_inertia.Clear();
    n_inertia.Clear();
    tangent_ref.Clear();
    tangent_ref_t.Clear();
    w_ref.Clear();
    NTN_gp.Clear();
    NsTNs_gp.Clear();
    CTC_gp.Clear();

    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //specific-for------------Smallest Rotation and Cross Product
    //Get hermite derivatives N'xi and N''xi
    DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    DRT::UTILS::shape_function_1D(C_i,xi,distype);

    AssembleShapefunctions(N_i, N_i_xi, C_i, jacobi_[numgp], N, N_s, C);

    //Calculation of r' and r'', gamma and gamma' at the gp
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<6*nnode+twistdofs; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag[j];
        r_s_t(i)+=N_s(i,j)*vel_totlag[j];
        r_s_tt(i)+=N_s(i,j)*acc_totlag[j];
        r_tt(i)+=N(i,j)*acc_totlag[j];
        r_t(i)+=N(i,j)*vel_totlag[j];
      }
    }
    for (int i=0; i<6*nnode+twistdofs; i++)
    {
      gamma_+=C(i)*disp_totlag[i];
      gamma_t+=C(i)*vel_totlag[i];
      gamma_tt+=C(i)*acc_totlag[i];
    }

    tangent=r_s;
    abs_r_s=Norm(r_s);
    tangent.Scale(1.0/abs_r_s);

    tangent_t.Update(1.0/abs_r_s,r_s_t,0.0);
    tangent_t.Update(-1.0/pow(abs_r_s,3)*ScalarProduct(r_s,r_s_t),r_s,1.0);

    tangent_tt_perp.Update(1.0/abs_r_s, r_s_tt, 0.0);
    tangent_tt_perp.Update(-2.0/pow(abs_r_s,3)*ScalarProduct(r_s,r_s_t), r_s_t, 1.0);

    w_perp=ScaleVector(VectorProduct(r_s,r_s_t),1.0/pow(abs_r_s,2));

    w_perp_t=ScaleVector(VectorProduct(r_s,r_s_tt),1.0/pow(Norm(r_s),2));
    w_perp_t.Update(-2.0*ScalarProduct(r_s,r_s_t)/pow(Norm(r_s),4),VectorProduct(r_s,r_s_t),1.0);

    for (int i=0;i<3;i++)
    {
      tangent_ref(i)=triads_mat_nodes(i,0);
      w_ref(i)=w_ref_nodes(0,i);
    }

    tangent_ref_t=VectorProduct(w_ref,tangent_ref);

    w_bar_gp=CalculateWISR(tangent, w_perp, tangent_ref, w_ref,false);

    w_bar_t_gp=CalculateWdotISR(tangent, tangent_t, w_perp, w_perp_t, tangent_ref, tangent_ref_t, w_ref, w_ref_perp_t_nodes[0], w_ref_parallel_t_nodes[0],false);

    theta_t_gp=(xi+1)/2.0*theta_t;

    theta_tt_gp=(xi+1)/2.0*theta_tt;

    alpha_t = gamma_t + theta_t_gp + w_bar_gp;

    alpha_tt = gamma_tt + theta_tt_gp + w_bar_t_gp;

    w.Update(1.0, w_perp, 0.0);
    w.Update(alpha_t, tangent , 1.0);

    w_t.Update(1.0, w_perp_t, 0.0);
    w_t.Update(alpha_tt, tangent, 1.0);
    w_t.Update(alpha_t, tangent_t, 1.0);

    //This is a simplified version for the inertia moments, using the assumption Izz_=Irr_
    m_inertia = ScaleVector(w_t,-density*Irr_);

    //the inertia forces
    n_inertia = ScaleVector(r_tt,-density*crosssec_);

    LINALG::TMatrix<FAD,3,1> r_ss_dummy;
    r_ss_dummy.Clear();
    FAD tau_bar_dummy=0.0;
    FAD tau_dummy=0.0;

    kin_energy_trans+=wgt*jacobi_[numgp]*density*crosssec_*ScalarProduct(r_t,r_t);
    kin_energy_rot+=wgt*jacobi_[numgp]*density*Irr_*ScalarProduct(w,w);
    kin_energy_alpha+=wgt*jacobi_[numgp]*density*Irr_*alpha_t*alpha_t;

    //Calculate the intermediate triad triad_bar_gp at the gauss points
    //Since we are not interested in the corresponding torsion here, we give the empty vector r_ss_dummy to the function
    //and receive the quantity tau_dummy, which has no mechanical meaning here
    CalculateIntermediateTriad(r_s, r_ss_dummy, numgp, triads_ref_nodes, triad_bar_gp, tau_bar_dummy);

    //calculate the material triad triad_mat_gp at current gp out of the reference system
    //again, the mechanical torsion is not needed here. Therefore we give only an empty quantity tau_dummy to the function
    CalculateMaterialTriad(xi, numgp, theta_nodes,triad_bar_gp, tau_bar_gp, gamma_, triad_mat_gp, tau_dummy);

    //Needed spinmatrix
    LARGEROTATIONS::computespin<FAD>(Stangent, tangent);

    M_aux.MultiplyTN(N_s,Stangent);
    M_aux.Scale(-1.0/abs_r_s);
    for (int i=0;i<6*nnode+twistdofs;i++)
    {
      for (int j=0;j<3;j++)
      {
        M_aux(i,j)+=C(i)*tangent(j);
      }
    }
    NTN_gp.MultiplyTN(N, N);
    NsTNs_gp.MultiplyTN(N_s, N_s);
    NTN.Update(density*crosssec_*wgt*jacobi_[numgp],NTN_gp,1.0);
    NsTNs.Update(density*Irr_*wgt*jacobi_[numgp],NsTNs_gp,1.0);
    CTC_gp.MultiplyTN(C, C);
    CTC.Update(density*Irr_*wgt*jacobi_[numgp],CTC_gp,1.0);

    //Calculation of the internal force vector (valid for all! triads)
    for (int row=0; row<6*nnode+twistdofs; row++)
    {
      for (int column=0; column<3; column++)
      {
        f_inertia(row)+=wgt*jacobi_[numgp]*(N(column, row)*n_inertia(column)+M_aux(row,column)*m_inertia(column));
        f_inertia_f(row)+=wgt*jacobi_[numgp]*N(column, row)*n_inertia(column);
        f_inertia_m(row)+=wgt*jacobi_[numgp]*M_aux(row,column)*m_inertia(column);
      }
    }

  }//end of gauss integration

  //Update stiffness matrix and force vector
  //if(*massmatrix != Teuchos::null)
  if(massmatrix != NULL)
  {
    double delta_t = params.get<double>("delta time");
    //TODO: get general time integration params here. So far we set a dserror()!!!
    //TODO: in the current implementation all contributions to the stiffness matrix (i.e linearization with respect to
    // d, d_t and d_tt) are multiplied with the factor (1-alpha_m) in strtimint_genalpha.cpp. It chould be checked, if it
    // is not more consistent to multiply the contributions of d and d_t with (1-alpha_f) as it is the case for internal
    // and damping forces!
    //TODO: Compare the update procedure of gamma_t and gamma_tt at the end of the time step (necessary due to the change of
    // the nodal reference triad in every time step) with the procedure suggested in "On the Use of Lie Group Time Integrators in Multibody Dynamics",
    //Olivier Brüls and Alberto Cardona, J. Comput. Nonlinear Dynam. 5(3), 031002 (May 14, 2010) (13 pages)
    dserror("Check, if we have the correct time integration params!!!");
    //old version
    double beta_newmark = params.get<double>("beta_newmark");
    double gamma_newmark = params.get<double>("gamma_newmark");
    //change it to:
    //double beta_newmark = params.get<double>("timintfac_dis");
    //double gamma_newmark = params.get<double>("timintfac_vel");

    //Calculating mass matrix with FAD
    for(int i = 0; i < 6*nnode+twistdofs; i++)
    {
      for(int j = 0; j < 6*nnode+twistdofs; j++)
      {
        (*massmatrix)(i,j)=beta_newmark*delta_t*delta_t*f_inertia(i).dx(j) + gamma_newmark*delta_t*f_inertia(i).dx(j+2*6+TWISTDOFS) + f_inertia(i).dx(j+2*(2*6+TWISTDOFS));
      }
    }

//    std::cout << "massmatrix gesamt: " << std::endl;
//    massmatrix->Print(std::cout);
  }

  if(force_inertia!=NULL)
  {
    //Calculating inertial force
    for(int i = 0; i < 6*nnode+twistdofs; i++)
    {
        (*force_inertia)(i)=f_inertia(i).val();
    }
    //force_inertia->Print(std::cout);
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition (public)                                  meier 09/12|
 *-----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3ebanisotrop::EvaluateNeumann(Teuchos::ParameterList& params,
                                               DRT::Discretization& discretization,
                                               DRT::Condition& condition,
                                               std::vector<int>& lm,
                                               Epetra_SerialDenseVector& elevec1,
                                               Epetra_SerialDenseMatrix* elemat1)
{

  const int twistdofs = TWISTDOFS;
  if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
  {
    dserror("Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  }
  //dimensions of freedom per node
  const int nnode=2;
  //const int dofpn = 6+twistdofpn;
  const int dofgamma = 6; //gamma is the seventh dof

  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

  //Apply DBC needs the nodal values as a FAD structure - currently, Apply_Dirichlet calculates internal values via FAD
  //if the linearizations are hard coded, this step is no longer necessary!
  std::vector<FAD> disp_totlag(6*nnode+twistdofs);
  for (int dof=0; dof<6*nnode+twistdofs; dof++)
  {
      int node=0;
      if(dof>=7)
      {
        node=1;
      }
      disp_totlag[dof]=mydisp[dof];

      //we have to add the values of the reference configuration:
      if(dof<3)
      {
        disp_totlag[dof]+=Nodes()[node]->X()[dof];
      }
      else if(dof>=3 && dof<6)
      {
        disp_totlag[dof]+=Tref_[node](dof-3);
      }
      else if(dof>=6 && dof<7)
      {
        //nothing has to be done here
      }
      else if(dof>=7 && dof<10)
      {
        disp_totlag[dof]+=Nodes()[node]->X()[dof- 7];
      }
      else if(dof>=10 && dof<13)
      {
        disp_totlag[dof]+=Tref_[node](dof-10);
      }
      else if(dof>=13)
      {
        //nothing has to be done here
      }
      //Last thing to do here - set variables for differentiation
      //FAD: disp_totlag[dof] is the dof-variable of nnode*dofpn variables to differentiate later
      disp_totlag[dof].diff(dof,6*nnode+twistdofs);
  }

  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
  if (vel==null) dserror("Cannot get state vectors 'velocity'");
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
  */

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

  // funct is related to the 6 "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

  //find out which node is correct
  const std::vector< int > * nodeids = condition.Nodes();


  //if a point neumann condition needs to be linearized
  if(condition.Type() == DRT::Condition::PointNeumannEB)
  {
    //find out local node number --> this is done since the first element of a neumann point condition is used for this function
    //in this case we do not know whether it is the left or the right node.
    int node = -1;

    if((*nodeids)[0] == Nodes()[0]->Id())
      node = 0;
    else if((*nodeids)[0] == Nodes()[1]->Id())
      node = 1;

    if (node == -1)
      dserror("\nNode could not be found on nodemap!\n");

    //matrix for current tangent, moment at node and crossproduct
    LINALG::Matrix<3,1> tangent;
    LINALG::Matrix<3,1> StM; // t x M
    LINALG::Matrix<3,1> moment;
    LINALG::Matrix<3,3> St; // t x ...
    LINALG::Matrix<3,3> SM; // M x
    double tTM=0.0;  // t^T/|t| M = r'^T/|r'| M
    double abs_tangent = 0.0;  //|t|

    //clear all matrices
    tangent.Clear();
    StM.Clear();
    moment.Clear();
    St.Clear();
    SM.Clear();

    //assemble r' and moment M at node
    for (int i=0; i<3; i++)
    {
      //get current tangent at nodes
      tangent(i)=Tref_[node](i)+mydisp[node*7+3+i];
      moment(i)=(*onoff)[i+3]*(*val)[i+3]*curvefac[i+3];
    }

    //calculate |t|=|r'| at the boundary
    abs_tangent=tangent.Norm2();
    //abs_tangent=1.0;


    //computespin = S ( tangent ) using the spinmatrix in namespace largerotations
    LARGEROTATIONS::computespin(St,tangent);

    //compute r' x M = t x M
    StM.Multiply(St,moment);

    //calculate the scalar value t^T M
    for(int i=0; i<3; i++)
    {
      tTM+=tangent(i)/abs_tangent*moment(i);
    }

    //Assemble force vector ----------------------------------------------------------------------------------------
    //Term - N^T F: add forces to the dofs d1 d2 d3 of the external force vector, fext is multiplied by (-1) in BACI
    for(int i = 0; i < 3 ; i++)
    {
      elevec1(node*7+i)+=-(*onoff)[i]*(*val)[i]*curvefac[i];
    }


    //Term - N'^T (r' x M / |r'|^2) add moments to the dofs t1 t2 t3 of the external force vector. fext is multiplied by (-1) in BACI
    for(int i=0; i<3 ; i++)
    {
      elevec1(node*7+3+i)+=StM(i)/pow(abs_tangent,2.0);
    }

    //Term - C^T r'^T/|r'| M add moment to the dof gamma of the external force vector. fext is multiplied by (-1) in BACI
    elevec1(node*7+dofgamma)+=-tTM;

    //Assemble stiffness matrix ------------------------------------------------------------------------------------

    //spinmatrix = S ( m )
    LARGEROTATIONS::computespin(SM,moment);

    //add the linearisation of f_ext to stiffness matrix
    //all parts have been evaluated at the boundaries which helps simplifying the matrices:
    //terms with N'^T N' do only affect the lines and colums of the dofs t1 t2 t3; N'T N' equals the unit matrix I there
    //terms with C'^T N' do only affect the lines of the dofs t1 t2 t3 and the column of the dof gamma
    //In contrast to the Neumann part of the residual force here is NOT a factor of (-1) needed, as elemat1 is directly added to the stiffness matrix
    //without sign change

    if(elemat1!=NULL)
    {
      for(int i=0; i<3 ; i++)
      {
        for(int j=0; j<3 ; j++)
        {
          (*elemat1)(node*7+3+i,node*7+3+j)+=SM(i,j)/pow(abs_tangent,2.0); //assemble M x into the stiffness matrix, columns: t1 t2 t3, lines: t1 t2 t3
          (*elemat1)(node*7+3+i,node*7+3+j)+=2*StM(i)*tangent(j)/pow(abs_tangent,4.0); //assemble t x M t^T into the stiffness matrix, columns: t1 t2 t3, lines: t1 t2 t3
        }
      }

      for(int j=0; j<3; j++)
      {
        (*elemat1)(node*7+dofgamma,node*7+3+j)+=moment(j)/abs_tangent; //assemble M^T into the stiffness matrix, columns: t1 t2 t3, lines: gamma
        (*elemat1)(node*7+dofgamma,node*7+3+j)+=-tTM*tangent(j)/pow(abs_tangent,2.0); //assemble M^T t t^T = tTM t^T into the stiffness matrix, columns: t1 t2 t3, lines: gamma
                                      //pow(...,2) is correct here, since tTM is already multiplied with 1/|t|
      }
    }
  }

  //if a line neumann condition needs to be linearized
    else if(condition.Type() == DRT::Condition::LineNeumann)
    {

      // gaussian points
      DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::mygaussruleebanisotrop);
      LINALG::Matrix<1,4> N_i;

      //integration loops
      for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
      {

        //cout << "numgp: " << numgp + 1 << endl;
        //integration points in parameter space and weights
        const double xi = gausspoints.qxg[numgp][0];
        const double wgt = gausspoints.qwgt[numgp];

        //Clear matrix for shape functions
        N_i.Clear();

        //evaluation of shape funcitons at Gauss points
        //Get hermite derivatives N'xi and N''xi (jacobi_*2.0 is length of the element)
        DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);

        //position vector at the gauss point at reference configuration needed for function evaluation
        std::vector<double> X_ref(3,0.0);
        //calculate coordinates of corresponding Guass point in reference configuration
        for (int node=0;node<2;node++)
        {
          for (int dof=0;dof<3;dof++)
          {
            X_ref[dof]+=Nodes()[node]->X()[dof]*N_i(2*node)+Tref_[node](dof)*N_i(2*node + 1);
          }
        }

        double fac=0.0;
        fac = wgt * jacobi_[numgp];

        // load vector ar
        double ar[6];

        // loop the dofs of a node
        for (int dof=0; dof<6; ++dof)
          ar[dof] = fac * (*onoff)[dof]*(*val)[dof]*curvefac[dof];
        double functionfac = 1.0;
        int functnum = -1;

        //Check if also moment line Neumann conditions are implemented accidentally and throw error
        for (int dof=3; dof<6; ++dof)
        {
          if (functions) functnum = (*functions)[dof];
          else functnum = -1;

          if (functnum>0)
          {
            dserror("Line Neumann conditions for distributed moments are not implemented for beam3eb so far! Only the function flag 1, 2 and 3 can be set!");
          }
        }

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

          for (int node=0; node<4; ++node)
          {
            if (node<2)
            elevec1[node*3 + dof] -= N_i(node) *ar[dof] *functionfac;
            else
            elevec1[node*3 + dof + 1] -= N_i(node) *ar[dof] *functionfac;
          }
        }
      } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)
    }

  return 0;

}  //DRT::ELEMENTS::Beam3ebanisotrop::EvaluateNeumann


/*-----------------------------------------------------------------------------------------------------------*
 |  Evaluate elastic energy of element (public)                                                   meier 09/12|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateEnergy(Epetra_SerialDenseVector& energy)
{
  energy[0]+=int_energy_;

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Update of the nodal triads at the end of time step (public)                                   meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateTriads( Teuchos::ParameterList& params,
                                              std::vector<double>& vel,
                                              std::vector<double>& disp)
{
  CalculateInternalForces(params,vel,disp,NULL,NULL,true);

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate nodal reference triad (public)                                                      meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::DetermineNodalTriads(std::vector<FAD> disp_totlag,
                                                       LINALG::TMatrix<FAD,3,6>& triads_mat_nodes,
                                                       std::vector<FAD>& theta_nodes,
                                                       bool refcalc,
                                                       std::vector<LINALG::TMatrix<FAD,2,3> >& triads_ref_nodes,
                                                       FAD& taubar0)
{

    const int twistdofs = TWISTDOFS;
    if(twistdofs!=2 and twistdofs!=3 and twistdofs!=4)
    {
      dserror("Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
    }

    //number of nodes fixed for these element
    const int nnode = 2;

    //number of values per node (functional value + derivative value)
    const int vnode = 2;

    //jacobi determinants
    FAD jacobi_local=0.0;
    FAD jacobi2_local=0.0;
    FAD jacobi3_local=0.0;
    LINALG::TMatrix<FAD,3,1> r0_xi;
    LINALG::TMatrix<FAD,3,1> r0_xixi;
    LINALG::TMatrix<FAD,3,1> r0_xixixi;

    //Matrix to store r'
    LINALG::TMatrix<FAD,3,1> r_s;
    LINALG::TMatrix<FAD,3,1> r_ss;
    LINALG::TMatrix<FAD,3,1> r_sss;
    //Matrix for |r'|
    FAD abs_r_s=0.0;
    LINALG::TMatrix<FAD,3,1> tangent;
    //intermediate triad at right element node
    LINALG::TMatrix<FAD,3,1> normal_bar;
    LINALG::TMatrix<FAD,3,1> binormal_bar;

    //the reference tangent at the left nodes equals the corresponding material tangent: t=g1
    //This matrix is only cleared once in the beginning but not in the node loop, since we want to store the value of the first node during the loop
    LINALG::TMatrix<FAD,3,1> tangent_leftnode;
    tangent_leftnode.Clear();

    //matrices holding the shape functions
    LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_s;
    LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_ss;
    LINALG::TMatrix<FAD,3,6*nnode+twistdofs> N_sss;

    LINALG::TMatrix<FAD,1,6*nnode+twistdofs> C;
    LINALG::TMatrix<FAD,1,6*nnode+twistdofs> C_s;

    LINALG::TMatrix<FAD,1,vnode*nnode> N_i;
    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xi;
    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xixi;
    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_xixixi;
    LINALG::TMatrix<FAD,1,twistdofs> C_i;
    LINALG::TMatrix<FAD,1,twistdofs> C_i_xi;

    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_s;
    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_ss;
    LINALG::TMatrix<FAD,1,vnode*nnode> N_i_sss;
    LINALG::TMatrix<FAD,1,twistdofs> C_i_s;

    //Get DiscretizationType of beam element
    const DRT::Element::DiscretizationType distype = Shape();

    for(int pos=0; pos < 2; pos++)
    {
      //Get location of node or GP in parameter space
      double xi=-1.0 + pos*2.0;

      r_s.Clear();
      r_ss.Clear();
      r_sss.Clear();
      abs_r_s=0.0;
      normal_bar.Clear();
      binormal_bar.Clear();
      tangent.Clear();

      jacobi_local=0.0;
      jacobi2_local=0.0;
      jacobi3_local=0.0;
      r0_xi.Clear();
      r0_xixi.Clear();
      r0_xixixi.Clear();

      N_i.Clear();
      N_i_xi.Clear();
      N_i_xixi.Clear();
      N_i_xixixi.Clear();
      C_i.Clear();
      C_i_xi.Clear();
      N_i_s.Clear();
      N_i_ss.Clear();
      N_i_sss.Clear();
      C_i_s.Clear();
      N_s.Clear();
      N_ss.Clear();
      N_sss.Clear();
      C.Clear();
      C_s.Clear();

      //Get hermite derivatives N'xi and N''xi
      DRT::UTILS::shape_function_hermite_1D(N_i,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi,xi,length_,line2);
      DRT::UTILS::shape_function_hermite_1D_deriv3(N_i_xixixi,xi,length_,line2);

      DRT::UTILS::shape_function_1D(C_i,xi,distype);

      DRT::UTILS::shape_function_1D_deriv1(C_i_xi,xi,distype);

      //Calculate local jacobi factors. We need this as all terms evaluated here are evaluated at the element nodes and not at the gauss points
      for (int i=0;i<3;i++)
      {
        r0_xi(i)+=Nodes()[0]->X()[i]*N_i_xi(0) + Nodes()[1]->X()[i]*N_i_xi(2) +Tref_[0](i)*N_i_xi(1) + Tref_[1](i)*N_i_xi(3);
        r0_xixi(i)+=Nodes()[0]->X()[i]*N_i_xixi(0) + Nodes()[1]->X()[i]*N_i_xixi(2) +Tref_[0](i)*N_i_xixi(1) + Tref_[1](i)*N_i_xixi(3);
        r0_xixixi(i)+=Nodes()[0]->X()[i]*N_i_xixixi(0) + Nodes()[1]->X()[i]*N_i_xixixi(2) +Tref_[0](i)*N_i_xixixi(1) + Tref_[1](i)*N_i_xixixi(3);
      }
      for (int i=0; i<3; i++)
      {
        jacobi_local+=pow(r0_xi(i),2.0);
      }
      jacobi_local=pow(jacobi_local,0.5);
      for (int i=0; i<3; i++)
      {
        jacobi2_local+=r0_xi(i)*r0_xixi(i);
      }
      for (int i=0; i<3; i++)
      {
        jacobi3_local+=r0_xixi(i)*r0_xixi(i)+r0_xi(i)*r0_xixixi(i);
      }

      AssembleShapefunctions(N_i, N_i_xi, N_i_xixi, N_i_xixixi, C_i, C_i_xi, jacobi_local, jacobi2_local, jacobi3_local, N_s, N_ss, N_sss, C, C_s);

      //Calculation of r' and r'', gamma and gamma' at the gp
      for (int i=0; i<3; i++)
      {
        for (int j=0; j<6*nnode+twistdofs; j++)
        {
          r_s(i)+=N_s(i,j)*disp_totlag[j];
          r_ss(i)+=N_ss(i,j)*disp_totlag[j];
          r_sss(i)+=N_sss(i,j)*disp_totlag[j];
        }
      }
      if (pos==0)
      {
        FAD abs_r_s=Norm(r_s);
        taubar0=-ScalarProduct(r_s,VectorProduct(r_ss,r_sss))/(pow(abs_r_s,3.0)+pow(abs_r_s,5.0))*pow(jacobi_local,2.0);
      }
      abs_r_s=Norm(r_s);
      for (int i=0; i<3; i++)
      {
        tangent(i)=r_s(i)/abs_r_s;
        triads_mat_nodes(i,3*pos)=tangent(i);
      }

      if (pos==0)
      {
        for (int i=0; i<3; i++)
        {
          tangent_leftnode(i)=tangent(i);
        }
      }

      LINALG::TMatrix<FAD,3,3> triad_ref;
      LINALG::TMatrix<FAD,3,3> triad;
      triad.Clear();
      triad_ref.Clear();
      for (int i=0;i<3;i++)
      {
        triad_ref(0,i)=t0_nodes_[pos](i);
        triad_ref(1,i)=n0_nodes_[pos](i);
        triad_ref(2,i)=b0_nodes_[pos](i);
      }

      CalculateSRTriads(r_s,triad_ref,triad);

      for (int i=0;i<3;i++)
      {
        triads_ref_nodes[pos](0,i) = triad(1,i);
        triads_ref_nodes[pos](1,i) = triad(2,i);
      }

      if (pos ==1)
      {
        for (int i=0;i<3;i++)
        {
          triad_ref(0,i)=tangent_leftnode(i);
          triad_ref(1,i)=triads_ref_nodes[0](0,i);
          triad_ref(2,i)=triads_ref_nodes[0](1,i);
        }

        CalculateSRTriads(r_s,triad_ref,triad);

        for (int i=0;i<3;i++)
        {
          normal_bar(i)=triad(1,i);
          binormal_bar(i)=triad(2,i);
        }
      }
    }//end of node loop

    //Calculation of angle theta between triads_ref_nodes and triads_bar_nodes at the right element node
    FAD cos_theta = 0.0;
    FAD sin_theta = 0.0;
    FAD sr_theta_nodes = 0.0;

    for (int i=0;i<3;i++)
    {
      cos_theta+=triads_ref_nodes[1](0,i)*normal_bar(i);
      sin_theta+=triads_ref_nodes[1](0,i)*binormal_bar(i);
    }
    if(cos_theta > 1/pow(2,0.5))
    {
      sr_theta_nodes=asin(sin_theta);
    }
    else if (cos_theta < -1/pow(2,0.5))
    {
      sr_theta_nodes=Signum(sin_theta)*M_PI-asin(sin_theta);
    }
    else
    {
      sr_theta_nodes=Signum(sin_theta)*acos(cos_theta);
    }


    //Shift intervall of theta_nodes
    double theta_min = 0.0;
    bool min = true;
    //Choose the singularity to occur at the farthest distance from the angle of the last time step
    if (theta_nodes_old_[1]<=0)
    {
      theta_min = theta_nodes_old_[1] + M_PI;
      min = false;
    }
    else
    {
      theta_min = theta_nodes_old_[1] - M_PI;
      min = true;
    }

    //so far theta_nodes is within ]-180°;180°]. In the next line the angle is shiftet into the interval [theta_min;theta_min + 360°[
    sr_theta_nodes=ShiftAngleIntervall(sr_theta_nodes, theta_min,min);

    //for the NSRISR method theta_nodes[0] is not needed, but for later developments (e.g. the NSRIFS method) two nodal angles could be needed
    #ifdef NSRISR
      theta_nodes[0]=0.0;
      theta_nodes[1]=sr_theta_nodes;
    #elif defined(SR1) or defined(VP)
      if (isinit_)
      {
        theta_nodes[0]=0.0;
        theta_nodes[1]=0.0;
      }
      else
      {
        theta_nodes[0]=0.0;
        theta_nodes[1]=sr_theta_nodes;
      }
    #else
      dserror("DetermineNodalTriads works so far only for the SR, VP and NSRISR method!");
    #endif

    //extract nodal relative angles
    FAD gamma1 = disp_totlag[6];
    FAD gamma2 = disp_totlag[13];

    //calculate nodal material base vectors
    for (int i=0;i<3;i++)
    {
      triads_mat_nodes(i,1)=triads_ref_nodes[0](0,i)*cos(gamma1)+triads_ref_nodes[0](1,i)*sin(gamma1);
      triads_mat_nodes(i,2)=-triads_ref_nodes[0](0,i)*sin(gamma1)+triads_ref_nodes[0](1,i)*cos(gamma1);
      triads_mat_nodes(i,4)=triads_ref_nodes[1](0,i)*cos(gamma2)+triads_ref_nodes[1](1,i)*sin(gamma2);
      triads_mat_nodes(i,5)=-triads_ref_nodes[1](0,i)*sin(gamma2)+triads_ref_nodes[1](1,i)*cos(gamma2);
    }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate reference angluar velocity (public)                                                 meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::DetermineWrefNodes(std::vector<FAD> disp_totlag,
                                                         std::vector<FAD> vel_totlag,
                                                         std::vector<FAD> acc_totlag,
                                                         LINALG::TMatrix<FAD,2,3>& w_ref_nodes,
                                                         std::vector<LINALG::TMatrix<FAD,3,1> >& w_ref_perp_t_nodes,
                                                         std::vector<FAD>& w_ref_parallel_t_nodes,
                                                         FAD& theta_t,
                                                         FAD& theta_tt)
{
  LINALG::TMatrix<FAD,3,1> r_s;
  LINALG::TMatrix<FAD,3,1> r_s_t;
  LINALG::TMatrix<FAD,3,1> r_s_tt;
  LINALG::TMatrix<FAD,3,1> tangent;
  LINALG::TMatrix<FAD,3,1> tangent_t;
  LINALG::TMatrix<FAD,3,1> tangent1;
  LINALG::TMatrix<FAD,3,1> tangent1_t;
  LINALG::TMatrix<FAD,3,1> tangent2;
  LINALG::TMatrix<FAD,3,1> tangent2_t;
  LINALG::TMatrix<FAD,3,1> w_perp;
  LINALG::TMatrix<FAD,3,1> w_perp_t;
  LINALG::TMatrix<FAD,3,1> w_perp2;
  LINALG::TMatrix<FAD,3,1> w_perp2_t;
  LINALG::TMatrix<FAD,3,1> tangent_ref;
  LINALG::TMatrix<FAD,3,1> tangent_ref_t;
  LINALG::TMatrix<FAD,3,1> w_ref;
  LINALG::TMatrix<FAD,3,1> w_ref_perp_t;
  LINALG::TMatrix<FAD,3,1> w_ref1_perp_t;
  FAD w_ref_parallel_t;
  LINALG::TMatrix<FAD,3,1> w_ref1;
  LINALG::TMatrix<FAD,2,1> w_ref_nodes_parallel;
  FAD w_ref1_parallel_t;

  r_s.Clear();
  r_s_t.Clear();
  r_s_tt.Clear();
  tangent.Clear();
  tangent_t.Clear();
  tangent1.Clear();
  tangent1_t.Clear();
  tangent2.Clear();
  tangent2_t.Clear();
  w_perp.Clear();
  w_perp_t.Clear();
  w_perp2.Clear();
  w_perp2_t.Clear();
  tangent_ref.Clear();
  tangent_ref_t.Clear();
  w_ref.Clear();
  w_ref_perp_t.Clear();
  w_ref1_perp_t.Clear();
  w_ref_parallel_t = 0.0;
  w_ref1.Clear();
  w_ref_nodes_parallel.Clear();
  w_ref1_parallel_t=0.0;

  //loop over the two element nodes
  for (int i=0;i<2;i++)
  {
    for (int j=0;j<3;j++)
    {
      r_s(j)=disp_totlag[3 + 7*i + j];
      r_s_t(j)=vel_totlag[3 + 7*i + j];
      r_s_tt(j)=acc_totlag[3 + 7*i + j];
      tangent_ref(j)=t0_nodes_[i](j);
      w_ref(j)=w0_nodes_[i](j);
      w_ref_perp_t(j)=dw0perpdt_nodes_[i](j);
    }

    w_ref_parallel_t=dw0paralleldt_nodes_[i];

    tangent_ref_t=VectorProduct(w_ref, tangent_ref);
    w_perp=ScaleVector(VectorProduct(r_s,r_s_t),1.0/pow(Norm(r_s),2));

    w_perp_t=ScaleVector(VectorProduct(r_s,r_s_tt),1.0/pow(Norm(r_s),2));
    w_perp_t.Update(-2.0*ScalarProduct(r_s,r_s_t)/pow(Norm(r_s),4),VectorProduct(r_s,r_s_t),1.0);

    tangent=ScaleVector(r_s,1.0/Norm(r_s));

    tangent_t=ScaleVector(r_s_t,1.0/Norm(r_s));
    tangent_t.Update(-ScalarProduct(r_s,r_s_t)/pow(Norm(r_s),3),r_s,1.0);



    if (i==0)
    {
      tangent1.Update(1.0,tangent,0.0);
      tangent1_t.Update(1.0,tangent_t,0.0);
    }
    else
    {
      tangent2.Update(1.0,tangent,0.0);
      w_perp2.Update(1.0,w_perp,0.0);
      tangent2_t.Update(1.0,tangent_t,0.0);
      w_perp2_t.Update(1.0,w_perp_t,0.0);
    }

    bool nodal=true;

    #if (NODALALPHAT==3)
      nodal = false;
    #endif

    w_ref_nodes_parallel(i)=CalculateWISR(tangent, w_perp, tangent_ref, w_ref, nodal);
    alphat_nodes_[i]=(w_ref_nodes_parallel(i)+vel_totlag[6+7*i]).val();

    w_ref_parallel_t_nodes[i]=CalculateWdotISR(tangent, tangent_t, w_perp, w_perp_t, tangent_ref, tangent_ref_t, w_ref, w_ref_perp_t, w_ref_parallel_t, nodal);
    alphatt_nodes_[i]=(w_ref_parallel_t_nodes[i]+acc_totlag[6+7*i]).val();

    for (int j=0;j<3;j++)
    {
      w_ref_nodes(i,j)+=w_ref_nodes_parallel(i)*tangent(j)+w_perp(j);
      w_ref_perp_t_nodes[i](j)+=w_perp_t(j);
    }
  }//end: loop over the two element nodes

  for (int i=0;i<3;i++)
  {
    w_ref1(i)=w_ref_nodes(0,i);
    w_ref1_perp_t(i)=w_ref_perp_t_nodes[0](i);
  }
  w_ref1_parallel_t=w_ref_parallel_t_nodes[0];

  //fill class variables
  for (int i=0;i<2;i++)
  {
    for (int j=0;j<3;j++)
    {
      w_nodes_[i](j)=w_ref_nodes(i,j).val();
      dwperpdt_nodes_[i](j)=w_ref_perp_t_nodes[i](j).val();
      dwparalleldt_nodes_[i]=w_ref_parallel_t_nodes[i].val();
    }
  }

  FAD w2_bar = CalculateWISR(tangent2, w_perp2, tangent1, w_ref1, false);
  theta_t=w_ref_nodes_parallel(1) - w2_bar;
  //cout << "theta_t: " << theta_t << endl;

  FAD w2_bar_t = CalculateWdotISR(tangent2, tangent2_t, w_perp2, w_perp2_t, tangent1, tangent1_t, w_ref1, w_ref1_perp_t, w_ref1_parallel_t, false);
  theta_tt=w_ref_parallel_t_nodes[1] - w2_bar_t;
  //cout << "theta_tt: " << theta_tt << endl;

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate signum function for a FAD value                                                     meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3ebanisotrop::Signum(FAD value)
{
  double sgn = 0.0;

  if (value.val() < 0) sgn = -1.0;
  else         sgn =  1.0;

  return sgn;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Shift definition interval of an given angle theta                                             meier 10/12|
 *-----------------------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::Beam3ebanisotrop::ShiftAngleIntervall(FAD theta, double theta_min, bool min)
{
  FAD theta_new=0.0;

  if (min)
  {
      while (theta.val()<theta_min)
      {
        theta += 2*M_PI;
      }
      theta_new=theta;
  }
  else
  {
      while (theta.val()>theta_min)
      {
        theta -= 2*M_PI;
      }
      theta_new=theta;
  }

  return theta_new;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate material triad out of reference triad                                               meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateMaterialTriad( double xi,
                                                              int numgp,
                                                              std::vector<FAD> theta_nodes,
                                                              LINALG::TMatrix<FAD,2,3> triad_bar_gp,
                                                              FAD tau_bar,
                                                              FAD gamma_gp,
                                                              LINALG::TMatrix<FAD,2,3>& triad_mat_gp,
                                                              FAD& tau_gp)
{

  FAD alpha_diff_gp_sr = (xi+1)/2.0*theta_nodes[1];
  FAD alpha_xi_sr = 1.0/2.0*theta_nodes[1];

  tau_gp = alpha_xi_sr/jacobi_[numgp] + tau_bar;

  for (int i=0; i<3; i++)
  {
    triad_mat_gp(0,i)=triad_bar_gp(0,i)*cos(alpha_diff_gp_sr + gamma_gp)+triad_bar_gp(1,i)*sin(alpha_diff_gp_sr + gamma_gp);
    triad_mat_gp(1,i)=-triad_bar_gp(0,i)*sin(alpha_diff_gp_sr + gamma_gp)+triad_bar_gp(1,i)*cos(alpha_diff_gp_sr + gamma_gp);
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate smallest rotation triad out ot a reference triad and a tangent vector               meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateSRTriads( LINALG::TMatrix<FAD,3,1> r_s,
                        LINALG::TMatrix<FAD,3,3> triad_ref,
                        LINALG::TMatrix<FAD,3,3>& triad)
{

  FAD temp_scalar1=0.0;
  FAD temp_scalar2=0.0;
  FAD temp_scalar3=0.0;
  FAD fac_n0=0.0;
  FAD fac_b0=0.0;
  FAD abs_r_s = Norm(r_s);

  for (int i=0; i<3; i++)
  {
    temp_scalar1+=triad_ref(1,i)*r_s(i);
    temp_scalar2+=triad_ref(2,i)*r_s(i);
    temp_scalar3+=triad_ref(0,i)*r_s(i);
  }

  fac_n0=temp_scalar1/(abs_r_s+temp_scalar3);
  fac_b0=temp_scalar2/(abs_r_s+temp_scalar3);

  if (pow(abs_r_s+temp_scalar3,2.0) < ROTATIONTOL)
  {
    dserror("Rotation to large!!!");
  }

  for (int i=0; i<3; i++)
  {
    triad(0,i)=r_s(i)/abs_r_s;
    triad(1,i)=triad_ref(1,i)-fac_n0*(r_s(i)/abs_r_s+triad_ref(0,i));
    triad(2,i)=triad_ref(2,i)-fac_b0*(r_s(i)/abs_r_s+triad_ref(0,i));
   }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate vector product triad out ot a reference triad and a tangent vector                  meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateVPTriads( LINALG::TMatrix<FAD,3,1> r_s,
                                                         LINALG::TMatrix<FAD,3,3> triad_ref,
                                                         LINALG::TMatrix<FAD,3,3>& triad)
{

  FAD abs_r_s = Norm(r_s);
  LINALG::TMatrix<FAD,3,1> tangent=r_s;
  tangent.Scale(1/abs_r_s);

  LINALG::TMatrix<FAD,3,1> b0;
  b0.Clear();
  for (int i=0; i<3; i++)
  {
    b0(i)=triad_ref(2,i);
  }

  LINALG::TMatrix<FAD,3,1> normal=VectorProduct(b0, tangent);

  FAD abs_b0_x_tangent = Norm(normal);

  if (abs_b0_x_tangent < ROTATIONTOL)
  {
    dserror("Rotation to large!!!");
  }

  normal.Scale(1/abs_b0_x_tangent);

  LINALG::TMatrix<FAD,3,1> binormal=VectorProduct(tangent, normal);
  binormal.Scale(1/Norm(binormal));

  for (int i=0; i<3; i++)
  {
    triad(0,i)=tangent(i);
    triad(1,i)=normal(i);
    triad(2,i)=binormal(i);
   }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate reference anglular velocity in case of ISR interpolation                            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::Beam3ebanisotrop::CalculateWISR(LINALG::TMatrix<FAD,3,1> tangent,
                                                   LINALG::TMatrix<FAD,3,1> w_perp,
                                                   LINALG::TMatrix<FAD,3,1> tangent_ref,
                                                   LINALG::TMatrix<FAD,3,1> w_ref,
                                                   bool nodal)
{
  FAD w_bar=0.0;
  LINALG::TMatrix<FAD,3,1> w_ref_perp;
  w_ref_perp.Clear();

  w_ref_perp.Update(1.0,w_ref,1.0);
  w_ref_perp.Update(-1.0,ScaleVector(tangent_ref,ScalarProduct(w_ref,tangent_ref)),1.0);

#if (NODALALPHAT!=1)
  w_bar = -ScalarProduct(w_perp,tangent_ref)/(1+ScalarProduct(tangent, tangent_ref));
#endif
  if (nodal == false)
  {
    w_bar += ScalarProduct(w_ref,tangent_ref) + ScalarProduct(w_ref_perp,tangent)/(1+ScalarProduct(tangent, tangent_ref));
  }

  return w_bar;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate time derivative reference anglular velocity in case of ISR interpolation            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::Beam3ebanisotrop::CalculateWdotISR(LINALG::TMatrix<FAD,3,1> tangent,
                                                      LINALG::TMatrix<FAD,3,1> tangent_t,
                                                      LINALG::TMatrix<FAD,3,1>  w_perp,
                                                      LINALG::TMatrix<FAD,3,1>  w_perp_t,
                                                      LINALG::TMatrix<FAD,3,1> tangent_ref,
                                                      LINALG::TMatrix<FAD,3,1>  tangent_ref_t,
                                                      LINALG::TMatrix<FAD,3,1> w_ref,
                                                      LINALG::TMatrix<FAD,3,1> w_ref_perp_t,
                                                      FAD w_ref_parallel_t,
                                                      bool nodal)
{

  FAD w_bar_t=0.0;
  LINALG::TMatrix<FAD,3,1> w_ref_perp;
  w_ref_perp.Clear();

  w_ref_perp.Update(1.0,w_ref,1.0);
  w_ref_perp.Update(-1.0,ScaleVector(tangent_ref,ScalarProduct(w_ref,tangent_ref)),1.0);

  if (nodal == false)
  {
    w_bar_t =  w_ref_parallel_t;

    w_bar_t += (ScalarProduct(tangent_t,w_ref_perp)+ScalarProduct(tangent,w_ref_perp_t)-ScalarProduct(tangent_ref_t,w_perp)-ScalarProduct(tangent_ref,w_perp_t))/(1+ScalarProduct(tangent, tangent_ref));

    w_bar_t -= (ScalarProduct(w_ref_perp,tangent)-ScalarProduct(w_perp,tangent_ref))*(ScalarProduct(tangent_t, tangent_ref)+ScalarProduct(tangent, tangent_ref_t))/pow((1+ScalarProduct(tangent, tangent_ref)),2);
  }
  else
  {
#if (NODALALPHAT!=1)
    w_bar_t = -ScalarProduct(tangent_ref,w_perp_t)/(1+ScalarProduct(tangent, tangent_ref)) + ScalarProduct(w_perp,tangent_ref)*(ScalarProduct(tangent_t, tangent_ref))/pow((1+ScalarProduct(tangent, tangent_ref)),2);
#endif
  }

  return w_bar_t;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble all shape functions                                                                   meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i,
                                                              LINALG::TMatrix<FAD,1,4> N_i_xi,
                                                              LINALG::TMatrix<FAD,1,4> N_i_xixi,
                                                              LINALG::TMatrix<FAD,1,4> N_i_xixixi,
                                                              LINALG::TMatrix<FAD,1,TWISTDOFS> C_i,
                                                              LINALG::TMatrix<FAD,1,TWISTDOFS> C_i_xi,
                                                              FAD jacobi_local,
                                                              FAD jacobi2_local,
                                                              FAD jacobi3_local,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_ss,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_sss,
                                                              LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C,
                                                              LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C_s)
{

  //Matrices for N_i,s and N_i,ss
  LINALG::TMatrix<FAD,1,4> N_i_s;
  LINALG::TMatrix<FAD,1,4> N_i_ss;
  LINALG::TMatrix<FAD,1,4> N_i_sss;
  LINALG::TMatrix<FAD,1,TWISTDOFS> C_i_s;

  N_i_s.Clear();
  N_i_ss.Clear();
  N_i_sss.Clear();
  C_i_s.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1/jacobi_local);
  C_i_s=C_i_xi;
  C_i_s.Scale(1/jacobi_local);

  //******end: Reorder shape functions***********************************************************
  for (int i=0; i<4; i++)
  {
    N_i_ss(i)=N_i_xixi(i)/pow(jacobi_local,2)-N_i_xi(i)*jacobi2_local/pow(jacobi_local,4.0);
    N_i_sss(i)=N_i_xixixi(i)/pow(jacobi_local,3.0)-3.0*N_i_xixi(i)*jacobi2_local/pow(jacobi_local,5.0)-N_i_xi(i)*jacobi3_local/pow(jacobi_local,5.0)+4.0*N_i_xi(i)*pow(jacobi2_local,2.0)/pow(jacobi_local,7.0);
  }

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C


  #if defined(TWISTDOFS) && (TWISTDOFS==2)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2};

  #elif defined(TWISTDOFS) && (TWISTDOFS==3)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3};

  #elif defined(TWISTDOFS) && (TWISTDOFS==4)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                         {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                         {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3,4};

  #else
  dserror("TWISTDOFS has to be defined. Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+TWISTDOFS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N_s(j,i)=0;
        N_ss(j,i)=0;
        N_sss(j,i)=0;
      }
      else
      {
        N_s(j,i)=N_i_s(assembly_N[j][i]-1);
        N_ss(j,i)=N_i_ss(assembly_N[j][i]-1);
        N_sss(j,i)=N_i_sss(assembly_N[j][i]-1);
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
 |  Assemble all shape functions                                                                  meier 12/14|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::AssembleShapefunctions( LINALG::Matrix<1,4> N_i_xi,
                                                              LINALG::Matrix<1,4> N_i_xixi,
                                                              double jacobi_local,
                                                              double jacobi2_local,
                                                              LINALG::Matrix<3,2*6+TWISTDOFS>& N_s,
                                                              LINALG::Matrix<3,2*6+TWISTDOFS>& N_ss)
{

  //Matrices for N_i,s and N_i,ss
  LINALG::Matrix<1,4> N_i_s;
  LINALG::Matrix<1,4> N_i_ss;

  N_i_s.Clear();
  N_i_ss.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1/jacobi_local);

  //******end: Reorder shape functions***********************************************************
  for (int i=0; i<4; i++)
  {
    N_i_ss(i)=N_i_xixi(i)/pow(jacobi_local,2)-N_i_xi(i)*jacobi2_local/pow(jacobi_local,4.0);
  }

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C


  #if defined(TWISTDOFS) && (TWISTDOFS==2)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

  #elif defined(TWISTDOFS) && (TWISTDOFS==3)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

  #elif defined(TWISTDOFS) && (TWISTDOFS==4)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                         {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                         {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

  #else
  dserror("TWISTDOFS has to be defined. Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+TWISTDOFS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N_s(j,i)=0;
        N_ss(j,i)=0;
      }
      else
      {
        N_s(j,i)=N_i_s(assembly_N[j][i]-1);
        N_ss(j,i)=N_i_ss(assembly_N[j][i]-1);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble some of the shape functions                                                          meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i,
                                                              LINALG::TMatrix<FAD,1,4> N_i_xi,
                                                              LINALG::TMatrix<FAD,1,TWISTDOFS> C_i,
                                                              FAD jacobi_local,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s,
                                                              LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C)
{

  LINALG::TMatrix<FAD,1,4> N_i_s;

  N_i_s.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1/jacobi_local);

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C

  #if defined(TWISTDOFS) && (TWISTDOFS==2)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2};

  #elif defined(TWISTDOFS) && (TWISTDOFS==3)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                      {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                      {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3};

  #elif defined(TWISTDOFS) && (TWISTDOFS==4)

  int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                         {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                         {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

  int assembly_C[2*6+TWISTDOFS]=      {0,0,0,0,0,0,1,0,0,0,0,0,0,2,3,4};

  #else
  dserror("TWISTDOFS has to be defined. Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
  #endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+TWISTDOFS; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N(j,i)=0;
        N_s(j,i)=0;
      }
      else
      {
        N(j,i)=N_i(assembly_N[j][i]-1);
        N_s(j,i)=N_i_s(assembly_N[j][i]-1);
      }
    }
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
void DRT::ELEMENTS::Beam3ebanisotrop::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i_xi,
                                                              FAD jacobi_local,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s)
{

  LINALG::TMatrix<FAD,1,4> N_i_s;

  N_i_s.Clear();

  //Calculate the derivatives in s
  N_i_s=N_i_xi;
  N_i_s.Scale(1/jacobi_local);

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C

#if defined(TWISTDOFS) && (TWISTDOFS==2)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

#elif defined(TWISTDOFS) && (TWISTDOFS==3)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

#elif defined(TWISTDOFS) && (TWISTDOFS==4)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                       {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                       {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

#else
dserror("TWISTDOFS has to be defined. Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
#endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+TWISTDOFS; i++)
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
 |  Assemble the N_s shape functions                                                              meier 12/14|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i,
                                                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N)
{

  //assembly_N respectively assembly_C is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N respectively C

#if defined(TWISTDOFS) && (TWISTDOFS==2)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0}};

#elif defined(TWISTDOFS) && (TWISTDOFS==3)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                    {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0},
                                    {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0}};

#elif defined(TWISTDOFS) && (TWISTDOFS==4)

int assembly_N[3][2*6+TWISTDOFS]=  {{1,0,0,2,0,0,0,3,0,0,4,0,0,0,0,0},
                                       {0,1,0,0,2,0,0,0,3,0,0,4,0,0,0,0},
                                       {0,0,1,0,0,2,0,0,0,3,0,0,4,0,0,0}};

#else
dserror("TWISTDOFS has to be defined. Only the values 2,3 and 4 are valid for TWISTDOFS!!!");
#endif

  //Assemble the matrices of the shape functions
  for (int i=0; i<2*6+TWISTDOFS; i++)
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
 |  Calculate the intermediate triad and corresponding torsion at gauss points                    meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateIntermediateTriad(LINALG::TMatrix<FAD,3,1> r_s,
                                LINALG::TMatrix<FAD,3,1> r_ss,
                                const int numgp,
                                std::vector<LINALG::TMatrix<FAD,2,3> > triads_ref_nodes,
                                LINALG::TMatrix<FAD,2,3>& triad_bar_gp,
                                FAD& tau_bar_gp)
{

  //Calculate intermediate triad triad_bar_gp
  LINALG::TMatrix<FAD,3,3> triad_ref;
  triad_ref.Clear();

  LINALG::TMatrix<FAD,3,3> triad_bar;
  triad_bar.Clear();

  FAD abs_r_s = Norm(r_s);

  //triad from which the sr or vp mapping is made
  LINALG::TMatrix<FAD,3,1> tangent_ref;
  LINALG::TMatrix<FAD,3,1> normal_ref;
  LINALG::TMatrix<FAD,3,1> binormal_ref;
  for (int i=0;i<3;i++)
  {
    if (!isinit_)
    {
      normal_ref(i)=triads_ref_nodes[0](0,i);
      binormal_ref(i)=triads_ref_nodes[0](1,i);
    }
    else
    {
      #if defined (NSRISR)
        normal_ref(i)=triads_ref_nodes[0](0,i);
        binormal_ref(i)=triads_ref_nodes[0](1,i);
      #elif defined (SR1) or defined (VP)
        normal_ref(i)=n0_[numgp](i);
        binormal_ref(i)=b0_[numgp](i);
      #else
        dserror("DetermineNodalTriads works so far only for the SR, VP and NSRISR method!");
      #endif
    }
  }
  tangent_ref=CompleteTriad(normal_ref, binormal_ref);

  //tangent onto which the sr or vp mapping is made
  LINALG::TMatrix<FAD,3,1> tangent=r_s;

  for (int i=0;i<3;i++)
  {
    tangent(i)=tangent(i)/abs_r_s;
    triad_ref(0,i)=tangent_ref(i);
    triad_ref(1,i)=normal_ref(i);
    triad_ref(2,i)=binormal_ref(i);
  }

  if (!isinit_)
  {
    CalculateSRTriads(r_s, triad_ref, triad_bar);
  }
  else
  {
    #if defined (NSRISR) or defined (SR1)
      CalculateSRTriads(r_s, triad_ref, triad_bar);
    #elif defined (VP)
      CalculateVPTriads(r_s, triad_ref, triad_bar);
    #else
      dserror("DetermineNodalTriads works so far only for the SR, VP and NSRISR method!");
    #endif
  }

  for (int i=0;i<3;i++)
  {
    triad_bar_gp(0,i)=triad_bar(1,i);
    triad_bar_gp(1,i)=triad_bar(2,i);
  }

  //Calculate intermediate torsion tau_bar_gp
  LINALG::TMatrix<FAD,3,1> kappa_vec=calculate_curvature(r_s, r_ss);
  LINALG::TMatrix<FAD,3,1> n0;
  LINALG::TMatrix<FAD,3,1> b0;
  LINALG::TMatrix<FAD,3,1> dt0ds;
  LINALG::TMatrix<FAD,3,1> db0ds;
  for (int i=0;i<3;i++)
  {
    n0(i)=n0_[numgp](i);
    b0(i)=b0_[numgp](i);
    dt0ds(i)=dt0ds_[numgp](i);
    db0ds(i)=db0ds_[numgp](i);
  }

  if (!isinit_)
  {
    tau_bar_gp = -ScalarProduct(tangent_ref, kappa_vec)/(1+ScalarProduct(tangent_ref,tangent));
  }
  else
  {
    #if defined (NSRISR)
      tau_bar_gp = -ScalarProduct(tangent_ref, kappa_vec)/(1+ScalarProduct(tangent_ref,tangent));
    #elif defined (SR1)
      tau_bar_gp=-ScalarProduct(db0ds, n0) -(ScalarProduct(tangent_ref, kappa_vec)-ScalarProduct(db0ds, tangent_ref)*ScalarProduct(n0, tangent)-ScalarProduct(n0, dt0ds)*ScalarProduct(b0, tangent))/(1+ScalarProduct(tangent_ref,tangent));
    #elif defined (VP)
      LINALG::TMatrix<FAD,3,1> tangent_x_b0=VectorProduct(tangent, b0);
      tau_bar_gp=(ScalarProduct(db0ds, tangent_x_b0) + ScalarProduct(b0, tangent)*ScalarProduct(b0, kappa_vec))/pow(Norm(tangent_x_b0),2.0);
    #else
      dserror("DetermineNodalTriads works so far only for the SR, VP and NSRISR method!");
    #endif
  }
  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate L2-Norm for a FAD vector                                                            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::Beam3ebanisotrop::Norm(LINALG::TMatrix<FAD,3,1> v)
{
  FAD norm = 0.0;
  for (int i=0;i<3;i++)
  {
    norm+=v(i)*v(i);
  }
  norm = pow(norm,0.5);

  return norm;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Complete two base vectors to an orthonormal triad                                             meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
LINALG::TMatrix<FAD,3,1> DRT::ELEMENTS::Beam3ebanisotrop::CompleteTriad(LINALG::TMatrix<FAD,3,1> first_basevector, LINALG::TMatrix<FAD,3,1> second_basevector)
{
  LINALG::TMatrix<FAD,3,3> first_basevector_matrix;
  LINALG::TMatrix<FAD,3,1> third_basevector;

  first_basevector_matrix.Clear();
  third_basevector.Clear();

  LARGEROTATIONS::computespin(first_basevector_matrix, first_basevector);
  third_basevector.Multiply(first_basevector_matrix,second_basevector);

  return third_basevector;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Scale FAD vector with a FAD factor                                                            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
LINALG::TMatrix<FAD,3,1> DRT::ELEMENTS::Beam3ebanisotrop::ScaleVector(LINALG::TMatrix<FAD,3,1> vector, FAD factor)
{
  vector.Scale(factor);

  return vector;
}

/*-----------------------------------------------------------------------------------------------------------*
 | Scalar Product of two FAD vectors                                                              meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::Beam3ebanisotrop::ScalarProduct(LINALG::TMatrix<FAD,3,1> first_vector, LINALG::TMatrix<FAD,3,1> second_vector)
{
  FAD result=0.0;

  for (int i=0;i<3;i++)
  {
    result+=first_vector(i)*second_vector(i);
  }

  return result;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Scalar Product of two Linalg Matrizes                                                         meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3ebanisotrop::ScalarProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector)
{
  double result=0.0;

  for (int i=0;i<3;i++)
  {
    result+=first_vector(i)*second_vector(i);
  }

  return result;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Vector product of two FAD vectors                                                             meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
LINALG::TMatrix<FAD,3,1> DRT::ELEMENTS::Beam3ebanisotrop::VectorProduct(LINALG::TMatrix<FAD,3,1> first_vector, LINALG::TMatrix<FAD,3,1> second_vector)
{
  LINALG::TMatrix<FAD,3,1> result_vector;
  result_vector.Clear();
  LINALG::TMatrix<FAD,3,3> S_first_vector;
  S_first_vector.Clear();
  LARGEROTATIONS::computespin(S_first_vector,first_vector);

  result_vector.Multiply(S_first_vector, second_vector);

  return result_vector;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate position vectors from displacement +  initial position                               meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateDispTotlag(std::vector<double> disp, std::vector<FAD>& disp_totlag)
{

  for (int dof=0; dof<2*6+TWISTDOFS; dof++)
    {
        int node=0;
        if(dof>=7)
        {
          node=1;
        }
        disp_totlag[dof]=disp[dof];

        //we have to add the values of the reference configuration:
        if(dof<3)
        {
          disp_totlag[dof]+=Nodes()[node]->X()[dof];
        }
        else if(dof>=3 && dof<6)
        {
          disp_totlag[dof]+=Tref_[node](dof-3);
        }
        else if(dof>=6 && dof<7)
        {
          //nothing has to be done here
        }
        else if(dof>=7 && dof<10)
        {
          disp_totlag[dof]+=Nodes()[node]->X()[dof- 7];
        }
        else if(dof>=10 && dof<13)
        {
          disp_totlag[dof]+=Tref_[node](dof-10);
        }
        else if(dof>=13)
        {
          //nothing has to be done here
        }
        //Last thing to do here - set variables for differentiation
        //FAD: disp_totlag[dof] is the dof-variable of nnode*dofpn variables to differentiate later
        disp_totlag[dof].diff(dof,2*6+TWISTDOFS);
    }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate position, displacement and velocity vectors in FAD style                            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateDispVelAccTotlag( std::vector<double> disp,
                                                              std::vector<double> vel,
                                                              std::vector<double> acc,
                                                              std::vector<FAD>& disp_totlag,
                                                              std::vector<FAD>& vel_totlag,
                                                              std::vector<FAD>& acc_totlag)
{

  for (int dof=0; dof<2*6+TWISTDOFS; dof++)
    {
      int node=0;
      if(dof>=7)
      {
        node=1;
      }
      disp_totlag[dof]=disp[dof];

      //we have to add the values of the reference configuration:
      if(dof<3)
      {
        disp_totlag[dof]+=Nodes()[node]->X()[dof];
      }
      else if(dof>=3 && dof<6)
      {
        disp_totlag[dof]+=Tref_[node](dof-3);
      }
      else if(dof>=6 && dof<7)
      {
        //nothing has to be done here
      }
      else if(dof>=7 && dof<10)
      {
        disp_totlag[dof]+=Nodes()[node]->X()[dof- 7];
      }
      else if(dof>=10 && dof<13)
      {
        disp_totlag[dof]+=Tref_[node](dof-10);
      }
      else if(dof>=13)
      {
        //nothing has to be done here
      }
      vel_totlag[dof]=vel[dof];
      acc_totlag[dof]=acc[dof];

      //Last thing to do here - set variables for differentiation
      //FAD: disp_totlag[dof] is the dof-variable of nnode*dofpn variables to differentiate later
      disp_totlag[dof].diff(dof,3*(2*6+TWISTDOFS));
      vel_totlag[dof].diff(dof+2*6+TWISTDOFS,3*(2*6+TWISTDOFS));
      acc_totlag[dof].diff(dof+2*(2*6+TWISTDOFS),3*(2*6+TWISTDOFS));
    }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate Frenet Serret triad at a curve point                                                meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
std::vector<LINALG::TMatrix<FAD,3,1> > DRT::ELEMENTS::Beam3ebanisotrop::calculate_fs_triad(LINALG::TMatrix<FAD,3,1>& r_s, LINALG::TMatrix<FAD,3,1>& r_ss)
{
  std::vector<LINALG::TMatrix<FAD,3,1> > triads;
  triads.resize(3);

  LINALG::TMatrix<FAD,3,1> tangent;
  LINALG::TMatrix<FAD,3,1> normal;
  LINALG::TMatrix<FAD,3,1> binormal;

  //calculate tangential vector
  for (int i=0; i<3; i++)
  {
    tangent(i)=r_s(i)/Norm(r_s);
  }

  FAD scalar1=0.0;
  FAD scalar2=0.0;

  //spinmatrix Sr' = r'x
  LINALG::TMatrix<FAD,3,3> Srx;
  LARGEROTATIONS::computespin(Srx,r_s);

  //cross-product r'xr''
  LINALG::TMatrix<FAD,3,1> Srxrxx;
  Srxrxx.Clear();
  Srxrxx.Multiply(Srx,r_ss);

  if (Norm(Srxrxx)<1.0e-12)
  {
    for (int i=0; i<3; i++)
    {
      binormal(i)=0.0;
      normal(i)=0.0;
    }
  }
  else
  {
    //calculate normal vector
    for (int i=0; i<3; i++)
    {
      scalar1+=r_s(i)*r_s(i);
      scalar2+=r_s(i)*r_ss(i);
    }
    for (int i=0; i<3; i++)
    {
      normal(i)=(scalar1*r_ss(i)-scalar2*r_s(i))/(Norm(r_s)*Norm(Srxrxx));
    }
    //calculate binormal vector
    for (int i=0; i<3; i++)
    {
      binormal(i)=Srxrxx(i)/Norm(Srxrxx);
    }
  }

  triads[0]=tangent;
  triads[1]=normal;
  triads[2]=binormal;

  return triads;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate tension at gauss point                                                              meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::CalculateEpsilonCP( std::vector<FAD> disp_totlag,
                                                          LINALG::TMatrix<FAD,3,1> r_s_midpoint,
                                                          LINALG::TMatrix<FAD,3,1>& epsilonCP,
                                                          LINALG::TMatrix<FAD,12+TWISTDOFS,3>& LINepsilonCP_T)
{

      LINALG::TMatrix<FAD,2,1> nodal_epsilon;
      nodal_epsilon.Clear();
      FAD midpoint_epsilon=0.0;

      for (int i=0;i<3;i++)
      {
        nodal_epsilon(0)+=pow(disp_totlag[3+i],2.0);
      }
      nodal_epsilon(0)=pow(nodal_epsilon(0),0.5)-1;

      for (int i=0;i<3;i++)
      {
        nodal_epsilon(1)+=pow(disp_totlag[10+i],2.0);
      }
      nodal_epsilon(1)=pow(nodal_epsilon(1),0.5)-1;

      midpoint_epsilon=Norm(r_s_midpoint)-1.0;

      //Here the order of the epsilonCP values is adapted to the nodal node numbering and the corresponding numbering of Larange shape functions for beam elements
      epsilonCP(0)=nodal_epsilon(0);
      epsilonCP(1)=midpoint_epsilon;
      epsilonCP(2)=nodal_epsilon(1);

      for (int i=0;i<3;i++)
      {
        const double xi=-1.0+i*1.0;

        LINALG::TMatrix<FAD,3,1> r0_xi;
        r0_xi.Clear();

        LINALG::TMatrix<FAD,3,1> t;
        t.Clear();

        LINALG::TMatrix<FAD,1,4> N_i_xi;
        N_i_xi.Clear();

        FAD jacobi_local=0.0;

        LINALG::TMatrix<FAD,3,12+TWISTDOFS> N_s;
        N_s.Clear();

        DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

        for (int n=0;n<3;n++)
        {
          r0_xi(n)+=Nodes()[0]->X()[n]*N_i_xi(0) + Nodes()[1]->X()[n]*N_i_xi(2) +Tref_[0](n)*N_i_xi(1) + Tref_[1](n)*N_i_xi(3);
        }

        jacobi_local=Norm(r0_xi);

        AssembleShapefunctions(N_i_xi, jacobi_local, N_s);
        for (int k=0;k<3;k++)
        {
          for (int l=0;l<12+TWISTDOFS;l++)
          {
            t(k)+=N_s(k,l)*disp_totlag[l];
          }
        }
        t=ScaleVector(t,1/Norm(t));
        for (int k=0;k<3;k++)
        {
          for (int l=0;l<12+TWISTDOFS;l++)
          {
            LINepsilonCP_T(l,i)+=N_s(k,l)*t(k);
          }
        }

      }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate r_s at specific point                                                               meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::Calculate_r_s(std::vector<FAD> disp_totlag, double xi, LINALG::TMatrix<FAD,3,1>& r_s_midpoint)
{

          FAD abs_r_s=0.0;
          LINALG::TMatrix<FAD,1,4> N_i_xi;
          LINALG::TMatrix<FAD,3,2*6+TWISTDOFS> N_s;
          N_i_xi.Clear();
          N_s.Clear();
          FAD jacobi_local=0.0;
          LINALG::TMatrix<FAD,3,1> r0_xi;
          r0_xi.Clear();

          //Get hermite derivatives N'xi and N''xi
          DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);

          //Calculate local jacobi factors. We need this as all terms evaluated here are evaluated at the local gauss point and not at the element gauss points
          for (int i=0;i<3;i++)
          {
            r0_xi(i)+=Nodes()[0]->X()[i]*N_i_xi(0) + Nodes()[1]->X()[i]*N_i_xi(2) +Tref_[0](i)*N_i_xi(1) + Tref_[1](i)*N_i_xi(3);
          }
          jacobi_local=Norm(r0_xi);

          AssembleShapefunctions(N_i_xi, jacobi_local, N_s);

          //Calculation of r' and r'', gamma and gamma' at the gp
          for (int i=0; i<3; i++)
          {
            for (int j=0; j<2*6+TWISTDOFS; j++)
            {
              r_s_midpoint(i)+=N_s(i,j)*disp_totlag[j];
            }
          }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Calculate r_s and r_ss at specific point xi                                                   meier 12/14|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::Calculate_r_s_and_r_ss(std::vector<FAD> disp_totlag, double xi, LINALG::TMatrix<FAD,3,1>& r_s, LINALG::TMatrix<FAD,3,1>& r_ss)
{

    LINALG::Matrix<1,4> N_i_xi(true);
    LINALG::Matrix<3,2*6+TWISTDOFS> N_s(true);
    LINALG::Matrix<1,4> N_i_xixi(true);
    LINALG::Matrix<3,2*6+TWISTDOFS> N_ss(true);
    LINALG::Matrix<1,4> N_i_xixixi(true);

    LINALG::Matrix<3,1> r0_xi(true);
    LINALG::Matrix<3,1> r0_xixi(true);
    LINALG::Matrix<3,1> r0_xixixi(true);

    //Get hermite derivatives N'xi and N''xi
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_xi,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_xixi,xi,length_,line2);
    DRT::UTILS::shape_function_hermite_1D_deriv3(N_i_xixixi,xi,length_,line2);


    //Calculate local jacobi factors. We need this as all terms evaluated here are evaluated at the local gauss point and not at the element gauss points
    for (int i=0;i<3;i++)
    {
      r0_xi(i)+=Nodes()[0]->X()[i]*N_i_xi(0) + Nodes()[1]->X()[i]*N_i_xi(2) +Tref_[0](i)*N_i_xi(1) + Tref_[1](i)*N_i_xi(3);
      r0_xixi(i)+=Nodes()[0]->X()[i]*N_i_xixi(0) + Nodes()[1]->X()[i]*N_i_xixi(2) +Tref_[0](i)*N_i_xixi(1) + Tref_[1](i)*N_i_xixi(3);
      r0_xixixi(i)+=Nodes()[0]->X()[i]*N_i_xixixi(0) + Nodes()[1]->X()[i]*N_i_xixixi(2) +Tref_[0](i)*N_i_xixixi(1) + Tref_[1](i)*N_i_xixixi(3);
    }

    //calculate jacobi factors
    std::vector<double> jacobi(3);
    jacobi=calculate_jacobi(r0_xi,r0_xixi,r0_xixixi);

    AssembleShapefunctions(N_i_xi,N_i_xixi, jacobi[0], jacobi[1], N_s, N_ss);

    //Calculation of r' and r'', gamma and gamma' at the gp
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<2*6+TWISTDOFS; j++)
      {
        r_s(i)+=N_s(i,j)*disp_totlag[j];
        r_ss(i)+=N_ss(i,j)*disp_totlag[j];
      }
    }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Update the reference triad at the gauss points at the end of time step                        meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateGPRefTriads( LINALG::TMatrix<FAD,3,1> tangent, LINALG::TMatrix<FAD,2,3> triad_mat_gp,
                                                            LINALG::TMatrix<FAD,2,3> triad_bar_gp, LINALG::TMatrix<FAD,3,1> kappa_vec,
                                                            FAD tau_gp, FAD gamma_s, int numgp)
{
  //*************************************************************************************************
  //Update of the old reference triads: from here From here on, the derivative and the basis
  //vector from the previous time-step are overwritten and therefore no longer available
  //*************************************************************************************************
  FAD kappa_g2 = 0.0;
  FAD kappa_g3 = 0.0;
  //Calculation of the projection of vec(kappa) onto g2 and g3
  for (int i=0; i<3; i++)
  {
    kappa_g3+=kappa_vec(i)*triad_mat_gp(1,i);
    kappa_g2+=kappa_vec(i)*triad_mat_gp(0,i);
  }
  #ifdef MATERIALREF
    for (int i=0; i<3; i++)
    {
      t0_[numgp](i)=tangent(i).val();
      n0_[numgp](i)=triad_mat_gp(0,i).val();
      b0_[numgp](i)=triad_mat_gp(1,i).val();

      dt0ds_[numgp](i)=(kappa_g3*triad_mat_gp(0,i) - kappa_g2*triad_mat_gp(1,i)).val();
      db0ds_[numgp](i)=(kappa_g2*tangent(i) - (tau_gp + gamma_s)*triad_mat_gp(0,i)).val();
    }
  #endif

  #ifndef MATERIALREF
    FAD kappa_normal = 0.0;
    FAD kappa_binormal = 0.0;
    //Calculation of the projection of vec(kappa) onto normal and binormal
    for (int i=0; i<3; i++)
    {
      kappa_binormal+=kappa_vec(i)*triad_bar_gp(1,i);
      kappa_normal+=kappa_vec(i)*triad_bar_gp(0,i);
    }

    for (int i=0; i<3; i++)
    {

      t0_[numgp](i)=tangent(i).val();
      n0_[numgp](i)=triad_bar_gp(0,i).val();
      b0_[numgp](i)=triad_bar_gp(1,i).val();

      dt0ds_[numgp](i)=(kappa_g3*triad_mat_gp(0,i) - kappa_g2*triad_mat_gp(1,i)).val();
      db0ds_[numgp](i)=(kappa_normal*tangent(i) - tau_gp*triad_bar_gp(0,i)).val();


    }
  #endif

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Update the reference triad at the element nodes at the end of time step                       meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateNodalRefTriads( LINALG::TMatrix<FAD,3,1> tangent,
                                                            LINALG::TMatrix<FAD,3,6> triads_mat_nodes,
                                                            std::vector<LINALG::TMatrix<FAD,2,3> > triads_ref_nodes,
                                                            std::vector<FAD> theta_nodes)
{

  theta_nodes_old_[0]=theta_nodes[0].val();
  theta_nodes_old_[1]=theta_nodes[1].val();

  LINALG::TMatrix<FAD,3,1> g2;
  LINALG::TMatrix<FAD,3,1> g3;

  for (int node=0;node<2;node++)
  {
    g2.Clear();
    g3.Clear();
    for (int i=0;i<3;i++)
    {
      tangent(i)=triads_mat_nodes(i,3*node);
      g2(i)=triads_mat_nodes(i,3*node+1);
      g3(i)=triads_mat_nodes(i,3*node+2);
    }
    #ifdef MATERIALREF
      for (int i=0; i<3; i++)
      {
        //Use smoothed Material triad as reference
        t0_nodes_[node](i)=tangent(i).val();
        n0_nodes_[node](i)=g2(i).val();
        b0_nodes_[node](i)=g3(i).val();
      }
    #else
      for (int i=0; i<3; i++)
      {
        //Use smoothed reference triad as reference
        t0_nodes_[node](i)=tangent(i).val();
        n0_nodes_[node](i)=triads_ref_nodes[node](0,i).val();
        b0_nodes_[node](i)=triads_ref_nodes[node](1,i).val();

      }
    #endif
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Update the reference angular velocity at the element nodes at the end of time step            meier 05/13|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ebanisotrop::UpdateNodalAngularVelocities(Epetra_SerialDenseVector* update_disp,
                                                                   Epetra_SerialDenseVector* update_vel,
                                                                   Epetra_SerialDenseVector* update_acc)
{

  //loop over the two element nodes and update of the nodal angular velocities and accelerations
  for (int i=0;i<2;i++)
  {
    for (int j=0;j<3;j++)
    {
      w0_nodes_[i](j)=w_nodes_[i](j);
      dw0perpdt_nodes_[i](j)=dwperpdt_nodes_[i](j);
      dw0paralleldt_nodes_[i]=dwparalleldt_nodes_[i];
    }
  }

#if (NODALALPHAT == 2)
  if (this->Id()==1)
    dserror("Velocity update currently only implemented for one element!!!");

  //update the velocities gamma_t which is necessary due to the update of the reference coordinate system
  if(update_vel!=NULL)
  {
    double w_parallel = 0.0;
    for(int i = 0; i < 2; i++)
    {
        w_parallel = ScalarProduct(w0_nodes_[i], t0_nodes_[i]);
        (*update_vel)(7*i +6)=w_parallel;
    }
  }

  //update the accelerations gamma_tt which is necessary due to the update of the reference coordinate system
  if(update_acc!=NULL)
  {

    for(int i = 0; i < 2; i++)
    {
        (*update_acc)(7*i +6)=dw0paralleldt_nodes_[i];
    }
  }
#endif


  return;
}
