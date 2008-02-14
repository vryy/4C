/*!----------------------------------------------------------------------
\file so_hex8_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_hex8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_hex8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_hex8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_hex8::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_homog_stressdens")     act = So_hex8::calc_homog_stressdens;
  else dserror("Unknown type of action for So_hex8");

  const double time = params.get("total time",-1.0);

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,time);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,time);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
      dserror("Case 'calc_struct_internalforce' not yet implemented");
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,time);
    }
    break;

    // evaluate stresses at nodes
    case calc_struct_stress: {
      if (stresstype_==soh8_stress_ndxyz) {
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        Epetra_SerialDenseMatrix stresses(NUMGPT_SOH8,NUMSTR_SOH8);
        vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stresses,time);

        // extrapolate stresses at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOH8,NUMSTR_SOH8);
        soh8_expol(stresses,nodalstresses);

        // average nodal stresses between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOH8);

        for (int i=0;i<NUMNOD_SOH8;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOH8;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOH8;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep: {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_homog_stressdens: {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8_homog(params, mydisp, time, myres);
    }
    break;

    default:
      dserror("Unknown type of action for Solid3");
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }

  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // compute the Jacobian matrix
    Epetra_SerialDenseMatrix jac(NUMDIM_SOH8,NUMDIM_SOH8);
    jac.Multiply('N','N',1.0,int_hex8.deriv_gp[gp],xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = int_hex8.weights(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
        elevec1[nodid*NUMDIM_SOH8+dim] += int_hex8.shapefct_gp[gp](nodid) * (*onoff)[dim] *\
        								  (*val)[dim] * fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::Shell8::s8_EvaluateNeumann

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* stresses,       // stresses at GP
      const double              time)           // current absolute time
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
  	#if 0 //NODE TRANSLATION AROUND ZERO TO GAIN ACCURACY
    xrefe(i,0) = Nodes()[i]->X()[0]-Nodes()[0]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1]-Nodes()[0]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2]-Nodes()[0]->X()[2];
	  #else
	  xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
    #endif

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas
  Epetra_SerialDenseMatrix* M_GP;   // EAS matrix M at all GPs
  Epetra_SerialDenseMatrix M;       // EAS matrix M at current GP
  Epetra_SerialDenseVector feas;    // EAS portion of internal forces
  Epetra_SerialDenseMatrix Kaa;     // EAS matrix Kaa
  Epetra_SerialDenseMatrix Kda;     // EAS matrix Kda
  double detJ0;                     // detJ(origin)
  Epetra_SerialDenseMatrix T0invT;  // trafo matrix
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history
  if (eastype_ != soh8_easnone) {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get old alpha
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
     Epetra_SerialDenseVector res_d(NUMDOF_SOH8);
     for (int i = 0; i < NUMDOF_SOH8; ++i) {
       res_d(i) = residual[i];
     }
     // add Kda . res_d to feas
     (*oldfeas).Multiply('N','N',1.0,(*oldKda),res_d,1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);
    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.Size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.Shape(neas_,neas_);

    // EAS matrix K_{d alpha}
    Kda.Shape(neas_,NUMDOF_SOH8);

    // transformation matrix T0, maps M-matrix evaluated at origin
    // between local element coords and global coords
    // here we already get the inverse transposed T0
    T0invT.Shape(NUMSTR_SOH8,NUMSTR_SOH8);

    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    soh8_eassetup(&M_GP,detJ0,T0invT,xrefe);
  } // -------------------------------------------------------------------- EAS

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    Epetra_SerialDenseMatrix jac(NUMDIM_SOH8,NUMDIM_SOH8);
    jac.Multiply('N','N',1.0,int_hex8.deriv_gp[gp],xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    #if 0 // SIMPLE DOUBLE DET COMPUTATION
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    #else
    double detJ= (long double)jac(0,0) * (long double)jac(1,1) * (long double)jac(2,2)
               + (long double)jac(0,1) * (long double)jac(1,2) * (long double)jac(2,0)
               + (long double)jac(0,2) * (long double)jac(1,0) * (long double)jac(2,1)
               - (long double)jac(0,0) * (long double)jac(1,2) * (long double)jac(2,1)
               - (long double)jac(0,1) * (long double)jac(1,0) * (long double)jac(2,2)
               - (long double)jac(0,2) * (long double)jac(1,1) * (long double)jac(2,0);
    #endif

    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_SOH8,NUMNOD_SOH8);
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    /* (unfortunately) we need a redundant copy of derivatives matrix
     *  as the solver factorisation changes it */
    Epetra_SerialDenseMatrix temp_deriv(int_hex8.deriv_gp[gp]);
    solve_for_inverseJac.SetVectors(N_XYZ,temp_deriv);// set X=N_XYZ, B=deriv_gp
    solve_for_inverseJac.FactorWithEquilibration(true);
    //solve_for_inverseJac.EquilibrateMatrix();
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if ((err != 0) && (err2!=0)) dserror("Inversion of Jacobian failed");

	/*cout << deriv_gp;
    cout << int_hex8.deriv_gp[gp];
    getchar();*/

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
    defgrd.Multiply('T','T',1.0,xcurr,N_XYZ,1.0);

    // Right Cauchy-Green tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOH8,NUMDIM_SOH8);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,1.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOH8);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

//    // debugging of EAS with hyperpoly
//    Epetra_SerialDenseMatrix cmat_disp(NUMSTR_SOH8,NUMSTR_SOH8);
//    Epetra_SerialDenseVector stress_disp(NUMSTR_SOH8);
//    double density_disp;
//    const int ele_ID_disp = Id();
//    soh8_mat_sel(&stress_disp,&cmat_disp,&density_disp,&glstrain, &defgrd, gp, ele_ID_disp, time);
//    Epetra_SerialDenseVector glstrain_disp(glstrain);
//    // debugging of EAS with hyperpoly

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone) {
      // get EAS matrix M at current gausspoint gp
      M.Shape(NUMSTR_SOH8,neas_);
      for (int m=0; m<NUMSTR_SOH8; ++m) {
        for (int n=0; n<neas_; ++n) {
          M(m,n)=(*M_GP)(NUMSTR_SOH8*gp+m,n);
        }
      }
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      M.Multiply('N','N',detJ0/detJ,T0invT,Mtemp,0.0);
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      glstrain.Multiply('N','N',1.0,M,(*alpha),1.0);
    } // ------------------------------------------------------------------ EAS


    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    Epetra_SerialDenseMatrix bop(NUMSTR_SOH8,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) {
      bop(0,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    double density;
    const int ele_ID = Id();
    soh8_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp, ele_ID, time);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (stresses != NULL){
      for (int i = 0; i < NUMSTR_SOH8; ++i) {
        (*stresses)(gp,i) = stress(i);
      }
    }

//    // debugging of EAS with hyperpoly
//    Epetra_SerialDenseVector DE(glstrain);
//    DE.Scale(-1.0);
//    DE += glstrain_disp;
//    Epetra_SerialDenseVector DS(stress);
//    DS.Scale(-1.0);
//    DS += stress_disp;
//    //cout << endl << "Delta E: " << DE;
//    //cout << endl << "Delta S: " << DS;

    if (force != NULL || stiffmatrix != NULL){
      if (force != NULL) {
        // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
        (*force).Multiply('T','N',detJ * int_hex8.weights(gp),bop,stress,1.0);
      }

      Epetra_SerialDenseMatrix cb(NUMSTR_SOH8,NUMDOF_SOH8);
      cb.Multiply('N','N',1.0,cmat,bop,1.0);          // temporary C . Bl

      if (stiffmatrix != NULL) {
        // integrate `elastic' and `initial-displacement' stiffness matrix
        // keu = keu + (B^T . C . B) * detJ * w(gp)
        (*stiffmatrix).Multiply('T','N',detJ * int_hex8.weights(gp),bop,cb,1.0);

        // integrate `geometric' stiffness matrix and add to keu *****************
        Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
        sfac.Scale(detJ * int_hex8.weights(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
        vector<double> SmB_L(NUMDIM_SOH8);     // intermediate Sm.B_L
        // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
        for (int inod=0; inod<NUMNOD_SOH8; ++inod){
          SmB_L[0] = sfac(0) * N_XYZ(0,inod) + sfac(3) * N_XYZ(1,inod) + sfac(5) * N_XYZ(2,inod);
          SmB_L[1] = sfac(3) * N_XYZ(0,inod) + sfac(1) * N_XYZ(1,inod) + sfac(4) * N_XYZ(2,inod);
          SmB_L[2] = sfac(5) * N_XYZ(0,inod) + sfac(4) * N_XYZ(1,inod) + sfac(2) * N_XYZ(2,inod);
          for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod){
            double bopstrbop = 0.0;            // intermediate value
            for (int idim=0; idim<NUMDIM_SOH8; ++idim) bopstrbop += N_XYZ(idim,jnod) * SmB_L[idim];
            (*stiffmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += bopstrbop;
            (*stiffmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += bopstrbop;
            (*stiffmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += bopstrbop;
          }
        } // end of integrate `geometric' stiffness******************************
      }

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone) {
        double integrationfactor = detJ * int_hex8.weights(gp);
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        Epetra_SerialDenseMatrix cM(NUMSTR_SOH8,neas_); // temporary c . M
        cM.Multiply('N','N',1.0,cmat,M,0.0);
        Kaa.Multiply('T','N',integrationfactor,M,cM,1.0);

        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        Kda.Multiply('T','N',integrationfactor,M,cb,1.0);

        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        feas.Multiply('T','N',integrationfactor,M,stress,1.0);
      } // ---------------------------------------------------------------- EAS
    }

    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate concistent mass matrix
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double massfactor = (int_hex8.shapefct_gp[gp])(inod) * density * (int_hex8.shapefct_gp[gp])(jnod)
                            * detJ * int_hex8.weights(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  if (force != NULL || stiffmatrix != NULL) {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ != soh8_easnone) {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();

      Epetra_SerialDenseMatrix KdaKaa(NUMDOF_SOH8,neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('T','N',1.0,Kda,Kaa,1.0);

      if (stiffmatrix != NULL) {
        // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
        (*stiffmatrix).Multiply('N','N',-1.0,KdaKaa,Kda,1.0);
      }

      if (force != NULL) {
        // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
        (*force).Multiply('N','N',-1.0,KdaKaa,feas,1.0);
      }

      // store current EAS data in history
      for (int i=0; i<neas_; ++i)
      {
        for (int j=0; j<neas_; ++j) (*oldKaainv)(i,j) = Kaa(i,j);
        for (int j=0; j<NUMDOF_SOH8; ++j) (*oldKda)(i,j) = Kda(i,j);
        (*oldfeas)(i,0) = feas(i);
      }
    } // -------------------------------------------------------------------- EAS
    if (stiffmatrix != NULL){
      SymmetriseMatrix(*stiffmatrix);
    }
  }
  return;
} // DRT::ELEMENTS::Shell8::s8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex8                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_shapederiv(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseMatrix** deriv,     // pointer to pointer of derivs
      Epetra_SerialDenseVector** weights)   // pointer to pointer of weights
{
  // static matrix objects, kept in memory
  static Epetra_SerialDenseMatrix  f(NUMNOD_SOH8,NUMGPT_SOH8);  // shape functions
  static Epetra_SerialDenseMatrix df(NUMDOF_SOH8,NUMNOD_SOH8);  // derivatives
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOH8);   // weights for each gp
  static bool fdf_eval;                      // flag for re-evaluate everything

  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double gpw      = 1.0;              // weight at every gp for linear fct

  if (fdf_eval==true){             // if true f,df already evaluated
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv    = &df;            // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    return;
  } else {
    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
    const double w[NUMGPT_SOH8] = {   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw};

    // fill up nodal f at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
      f(0,i) = (1.0-r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      f(1,i) = (1.0+r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      f(2,i) = (1.0+r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      f(3,i) = (1.0-r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      f(4,i) = (1.0-r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      f(5,i) = (1.0+r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      f(6,i) = (1.0+r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      f(7,i) = (1.0-r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      weightfactors[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
    }

    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
        // df wrt to r "+0" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;

        // df wrt to s "+1" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;

        // df wrt to t "+2" for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }

    // return adresses of just evaluated matrices
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv = &df;               // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = true;               // now all arrays are filled statically
  }
  return;
}  // of soh8_shapederiv


bool DRT::ELEMENTS::So_hex8::soh8_checkRewinding()
{
    const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(DRT::UTILS::intrule_hex_1point);
    const double r = intpoints.qxg[0][0];
    const double s = intpoints.qxg[0][1];
    const double t = intpoints.qxg[0][2];

    Epetra_SerialDenseMatrix deriv(NUMDIM_SOH8, NUMNOD_SOH8);
    DRT::UTILS::shape_function_3D_deriv1(deriv, r, s, t, hex8);

    // update element geometry
    Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
    for (int i=0; i<NUMNOD_SOH8; ++i){
      xrefe(i,0) = Nodes()[i]->X()[0];
      xrefe(i,1) = Nodes()[i]->X()[1];
      xrefe(i,2) = Nodes()[i]->X()[2];
    }

      /* compute the Jacobian matrix which looks like:
      **         [ x_,r  y_,r  z_,r ]
      **     J = [ x_,s  y_,s  z_,s ]
      **         [ x_,t  y_,t  z_,t ]
      */
      Epetra_SerialDenseMatrix jac(NUMDIM_SOH8,NUMDIM_SOH8);
      jac.Multiply('N','N',1.0,deriv,xrefe,1.0);

      // compute determinant of Jacobian by Sarrus' rule
      double detJ= jac(0,0) * jac(1,1) * jac(2,2)
                 + jac(0,1) * jac(1,2) * jac(2,0)
                 + jac(0,2) * jac(1,0) * jac(2,1)
                 - jac(0,0) * jac(1,2) * jac(2,1)
                 - jac(0,1) * jac(1,0) * jac(2,2)
                 - jac(0,2) * jac(1,1) * jac(2,0);
      if (abs(detJ) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
      else if (detJ < 0.0) return true;
      else if (detJ > 0.0) return false;
      dserror("chekRewinding failed!");
      return false;
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Soh8Register::Initialize(DRT::Discretization& dis)
{
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_hex8) continue;
    DRT::ELEMENTS::So_hex8* actele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8* failed");

    if (!actele->donerewinding_) {
      actele->rewind_ = actele->soh8_checkRewinding();

      if (actele->rewind_) {
        int new_nodeids[NUMNOD_SOH8];
        const int* old_nodeids;
        old_nodeids = actele->NodeIds();
        // rewinding of nodes to arrive at mathematically positive element
        new_nodeids[0] = old_nodeids[4];
        new_nodeids[1] = old_nodeids[5];
        new_nodeids[2] = old_nodeids[6];
        new_nodeids[3] = old_nodeids[7];
        new_nodeids[4] = old_nodeids[0];
        new_nodeids[5] = old_nodeids[1];
        new_nodeids[6] = old_nodeids[2];
        new_nodeids[7] = old_nodeids[3];
        actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
      }
      // process of rewinding done
      actele->donerewinding_ = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false,false,false);

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
