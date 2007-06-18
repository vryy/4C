/*!----------------------------------------------------------------------
\file so_sh8_evaluate.cpp
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
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_sh8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


extern "C"
{
#include "../headers/standardtypes.h"
// see if we can avoid this #include "../shell8/shell8.h"
}
#include "../drt_lib/dstrc.H"
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_sh8::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("So_sh8::Evaluate");

  // start with "none"
  DRT::Elements::So_hex8::ActionType act = So_hex8::none;

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
  else dserror("Unknown type of action for So_hex8");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);
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
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;

    // evaluate stresses
    case calc_struct_stress: {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      Epetra_SerialDenseMatrix stresses(NUMGPT_SOH8,NUMSTR_SOH8);
      dserror("no stress evaluation yet");
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

    default:
      dserror("Unknown type of action for Solid3");
  }
  return 0;
}




/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_sh8::sosh8_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      MATERIAL*                 material)       // element material data
{
  DSTraceHelper dst("So_sh8::sosh8_nlnstiffmass");

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_SOH8][NUMGPT_SOH8]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_SOH8*NUMDIM][NUMNOD_SOH8]
/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOH8]
  soh8_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

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
  if (eastype_ == soh8_eassosh8) {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    //(*alpha).Shape(neas_,1);
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get old alpha
//    // evaluate current (updated) EAS alphas (from history variables)
//    soh8_easupdate(alpha,disp,residual);
    // get stored EAS history
    //(*oldfeas).Shape(neas_,1);
    //(*oldKaainv).Shape(neas_,neas_);
    //(*oldKda).Shape(neas_,NUMDOF_SOH8);
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");
  
    // we need the displacement at the previous step
    Epetra_SerialDenseVector old_d(NUMDOF_SOH8);
    for (int i=0; i<NUMDOF_SOH8; ++i) old_d(i) = disp[i] - residual[i];
  
    // add Kda . old_d to feas
    (*oldfeas).Multiply('N','N',1.0,(*oldKda),old_d,1.0);
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
  } else dserror("Solid-Shell8 only with eas_sosh8");// ------------------- EAS
  
  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space
  const int num_sp = 8;       // number of ANS sampling points
  const int num_ans = 3;      // number of modified ANS strains (E_rt,E_st,E_tt)
  Epetra_SerialDenseMatrix B_ans_loc;  
  B_ans_loc.Shape(num_ans*num_sp,NUMDOF_SOH8);
  sosh8_anssetup(num_sp,num_ans,xcurr,B_ans_loc);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // get submatrix of deriv at actual gp
    Epetra_SerialDenseMatrix deriv_gp(NUMDIM_SOH8,NUMGPT_SOH8);
    for (int m=0; m<NUMDIM_SOH8; ++m) {
      for (int n=0; n<NUMGPT_SOH8; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_SOH8*gp+m,n);
      }
    }

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    Epetra_SerialDenseMatrix jac(NUMDIM_SOH8,NUMDIM_SOH8);
    jac.Multiply('N','N',1.0,deriv_gp,xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    Epetra_SerialDenseMatrix jac_cur(NUMDIM_SOH8,NUMDIM_SOH8);
    jac.Multiply('N','N',1.0,deriv_gp,xcurr,1.0);
    
    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_SOH8,NUMNOD_SOH8);
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(N_XYZ,deriv_gp);// set X=N_XYZ, B=deriv_gp
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if (err != 0) dserror("Inversion of Jacobian failed");

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
    soh8_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, material, gp);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_SOH8,NUMDOF_SOH8);
    cb.Multiply('N','N',1.0,cmat,bop,1.0);          // temporary C . B
    (*stiffmatrix).Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);

    // intergrate `geometric' stiffness matrix and add to keu *****************
    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
    sfac.Scale(detJ * (*weights)(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
    } // end of intergrate `geometric' stiffness ******************************

    // EAS technology: integrate matrices --------------------------------- EAS
    if (eastype_ != soh8_easnone) {
      double integrationfactor = detJ * (*weights)(gp);
      // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
      Epetra_SerialDenseMatrix cM(NUMSTR_SOH8,neas_); // temporary c . M
      cM.Multiply('N','N',1.0,cmat,M,0.0);
      Kaa.Multiply('T','N',integrationfactor,M,cM,1.0);

      // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
      Kda.Multiply('T','N',integrationfactor,M,cb,1.0);

      // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
      feas.Multiply('T','N',integrationfactor,M,stress,1.0);
    } // ------------------------------------------------------------------ EAS

    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate concistent mass matrix
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double massfactor = (*shapefct)(inod,gp) * density * (*shapefct)(jnod,gp)
                            * detJ * (*weights)(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  // EAS technology: ------------------------------------------------------ EAS
  // subtract EAS matrices from disp-based Kdd to "soften" element
  if (eastype_ != soh8_easnone) {
    // we need the inverse of Kaa
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(Kaa);
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix KdaKaa(NUMDOF_SOH8,neas_); // temporary Kda.Kaa^{-1}
    KdaKaa.Multiply('T','N',1.0,Kda,Kaa,0.0);

    // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
    (*stiffmatrix).Multiply('N','N',-1.0,KdaKaa,Kda,1.0);

    // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
    (*force).Multiply('N','N',-1.0,KdaKaa,feas,1.0);

    // store current EAS data in history
    for (int i=0; i<neas_; ++i)
    {
      for (int j=0; j<neas_; ++j) (*oldKaainv)(i,j) = Kaa(i,j);
      for (int j=0; j<NUMDOF_SOH8; ++j) (*oldKda)(i,j) = Kda(i,j);
      (*oldfeas)(i,0) = feas(i);
    }
  } // -------------------------------------------------------------------- EAS
  return;
} // DRT::Elements::Shell8::s8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_sh8::sosh8_anssetup(
          const int numsp,              // number of sampling points
          const int numans,             // number of ans strains
          const Epetra_SerialDenseMatrix& xcurr, // current element coords
          Epetra_SerialDenseMatrix&  B_ans_loc) // modified B 
{
  // static matrix object of derivs at sampling points, kept in memory
  static Epetra_SerialDenseMatrix df_sp(NUMDIM_SOH8*numsp,NUMNOD_SOH8);
  static bool dfsp_eval;                      // flag for re-evaluate everything
  
  Epetra_SerialDenseMatrix* deriv_sp;

  if (dfsp_eval!=0){             // if true f,df already evaluated
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
    return;
  } else {
  /*====================================================================*/
  /* 8-node hexhedra Solid-Shell node topology
   * and location of sampling points A to H                             */
  /*--------------------------------------------------------------------*/
  /*                      t
   *                      |
   *             4========|================7
   *          // |        |              //||
   *        //   |        |            //  ||
   *      //     |        |   D      //    ||
   *     5=======E=================6       H 
   *    ||       |        |        ||      ||
   *    ||   A   |        o--------||-- C -------s
   *    ||       |       /         ||      ||
   *    F        0----- B ---------G ------3
   *    ||     //     /            ||    //
   *    ||   //     /              ||  //
   *    || //     r                ||//
   *     1=========================2
   *
   */
  /*====================================================================*/
    // (r,s,t) gp-locations of sampling points A,B,C,D,E,F,G,H
    double r[8] = { 0.0, 1.0, 0.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double s[8] = {-1.0, 0.0, 1.0, 0.0,-1.0,-1.0, 1.0, 1.0};
    double t[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i=0; i<numsp; ++i) {
        // df wrt to r "+0" for each node(0..7) at each sp [i]
        df_sp(NUMDIM_SOH8*i+0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;

        // df wrt to s "+1" for each node(0..7) at each sp [i]
        df_sp(NUMDIM_SOH8*i+1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        df_sp(NUMDIM_SOH8*i+1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;

        // df wrt to t "+2" for each node(0..7) at each sp [i]
        df_sp(NUMDIM_SOH8*i+2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        df_sp(NUMDIM_SOH8*i+2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
    dfsp_eval = 1;               // now all arrays are filled statically
  }
  
  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  for (int sp = 0; sp < numsp; ++sp) {
    // get submatrix of deriv_sp at actual sp
    Epetra_SerialDenseMatrix deriv_asp(NUMDIM_SOH8,numsp);
    for (int m=0; m<NUMDIM_SOH8; ++m) {
      for (int n=0; n<numsp; ++n) {
        deriv_asp(m,n)=(*deriv_sp)(NUMDIM_SOH8*sp+m,n);
      }
    }
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    Epetra_SerialDenseMatrix jac_cur(NUMDIM_SOH8,NUMDIM_SOH8);
    jac_cur.Multiply('N','N',1.0,deriv_asp,xcurr,1.0);
    
    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_SOH8; ++inode) {
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp*numans+0,inode*3+dim) = deriv_asp(3,inode)*jac_cur(3,dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp*numans+1,inode*3+dim) = deriv_asp(1,inode)*jac_cur(3,dim)
                                            +deriv_asp(3,inode)*jac_cur(1,dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s                                          
        B_ans_loc(sp*numans+2,inode*3+dim) = deriv_asp(2,inode)*jac_cur(3,dim)
                                            +deriv_asp(3,inode)*jac_cur(2,dim);
      }
    }
  }    
    
  
  return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
