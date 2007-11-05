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

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1);
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
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1);
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
      sosh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1);
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
      Epetra_SerialDenseVector* force)          // element internal force vector
{
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
  } else if (eastype_ == soh8_easnone){
    cout << "Warning: Solid-Shell8 without EAS" << endl;
  } else dserror("Solid-Shell8 only with eas_sosh8");// ------------------- EAS

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space
  const int num_sp = 8;       // number of ANS sampling points
  const int num_ans = 3;      // number of modified ANS strains (E_rt,E_st,E_tt)
  // ANS modified rows of bop in local(parameter) coords
  Epetra_SerialDenseMatrix B_ans_loc(num_ans*num_sp,NUMDOF_SOH8);
  // Jacobian evaluated at all ANS sampling points
  Epetra_SerialDenseMatrix jac_sps(NUMDIM_SOH8*num_sp,NUMDIM_SOH8);
  // CURRENT Jacobian evaluated at all ANS sampling points
  Epetra_SerialDenseMatrix jac_cur_sps(NUMDIM_SOH8*num_sp,NUMDIM_SOH8);
  // pointer to derivs evaluated at all sampling points
  Epetra_SerialDenseMatrix* deriv_sp; //[NUMDIM_SOH8*numsp][NUMNOD_SOH8]
  // evaluate all necessary variables for ANS
  sosh8_anssetup(num_sp,num_ans,xrefe,xcurr,&deriv_sp,jac_sps,jac_cur_sps,B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
  const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};

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
    jac_cur.Multiply('N','N',1.0,deriv_gp,xcurr,1.0);

    // set up B-Operator in local(parameter) element space including ANS
    Epetra_SerialDenseMatrix bop_loc(NUMSTR_SOH8,NUMDOF_SOH8);
    for (int inode = 0; inode < NUMNOD_SOH8; ++inode) {
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
        // B_loc_rr = N_r.X_r
        bop_loc(0,inode*3+dim) = deriv_gp(0,inode) * jac_cur(0,dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1,inode*3+dim) = deriv_gp(1,inode) * jac_cur(1,dim);
        // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
        //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
        //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
        bop_loc(2,inode*3+dim) = 0.25*(1-r[gp])*(1-s[gp]) * B_ans_loc(0+4*num_ans,inode*3+dim)
                                +0.25*(1+r[gp])*(1-s[gp]) * B_ans_loc(0+5*num_ans,inode*3+dim)
                                +0.25*(1+r[gp])*(1+s[gp]) * B_ans_loc(0+6*num_ans,inode*3+dim)
                                +0.25*(1-r[gp])*(1+s[gp]) * B_ans_loc(0+7*num_ans,inode*3+dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3,inode*3+dim) = deriv_gp(0,inode) * jac_cur(1,dim)
                                +deriv_gp(1,inode) * jac_cur(0,dim);
        // B_loc_st = interpolation along r of ANS B_loc_st
        //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
        bop_loc(4,inode*3+dim) = 0.5*(1.0+r[gp]) * B_ans_loc(1+1*num_ans,inode*3+dim)
                                +0.5*(1.0-r[gp]) * B_ans_loc(1+3*num_ans,inode*3+dim);
        // B_loc_rt = interpolation along s of ANS B_loc_rt
        //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
        bop_loc(5,inode*3+dim) = 0.5*(1.0-s[gp]) * B_ans_loc(2+0*num_ans,inode*3+dim)
                                +0.5*(1.0+s[gp]) * B_ans_loc(2+2*num_ans,inode*3+dim);
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    Epetra_SerialDenseMatrix TinvT(NUMSTR_SOH8,NUMSTR_SOH8);
    sosh8_evaluateT(jac,TinvT);
    Epetra_SerialDenseMatrix bop(NUMSTR_SOH8,NUMDOF_SOH8);
    bop.Multiply('N','N',1.0,TinvT,bop_loc,1.0);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    // but with modified ANS strains E33, E23 and E13
    Epetra_SerialDenseVector lstrain(NUMSTR_SOH8);
    // evaluate glstrains in local(parameter) coords
    // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
    lstrain(0)= 0.5 * (
       +(jac_cur(0,0)*jac_cur(0,0) + jac_cur(0,1)*jac_cur(0,1) + jac_cur(0,2)*jac_cur(0,2))
       -(jac(0,0)*jac(0,0)         + jac(0,1)*jac(0,1)         + jac(0,2)*jac(0,2)));
    // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
    lstrain(1)= 0.5 * (
       +(jac_cur(1,0)*jac_cur(1,0) + jac_cur(1,1)*jac_cur(1,1) + jac_cur(1,2)*jac_cur(1,2))
       -(jac(1,0)*jac(1,0)         + jac(1,1)*jac(1,1)         + jac(1,2)*jac(1,2)));
    // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
    lstrain(3)= (
       +(jac_cur(0,0)*jac_cur(1,0) + jac_cur(0,1)*jac_cur(1,1) + jac_cur(0,2)*jac_cur(1,2))
       -(jac(0,0)*jac(1,0)         + jac(0,1)*jac(1,1)         + jac(0,2)*jac(1,2)));

    // ANS modification of strains ************************************** ANS
    double dydt_A = 0.0; double dYdt_A = 0.0;
    double dxdt_B = 0.0; double dXdt_B = 0.0;
    double dydt_C = 0.0; double dYdt_C = 0.0;
    double dxdt_D = 0.0; double dXdt_D = 0.0;
    double dzdt_E = 0.0; double dZdt_E = 0.0;
    double dzdt_F = 0.0; double dZdt_F = 0.0;
    double dzdt_G = 0.0; double dZdt_G = 0.0;
    double dzdt_H = 0.0; double dZdt_H = 0.0;

    // vector product of rows of jacobians at corresponding sampling point    cout << jac_cur_sps;
    for (int dim = 0; dim < NUMDIM_SOH8; ++dim) {
      dydt_A += jac_cur_sps(0+0*NUMDIM_SOH8,dim) * jac_cur_sps(2+0*NUMDIM_SOH8,dim);
      dYdt_A += jac_sps(0+0*NUMDIM_SOH8,dim)     * jac_sps(2+0*NUMDIM_SOH8,dim);
      dxdt_B += jac_cur_sps(1+1*NUMDIM_SOH8,dim) * jac_cur_sps(2+1*NUMDIM_SOH8,dim);
      dXdt_B += jac_sps(1+1*NUMDIM_SOH8,dim)     * jac_sps(2+1*NUMDIM_SOH8,dim);
      dydt_C += jac_cur_sps(0+2*NUMDIM_SOH8,dim) * jac_cur_sps(2+2*NUMDIM_SOH8,dim);
      dYdt_C += jac_sps(0+2*NUMDIM_SOH8,dim)     * jac_sps(2+2*NUMDIM_SOH8,dim);
      dxdt_D += jac_cur_sps(1+3*NUMDIM_SOH8,dim) * jac_cur_sps(2+3*NUMDIM_SOH8,dim);
      dXdt_D += jac_sps(1+3*NUMDIM_SOH8,dim)     * jac_sps(2+3*NUMDIM_SOH8,dim);

      dzdt_E += jac_cur_sps(2+4*NUMDIM_SOH8,dim) * jac_cur_sps(2+4*NUMDIM_SOH8,dim);
      dZdt_E += jac_sps(2+4*NUMDIM_SOH8,dim)     * jac_sps(2+4*NUMDIM_SOH8,dim);
      dzdt_F += jac_cur_sps(2+5*NUMDIM_SOH8,dim) * jac_cur_sps(2+5*NUMDIM_SOH8,dim);
      dZdt_F += jac_sps(2+5*NUMDIM_SOH8,dim)     * jac_sps(2+5*NUMDIM_SOH8,dim);
      dzdt_G += jac_cur_sps(2+6*NUMDIM_SOH8,dim) * jac_cur_sps(2+6*NUMDIM_SOH8,dim);
      dZdt_G += jac_sps(2+6*NUMDIM_SOH8,dim)     * jac_sps(2+6*NUMDIM_SOH8,dim);
      dzdt_H += jac_cur_sps(2+7*NUMDIM_SOH8,dim) * jac_cur_sps(2+7*NUMDIM_SOH8,dim);
      dZdt_H += jac_sps(2+7*NUMDIM_SOH8,dim)     * jac_sps(2+7*NUMDIM_SOH8,dim);
    }
    // E33: remedy of curvature thickness locking
    // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
    lstrain(2) = 0.5 * (
       0.25*(1-r[gp])*(1-s[gp]) * (dzdt_E - dZdt_E)
      +0.25*(1+r[gp])*(1-s[gp]) * (dzdt_F - dZdt_F)
      +0.25*(1+r[gp])*(1+s[gp]) * (dzdt_G - dZdt_G)
      +0.25*(1-r[gp])*(1+s[gp]) * (dzdt_H - dZdt_H));
    // E23: remedy of transverse shear locking
    // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
    lstrain(4) = 0.5*(1+r[gp]) * (dxdt_B - dXdt_B) + 0.5*(1-r[gp]) * (dxdt_D - dXdt_D);
    // E13: remedy of transverse shear locking
    // Ert = (1-s)/2 * Est(SP A) + (1+s)/2 * Est(SP C)
    lstrain(5) = 0.5*(1-s[gp]) * (dydt_A - dYdt_A) + 0.5*(1+s[gp]) * (dydt_C - dYdt_C);
    // ANS modification of strains ************************************** ANS

    // transformation of local glstrains 'back' to global(material) space
    Epetra_SerialDenseVector glstrain(NUMSTR_SOH8);
    glstrain.Multiply('N','N',1.0,TinvT,lstrain,1.0);

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



    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    double density;
    // Caution!! the defgrd can not be modified with ANS to remedy locking
    // therefore it is empty and passed only for compatibility reasons
    Epetra_SerialDenseMatrix defgrd; // Caution!! empty!!
    const int ele_ID = Id();
    double time = 0.;                // set to 0. because time is
                                     // not needed currently in solid shell
                                     // but in solid for time-dependent
                                     // complex material behavior

    soh8_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp, ele_ID, time);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_SOH8,NUMDOF_SOH8);
    cb.Multiply('N','N',1.0,cmat,bop,1.0);          // temporary C . B
    (*stiffmatrix).Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);

    // intergrate `geometric' stiffness matrix and add to keu *****************
//    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
//    sfac.Scale(detJ * (*weights)(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
//    vector<double> SmB_L(NUMDIM_SOH8);     // intermediate Sm.B_L
//    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<NUMNOD_SOH8; ++inod){
      for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod){
        Epetra_SerialDenseVector G_ij(NUMSTR_SOH8);
        G_ij(0) = deriv_gp(0,inod) * deriv_gp(0,jnod); // rr-dir
        G_ij(1) = deriv_gp(1,inod) * deriv_gp(1,jnod); // ss-dir
        G_ij(3) = deriv_gp(0,inod) * deriv_gp(1,jnod) + deriv_gp(1,inod) * deriv_gp(0,jnod); //rs-dir
        // ANS modification in tt-dir
        G_ij(2) = 0.25*(1-r[gp])*(1-s[gp]) * (*deriv_sp)(2+4*NUMDIM_SOH8,inod) * (*deriv_sp)(2+4*NUMDIM_SOH8,jnod)
                 +0.25*(1+r[gp])*(1-s[gp]) * (*deriv_sp)(2+5*NUMDIM_SOH8,inod) * (*deriv_sp)(2+5*NUMDIM_SOH8,jnod)
                 +0.25*(1+r[gp])*(1+s[gp]) * (*deriv_sp)(2+6*NUMDIM_SOH8,inod) * (*deriv_sp)(2+6*NUMDIM_SOH8,jnod)
                 +0.25*(1-r[gp])*(1+s[gp]) * (*deriv_sp)(2+7*NUMDIM_SOH8,inod) * (*deriv_sp)(2+7*NUMDIM_SOH8,jnod);
        // ANS modification in st-dir
        G_ij(4) = 0.5*((1+r[gp]) * ((*deriv_sp)(1+1*NUMDIM_SOH8,inod) * (*deriv_sp)(2+1*NUMDIM_SOH8,jnod)
                                   +(*deriv_sp)(2+1*NUMDIM_SOH8,inod) * (*deriv_sp)(1+1*NUMDIM_SOH8,jnod))
                      +(1-r[gp]) * ((*deriv_sp)(1+3*NUMDIM_SOH8,inod) * (*deriv_sp)(2+3*NUMDIM_SOH8,jnod)
                                   +(*deriv_sp)(2+3*NUMDIM_SOH8,inod) * (*deriv_sp)(1+3*NUMDIM_SOH8,jnod)));
        // ANS modification in rt-dir
        G_ij(5) = 0.5*((1-s[gp]) * ((*deriv_sp)(0+0*NUMDIM_SOH8,inod) * (*deriv_sp)(2+0*NUMDIM_SOH8,jnod)
                                   +(*deriv_sp)(2+0*NUMDIM_SOH8,inod) * (*deriv_sp)(0+0*NUMDIM_SOH8,jnod))
                      +(1+s[gp]) * ((*deriv_sp)(0+2*NUMDIM_SOH8,inod) * (*deriv_sp)(2+2*NUMDIM_SOH8,jnod)
                                   +(*deriv_sp)(2+2*NUMDIM_SOH8,inod) * (*deriv_sp)(0+2*NUMDIM_SOH8,jnod)));
        // transformation of local(parameter) space 'back' to global(material) space
        Epetra_SerialDenseVector G_ij_glob(NUMSTR_SOH8);
        G_ij_glob.Multiply('N','N',1.0,TinvT,G_ij,0.0);

        // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
        Epetra_SerialDenseVector Gij(1); // this is a scalar
        Gij.Multiply('T','N',detJ * (*weights)(gp),stress,G_ij_glob,0.0);

        // add "geometric part" Gij times detJ*weights to stiffness matrix
        (*stiffmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += Gij(0);
        (*stiffmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += Gij(0);
        (*stiffmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += Gij(0);
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
    SymmetriseMatrix(Kaa);
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(Kaa);
    solve_for_inverseKaa.Invert();

//    cout << "Kda" << Kda;
    Epetra_SerialDenseMatrix KdaTKaa(NUMDOF_SOH8,neas_); // temporary Kda^T.Kaa^{-1}
    KdaTKaa.Multiply('T','N',1.0,Kda,Kaa,1.0);

    Epetra_SerialDenseVector L6(NUMDOF_SOH8);
    //SymmetricEigen(KdaTKaaKda,L6,NUMDOF_SOH8,'N');
    //cout << setprecision(16) << KdaTKaaKda;
    //cout << "eigen(KEAS): " << L6;
    
    // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
    (*stiffmatrix).Multiply('N','N',-1.0,KdaTKaa,Kda,1.0);
    //cout << setprecision(16) << *stiffmatrix;
    
    // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
    (*force).Multiply('N','N',-1.0,KdaTKaa,feas,1.0);

    // store current EAS data in history
    for (int i=0; i<neas_; ++i)
    {
      for (int j=0; j<neas_; ++j) (*oldKaainv)(i,j) = Kaa(i,j);
      for (int j=0; j<NUMDOF_SOH8; ++j) (*oldKda)(i,j) = Kda(i,j);
      (*oldfeas)(i,0) = feas(i);
    }
  } // -------------------------------------------------------------------- EAS
  /*------------------------------------- make estif absolute symmetric */
  SymmetriseMatrix(*stiffmatrix);
//  Epetra_SerialDenseVector L5(NUMDOF_SOH8);
//  SymmetricEigen(*stiffmatrix,L5,NUMDOF_SOH8,'N');
//  cout << "eigen(K): " << L5;
  
//  SymmetriseMatrix(KdaTKaaKda);
//  Epetra_SerialDenseVector L16(NUMDOF_SOH8);
//  Epetra_SerialDenseMatrix newstiff2(newstiff);
//  cout << "newstiff2"<<setprecision(16) << newstiff2;
// 
//  SymmetricEigen(newstiff,L16,NUMDOF_SOH8,'N');
//  cout << "eigen(newstiffpreEAS): " << L16;
//  double normKeas = KdaTKaaKda.NormOne();
//  cout << "KdaTKda"<<setprecision(16) << KdaTKaaKda;
//  Epetra_SerialDenseMatrix KdaTKaaKda2(KdaTKaaKda);
//  Epetra_SerialDenseVector L36(NUMDOF_SOH8);
//  SymmetricEigen(KdaTKaaKda2,L36,NUMDOF_SOH8,'N');
//  cout << "eigen(KEAS): " << L36;
//  newstiff2 += KdaTKaaKda;
//  Epetra_SerialDenseVector L26(NUMDOF_SOH8);
//  SymmetricEigen(newstiff2,L26,NUMDOF_SOH8,'N');
//  cout << "eigen(newstiff2): " << L26;
  
  return;
} // DRT::Elements::Shell8::s8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_sh8::sosh8_anssetup(
          const int numsp,              // number of sampling points
          const int numans,             // number of ans strains
          const Epetra_SerialDenseMatrix& xrefe, // material element coords
          const Epetra_SerialDenseMatrix& xcurr, // current element coords
          Epetra_SerialDenseMatrix** deriv_sp,   // derivs eval. at all sampling points
          Epetra_SerialDenseMatrix& jac_sps,     // jac at all sampling points
          Epetra_SerialDenseMatrix& jac_cur_sps, // current jac at all sampling points
          Epetra_SerialDenseMatrix& B_ans_loc) // modified B
{
  // static matrix object of derivs at sampling points, kept in memory
  static Epetra_SerialDenseMatrix df_sp(NUMDIM_SOH8*numsp,NUMNOD_SOH8);
  static bool dfsp_eval;                      // flag for re-evaluate everything

  if (dfsp_eval!=0){             // if true f,df already evaluated
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
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
    // numsp = 8 here set explicitly to allow direct initializing
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

  // compute Jacobian matrix at all sampling points
  jac_sps.Multiply('N','N',1.0,df_sp,xrefe,1.0);

  // compute CURRENT Jacobian matrix at all sampling points
  jac_cur_sps.Multiply('N','N',1.0,df_sp,xcurr,1.0);

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
        deriv_asp(m,n)=df_sp(NUMDIM_SOH8*sp+m,n);
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
        B_ans_loc(sp*numans+0,inode*3+dim) = deriv_asp(2,inode)*jac_cur(2,dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp*numans+1,inode*3+dim) = deriv_asp(1,inode)*jac_cur(2,dim)
                                            +deriv_asp(2,inode)*jac_cur(1,dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp*numans+2,inode*3+dim) = deriv_asp(0,inode)*jac_cur(2,dim)
                                            +deriv_asp(2,inode)*jac_cur(0,dim);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  evaluate 'T'-transformation matrix )                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_sh8::sosh8_evaluateT(const Epetra_SerialDenseMatrix jac,
                                            Epetra_SerialDenseMatrix& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0,0) = jac(0,0) * jac(0,0);
  TinvT(1,0) = jac(1,0) * jac(1,0);
  TinvT(2,0) = jac(2,0) * jac(2,0);
  TinvT(3,0) = 2 * jac(0,0) * jac(1,0);
  TinvT(4,0) = 2 * jac(1,0) * jac(2,0);
  TinvT(5,0) = 2 * jac(0,0) * jac(2,0);
  
  TinvT(0,1) = jac(0,1) * jac(0,1);
  TinvT(1,1) = jac(1,1) * jac(1,1);
  TinvT(2,1) = jac(2,1) * jac(2,1);
  TinvT(3,1) = 2 * jac(0,1) * jac(1,1);
  TinvT(4,1) = 2 * jac(1,1) * jac(2,1);
  TinvT(5,1) = 2 * jac(0,1) * jac(2,1);

  TinvT(0,2) = jac(0,2) * jac(0,2);
  TinvT(1,2) = jac(1,2) * jac(1,2);
  TinvT(2,2) = jac(2,2) * jac(2,2);
  TinvT(3,2) = 2 * jac(0,2) * jac(1,2);
  TinvT(4,2) = 2 * jac(1,2) * jac(2,2);
  TinvT(5,2) = 2 * jac(0,2) * jac(2,2);
  
  TinvT(0,3) = jac(0,0) * jac(0,1);
  TinvT(1,3) = jac(1,0) * jac(1,1);
  TinvT(2,3) = jac(2,0) * jac(2,1);
  TinvT(3,3) = jac(0,0) * jac(1,1) + jac(1,0) * jac(0,1);
  TinvT(4,3) = jac(1,0) * jac(2,1) + jac(2,0) * jac(1,1);
  TinvT(5,3) = jac(0,0) * jac(2,1) + jac(2,0) * jac(0,1);


  TinvT(0,4) = jac(0,1) * jac(0,2);
  TinvT(1,4) = jac(1,1) * jac(1,2);
  TinvT(2,4) = jac(2,1) * jac(2,2);
  TinvT(3,4) = jac(0,1) * jac(1,2) + jac(1,1) * jac(0,2);
  TinvT(4,4) = jac(1,1) * jac(2,2) + jac(2,1) * jac(1,2);
  TinvT(5,4) = jac(0,1) * jac(2,2) + jac(2,1) * jac(0,2);
  
  TinvT(0,5) = jac(0,0) * jac(0,2);
  TinvT(1,5) = jac(1,0) * jac(1,2);
  TinvT(2,5) = jac(2,0) * jac(2,2);
  TinvT(3,5) = jac(0,0) * jac(1,2) + jac(1,0) * jac(0,2);
  TinvT(4,5) = jac(1,0) * jac(2,2) + jac(2,0) * jac(1,2);
  TinvT(5,5) = jac(0,0) * jac(2,2) + jac(2,0) * jac(0,2);

  // now evaluate T^{-T} with solver
  Epetra_SerialDenseSolver solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();        
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Tinv (Jacobian) failed");
  return;
}

/*----------------------------------------------------------------------*
 |  find shell-thickness direction via Jacobian                maf 07/07|
 *----------------------------------------------------------------------*/
DRT::Elements::So_sh8::ThicknessDirection DRT::Elements::So_sh8::sosh8_findthickdir() 
{
  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8); // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i) {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }

  // vector of df(origin)
  double df0_vector[NUMDOF_SOH8*NUMNOD_SOH8] =
               {-0.125,-0.125,-0.125,
                +0.125,-0.125,-0.125,
                +0.125,+0.125,-0.125,
                -0.125,+0.125,-0.125,
                -0.125,-0.125,+0.125,
                +0.125,-0.125,+0.125,
                +0.125,+0.125,+0.125,
                -0.125,+0.125,+0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  Epetra_DataAccess CV = Copy;
  Epetra_SerialDenseMatrix df0(CV, df0_vector, NUMDIM_SOH8, NUMDIM_SOH8, 
  NUMNOD_SOH8);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  Epetra_SerialDenseMatrix jac0(NUMDIM_SOH8,NUMDIM_SOH8);
  jac0.Multiply('N', 'N', 1.0, df0, xrefe, 0.0);

  // compute norm of dX_i/dr and dX_i/ds and dX_i/dt
  double dX_dr = 0.0;
  double dX_ds = 0.0;
  double dX_dt = 0.0;
  for (int col = 0; col < NUMDIM_SOH8; ++col) {
    dX_dr += jac0(0, col) * jac0(0, col);
    dX_ds += jac0(1, col) * jac0(1, col);
    dX_dt += jac0(2, col) * jac0(2, col);
  }
  double min_ar = min(dX_dt, min(dX_dr, dX_ds)); // minimum aspect ratio
  
  ThicknessDirection thickdir; // of actual element
  
  if (min_ar == dX_dr) {
    if ((min_ar / dX_ds >= 0.5) || (min_ar / dX_dt >=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autox;
  }
  else if (min_ar == dX_ds) {
    if ((min_ar / dX_dr >= 0.5) || (min_ar / dX_dt >=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autoy;
  }
  else if (min_ar == dX_dt) {
    if ((min_ar / dX_dr >= 0.5) || (min_ar / dX_ds >=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autoz;
  }
  return thickdir;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Sosh8Register::Initialize(DRT::Discretization& dis)
//int DRT::Elements::So_sh8::Initialize_numbers(DRT::Discretization& dis)
{
  //-------------------- loop all my column elements and define thickness direction
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::Elements::So_sh8* actele = dynamic_cast<DRT::Elements::So_sh8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");
    
    if (!actele->nodes_rearranged_){
    // check for automatic definition of thickness direction
    if (actele->thickdir_ == DRT::Elements::So_sh8::autoj) {
      actele->thickdir_ = actele->sosh8_findthickdir();
    }
    
    //int orig_nodeids[NUMNOD_SOH8];
    int new_nodeids[NUMNOD_SOH8];

    switch (actele->thickdir_) {
      case DRT::Elements::So_sh8::autox:
      case DRT::Elements::So_sh8::globx:{
//        cout << endl << "thickness direction is X" << endl;
//        for (int node = 0; node < NUMNOD_SOH8; ++node) {
//          orig_nodeids[node] = actele->NodeIds()[node];
//          cout << node << ": " << orig_nodeids[node] << " inputID: " << actele->inp_nodeIds_[node]<< endl;
//        }
        // resorting of nodes to arrive at local t-dir for global x-dir
        new_nodeids[0] = actele->inp_nodeIds_[7];
        new_nodeids[1] = actele->inp_nodeIds_[4];
        new_nodeids[2] = actele->inp_nodeIds_[0];
        new_nodeids[3] = actele->inp_nodeIds_[3];
        new_nodeids[4] = actele->inp_nodeIds_[6];
        new_nodeids[5] = actele->inp_nodeIds_[5];
        new_nodeids[6] = actele->inp_nodeIds_[1];
        new_nodeids[7] = actele->inp_nodeIds_[2];
        break;
      }
      case DRT::Elements::So_sh8::autoy:
      case DRT::Elements::So_sh8::globy:{
//        for (int node = 0; node < NUMNOD_SOH8; ++node) {
//          orig_nodeids[node] = actele->NodeIds()[node];
//          cout << node << ": " << orig_nodeids[node] << " inputID: " << actele->inp_nodeIds_[node]<< endl;
//        }
        // resorting of nodes to arrive at local t-dir for global y-dir
        new_nodeids[0] = actele->inp_nodeIds_[4];
        new_nodeids[1] = actele->inp_nodeIds_[5];
        new_nodeids[2] = actele->inp_nodeIds_[1];
        new_nodeids[3] = actele->inp_nodeIds_[0];
        new_nodeids[4] = actele->inp_nodeIds_[7];
        new_nodeids[5] = actele->inp_nodeIds_[6];
        new_nodeids[6] = actele->inp_nodeIds_[2];
        new_nodeids[7] = actele->inp_nodeIds_[3];
        break;
      }
      case DRT::Elements::So_sh8::autoz:
      case DRT::Elements::So_sh8::globz:{ 
        // no resorting necessary
        for (int node = 0; node < 8; ++node) {
          new_nodeids[node] = actele->inp_nodeIds_[node];
//          orig_nodeids[node] = actele->NodeIds()[node];
//          cout << node << ": " << orig_nodeids[node] << " inputID: " << actele->inp_nodeIds_[node]<< endl;
        }
        break;
      }
      default:
      dserror("no thickness direction for So_sh8");
    }
    actele->SetNodeIds(NUMNOD_SOH8,new_nodeids);
    actele->nodes_rearranged_ = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false,false,false);
  
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
