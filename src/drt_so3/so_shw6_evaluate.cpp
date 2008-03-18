/*!----------------------------------------------------------------------
\file so_shw6_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_shw6.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_shw6::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::So_weg6::ActionType act = So_weg6::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_weg6::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_weg6::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_weg6::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_weg6::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_weg6::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_weg6::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_weg6::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_weg6::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_weg6::calc_struct_update_istep;
  else if (action=="calc_struct_update_genalpha_imrlike")  act = So_weg6::calc_struct_update_genalpha_imrlike;
  else if (action=="postprocess_stress")        act = So_weg6::postprocess_stress;
  else dserror("Unknown type of action for So_weg6");

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
      soshw6_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,time);
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
      soshw6_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,time);
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
      soshw6_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,time);
    }
    break;

    // evaluate stresses at gauss points
    case calc_struct_stress:{
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Epetra_SerialDenseMatrix stress(NUMGPT_WEG6,NUMSTR_WEG6);
      Epetra_SerialDenseMatrix strain(NUMGPT_WEG6,NUMSTR_WEG6);
      soshw6_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,time);
      cout << "gpstress: " << stress << endl;
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses and strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:{

      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_WEG6,NUMSTR_WEG6);
        soweg6_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_WEG6);

        for (int i=0;i<NUMNOD_WEG6;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz") {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_WEG6; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_WEG6; ++j) {
              (*((*elestress)(i)))[lid] += 0.125 * (*gpstress)(j,i);
            }
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_WEG6,NUMSTR_WEG6);
        soweg6_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_WEG6);

        for (int i=0;i<NUMNOD_WEG6;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_WEG6;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_WEG6; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_WEG6; ++j) {
              (*((*elestress)(i)))[lid] += 0.125 * (*gpstress)(j,i);
            }
          }
        }
      }
      else{
        dserror("unknown type of stress/strain output on element level");
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

    case calc_struct_update_genalpha_imrlike: {
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
void DRT::ELEMENTS::So_shw6::soshw6_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // element stresses
      Epetra_SerialDenseMatrix* elestrain,      // element stresses
      const double              time)           // current absolute time
{

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_WEG6][NUMGPT_WEG6]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_WEG6*NUMDIM][NUMNOD_WEG6]
/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_WEG6]
  sow6_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_WEG6,NUMDIM_WEG6);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_WEG6,NUMDIM_WEG6);  // current  coord. of element
  for (int i=0; i<NUMNOD_WEG6; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_WEG6+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_WEG6+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_WEG6+2];
  }

  for (int i=0; i<NUMDOF_WEG6; ++i) cout << disp[i] << ",";
  cout << endl;

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space
  const int num_sp = 5;       // number of ANS sampling points
  const int num_ans = 3;      // number of modified ANS strains (E_rt,E_st,E_tt)
  // ANS modified rows of bop in local(parameter) coords
  Epetra_SerialDenseMatrix B_ans_loc(num_ans*num_sp,NUMDOF_WEG6);
  // Jacobian evaluated at all ANS sampling points
  Epetra_SerialDenseMatrix jac_sps(NUMDIM_WEG6*num_sp,NUMDIM_WEG6);
  // CURRENT Jacobian evaluated at all ANS sampling points
  Epetra_SerialDenseMatrix jac_cur_sps(NUMDIM_WEG6*num_sp,NUMDIM_WEG6);
  // pointer to derivs evaluated at all sampling points
  Epetra_SerialDenseMatrix* deriv_sp; //[NUMDIM_SOH8*numsp][NUMNOD_SOH8]
  // evaluate all necessary variables for ANS
  soshw6_anssetup(num_sp,num_ans,xrefe,xcurr,&deriv_sp,jac_sps,jac_cur_sps,B_ans_loc);
  const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_wedge_6point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule_);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_WEG6; ++gp) {

    // get submatrix of deriv at actual gp
    Epetra_SerialDenseMatrix deriv_gp(NUMDIM_WEG6,NUMGPT_WEG6);
    for (int m=0; m<NUMDIM_WEG6; ++m) {
      for (int n=0; n<NUMGPT_WEG6; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_WEG6*gp+m,n);
      }
    }

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    Epetra_SerialDenseMatrix jac(NUMDIM_WEG6,NUMDIM_WEG6);
    jac.Multiply('N','N',1.0,deriv_gp,xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    if (abs(detJ) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    Epetra_SerialDenseMatrix jac_cur(NUMDIM_WEG6,NUMDIM_WEG6);
    jac_cur.Multiply('N','N',1.0,deriv_gp,xcurr,1.0);

    // need gp-locations for ANS interpolation
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    //const double t = intpoints.qxg[gp][2]; // not needed
    // set up B-Operator in local(parameter) element space including ANS
    Epetra_SerialDenseMatrix bop_loc(NUMSTR_WEG6,NUMDOF_WEG6);
    for (int inode = 0; inode < NUMNOD_WEG6; ++inode) {
      for (int dim = 0; dim < NUMDIM_WEG6; ++dim) {
        // B_loc_rr = N_r.X_r
        bop_loc(0,inode*3+dim) = deriv_gp(0,inode) * jac_cur(0,dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1,inode*3+dim) = deriv_gp(1,inode) * jac_cur(1,dim);
//        // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
//        //          = (1-r-s) * B_ans(SP C) + r * B_ans(SP D) + s * B_ans(SP E)
//        bop_loc(2,inode*3+dim) = (1-r-s) * B_ans_loc(0+2*num_ans,inode*3+dim)
//                                + r      * B_ans_loc(0+3*num_ans,inode*3+dim)
//                                + s      * B_ans_loc(0+4*num_ans,inode*3+dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3,inode*3+dim) = deriv_gp(0,inode) * jac_cur(1,dim)
                                +deriv_gp(1,inode) * jac_cur(0,dim);
//        // B_loc_st = interpolation along r of ANS B_loc_st
//        //          = r * B_ans(SP B)
//        bop_loc(4,inode*3+dim) = r * B_ans_loc(1+1*num_ans,inode*3+dim);
//        // B_loc_rt = interpolation along s of ANS B_loc_rt
//        //          = s * B_ans(SP A)
//        bop_loc(5,inode*3+dim) = s * B_ans_loc(2+0*num_ans,inode*3+dim);

        // testing without ans:
        bop_loc(2,inode*3+dim) = deriv_gp(2,inode) * jac_cur(2,dim);
        bop_loc(4,inode*3+dim) = deriv_gp(1,inode) * jac_cur(2,dim)
                                +deriv_gp(2,inode) * jac_cur(1,dim);
        bop_loc(5,inode*3+dim) = deriv_gp(0,inode) * jac_cur(2,dim)
                                +deriv_gp(2,inode) * jac_cur(0,dim);
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    Epetra_SerialDenseMatrix TinvT(NUMSTR_WEG6,NUMSTR_WEG6);
    soshw6_evaluateT(jac,TinvT);
    Epetra_SerialDenseMatrix bop(NUMSTR_WEG6,NUMDOF_WEG6);
    bop.Multiply('N','N',1.0,TinvT,bop_loc,1.0);

    // local GL strain vector lstrain={Err,Ess,Ett,2*Ers,2*Est,2*Ert}
    // but with modified ANS strains Ett, Est and Ert
    Epetra_SerialDenseVector lstrain(NUMSTR_WEG6);
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

    // testing without ans:
    lstrain(2)= 0.5 * (
       +(jac_cur(2,0)*jac_cur(2,0) + jac_cur(2,1)*jac_cur(2,1) + jac_cur(2,2)*jac_cur(2,2))
       -(jac(2,0)*jac(2,0)         + jac(2,1)*jac(2,1)         + jac(2,2)*jac(2,2)));
    lstrain(4)= (
       +(jac_cur(1,0)*jac_cur(2,0) + jac_cur(1,1)*jac_cur(2,1) + jac_cur(1,2)*jac_cur(2,2))
       -(jac(1,0)*jac(2,0)         + jac(1,1)*jac(2,1)         + jac(1,2)*jac(2,2)));
    lstrain(5)= (
       +(jac_cur(0,0)*jac_cur(2,0) + jac_cur(0,1)*jac_cur(2,1) + jac_cur(0,2)*jac_cur(2,2))
       -(jac(0,0)*jac(2,0)         + jac(0,1)*jac(2,1)         + jac(0,2)*jac(2,2)));

    // ANS modification of strains ************************************** ANS
    double dydt_A = 0.0; double dYdt_A = 0.0; const int spA = 0;
    double dxdt_B = 0.0; double dXdt_B = 0.0; const int spB = 1;
    double dzdt_C = 0.0; double dZdt_C = 0.0; const int spC = 2;
    double dzdt_D = 0.0; double dZdt_D = 0.0; const int spD = 3;
    double dzdt_E = 0.0; double dZdt_E = 0.0; const int spE = 4;

    const int xdir = 0; // index to matrix x-row, r-row respectively
    const int ydir = 1; // index to matrix y-row, s-row respectively
    const int zdir = 2; // index to matrix z-row, t-row respectively

    // vector product of rows of jacobians at corresponding sampling point
    for (int dim = 0; dim < NUMDIM_WEG6; ++dim) {
      dydt_A += jac_cur_sps(ydir+spA*NUMDIM_WEG6,dim) * jac_cur_sps(dim+spA*NUMDIM_WEG6,zdir);
      dYdt_A += jac_sps(ydir+spA*NUMDIM_WEG6,dim)     * jac_sps(dim+spA*NUMDIM_WEG6,zdir);
      dxdt_B += jac_cur_sps(xdir+spB*NUMDIM_WEG6,dim) * jac_cur_sps(dim+spB*NUMDIM_WEG6,zdir);
      dXdt_B += jac_sps(xdir+spB*NUMDIM_WEG6,dim)     * jac_sps(dim+spB*NUMDIM_WEG6,zdir);

      dzdt_C += jac_cur_sps(zdir+spC*NUMDIM_WEG6,dim) * jac_cur_sps(dim+spC*NUMDIM_WEG6,zdir);
      dZdt_C += jac_sps(zdir+spC*NUMDIM_WEG6,dim)     * jac_sps(dim+spC*NUMDIM_WEG6,zdir);
      dzdt_D += jac_cur_sps(zdir+spD*NUMDIM_WEG6,dim) * jac_cur_sps(dim+spD*NUMDIM_WEG6,zdir);
      dZdt_D += jac_sps(zdir+spD*NUMDIM_WEG6,dim)     * jac_sps(dim+spD*NUMDIM_WEG6,zdir);
      dzdt_E += jac_cur_sps(zdir+spE*NUMDIM_WEG6,dim) * jac_cur_sps(dim+spE*NUMDIM_WEG6,zdir);
      dZdt_E += jac_sps(zdir+spE*NUMDIM_WEG6,dim)     * jac_sps(dim+spE*NUMDIM_WEG6,zdir);
    }
//    // E33: remedy of curvature thickness locking
//    // Ett = 0.5* ( (1-r-s) * Ett(SP C) + r * Ett(SP D) + s * Ett(SP E) )
//    lstrain(2) = 0.5 * ( (1-r-s) * (dzdt_C - dZdt_C)
//                        + r * (dzdt_D - dZdt_D)
//                        + s * (dzdt_E - dZdt_E));
//    // E23: remedy of transverse shear locking
//    // Est = r * Est(SP B)
//    lstrain(4) = r * (dxdt_B - dXdt_B);
//    // E13: remedy of transverse shear locking
//    // Ert = s * Est(SP A)
//    lstrain(5) = s * (dydt_A - dYdt_A);
//    // ANS modification of strains ************************************** ANS

    // transformation of local glstrains 'back' to global(material) space
    Epetra_SerialDenseVector glstrain(NUMSTR_WEG6);
    glstrain.Multiply('N','N',1.0,TinvT,lstrain,1.0);

    // return gp strains (only in case of stress/strain output)
    if (elestress != NULL){
      for (int i = 0; i < NUMSTR_WEG6; ++i) {
        (*elestrain)(gp,i) = glstrain(i);
      }
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_WEG6,NUMSTR_WEG6);
    Epetra_SerialDenseVector stress(NUMSTR_WEG6);
    double density;
    sow6_mat_sel(&stress,&cmat,&density,&glstrain, time);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL){
      for (int i = 0; i < NUMSTR_WEG6; ++i) {
        (*elestress)(gp,i) = stress(i);
        //(*elestress)(gp,i) = glstrain(i);
      }
    }

    if (force != NULL && stiffmatrix != NULL) {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T', 'N', detJ * (*weights)(gp), bop, stress, 1.0);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      Epetra_SerialDenseMatrix cb(NUMSTR_WEG6, NUMDOF_WEG6);
      cb.Multiply('N', 'N', 1.0, cmat, bop, 1.0); // temporary C . B
      (*stiffmatrix).Multiply('T', 'N', detJ * (*weights)(gp), bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      for (int inod=0; inod<NUMNOD_WEG6; ++inod) {
        for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod) {
          Epetra_SerialDenseVector G_ij(NUMSTR_WEG6);
          G_ij(0) = deriv_gp(0, inod) * deriv_gp(0, jnod); // rr-dir
          G_ij(1) = deriv_gp(1, inod) * deriv_gp(1, jnod); // ss-dir
          G_ij(3) = deriv_gp(0, inod) * deriv_gp(1, jnod) + deriv_gp(1, inod)
              * deriv_gp(0, jnod); //rs-dir
          // ANS modification in tt-dir
          G_ij(2) = (1-r-s) * (*deriv_sp)(zdir+spC*NUMDIM_WEG6, inod)
              * (*deriv_sp)(zdir+spC*NUMDIM_WEG6, jnod) + r * (*deriv_sp)(zdir
              +spD*NUMDIM_WEG6, inod) * (*deriv_sp)(zdir+spD*NUMDIM_WEG6, jnod)
              + s * (*deriv_sp)(zdir+spE*NUMDIM_WEG6, inod) * (*deriv_sp)(zdir
                  +spE*NUMDIM_WEG6, jnod);
          // ANS modification in st-dir
          G_ij(4) = 0.5*(r * ((*deriv_sp)(ydir+spB*NUMDIM_WEG6, inod)
              * (*deriv_sp)(zdir+spB*NUMDIM_WEG6, jnod) +(*deriv_sp)(zdir+spB
              *NUMDIM_WEG6, inod) * (*deriv_sp)(ydir+spB*NUMDIM_WEG6, jnod)));
          // ANS modification in rt-dir
          G_ij(5) = 0.5*(s * ((*deriv_sp)(xdir+spA*NUMDIM_WEG6, inod)
              * (*deriv_sp)(zdir+spA*NUMDIM_WEG6, jnod) +(*deriv_sp)(zdir+spA
              *NUMDIM_WEG6, inod) * (*deriv_sp)(xdir+spA*NUMDIM_WEG6, jnod)));
          // transformation of local(parameter) space 'back' to global(material) space
          Epetra_SerialDenseVector G_ij_glob(NUMSTR_WEG6);
          G_ij_glob.Multiply('N', 'N', 1.0, TinvT, G_ij, 0.0);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          Epetra_SerialDenseVector Gij(1); // this is a scalar
          Gij.Multiply('T', 'N', detJ * (*weights)(gp), stress, G_ij_glob, 0.0);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_WEG6*inod+0, NUMDIM_WEG6*jnod+0) += Gij(0);
          (*stiffmatrix)(NUMDIM_WEG6*inod+1, NUMDIM_WEG6*jnod+1) += Gij(0);
          (*stiffmatrix)(NUMDIM_WEG6*inod+2, NUMDIM_WEG6*jnod+2) += Gij(0);
        }
      } // end of intergrate `geometric' stiffness ******************************
    }


    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate concistent mass matrix
      for (int inod=0; inod<NUMNOD_WEG6; ++inod) {
        for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod) {
          double massfactor = (*shapefct)(inod,gp) * density * (*shapefct)(jnod,gp)
                            * detJ * (*weights)(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_WEG6*inod+0,NUMDIM_WEG6*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+1,NUMDIM_WEG6*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+2,NUMDIM_WEG6*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::Shell8::s8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_shw6::soshw6_anssetup(
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
  static Epetra_SerialDenseMatrix df_sp(NUMDIM_WEG6*numsp,NUMNOD_WEG6);
  static bool dfsp_eval;                      // flag for re-evaluate everything

  if (dfsp_eval!=0){             // if true f,df already evaluated
    *deriv_sp = &df_sp;         // return adress of static object to target of pointer
  } else {
  /*====================================================================*/
  /* 6-node wedge Solid-Shell node topology
   * and location of sampling points A to E                             */
  /*--------------------------------------------------------------------*/
  /*
   *                             s
   *                   6        /
   *                 //||\\   /
   *      t        //  ||   \\
   *      ^      //    || /    \\
   *      |    //      E          \\
   *      |  //        ||            \\
   *      |//       /  ||               \\
   *      5===============================6
   *     ||      B      3                 ||
   *     ||    /      // \\               ||
   *     || /       //      \\            ||
   *   - C -  -  -// -  A  -  -\\ -  -   -D  ----> r
   *     ||     //                \\      ||
   *  /  ||   //                     \\   ||
   *     || //                          \\||
   *      1================================2
   *
   */
  /*====================================================================*/
    // (r,s,t) gp-locations of sampling points A,B,C,D,E
    // numsp = 5 here set explicitly to allow direct initializing
    double r[5] = { 0.5, 0.0, 0.0, 1.0, 0.0};
    double s[5] = { 0.0, 0.5, 0.0, 0.0, 1.0};
    double t[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i=0; i<numsp; ++i) {
        Epetra_SerialDenseMatrix deriv(NUMDIM_WEG6, NUMNOD_WEG6);
        DRT::UTILS::shape_function_3D_deriv1(deriv, r[i], s[i], t[i], wedge6);
        for (int inode = 0; inode < NUMNOD_WEG6; ++inode) {
          df_sp(i*NUMDIM_WEG6+0, inode) = deriv(0, inode);
          df_sp(i*NUMDIM_WEG6+1, inode) = deriv(1, inode);
          df_sp(i*NUMDIM_WEG6+2, inode) = deriv(2, inode);
        }
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
    Epetra_SerialDenseMatrix deriv_asp(NUMDIM_WEG6,NUMNOD_WEG6);
    for (int m=0; m<NUMDIM_WEG6; ++m) {
      for (int n=0; n<NUMNOD_WEG6; ++n) {
        deriv_asp(m,n)=df_sp(NUMDIM_WEG6*sp+m,n);
      }
    }
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    Epetra_SerialDenseMatrix jac_cur(NUMDIM_WEG6,NUMDIM_WEG6);
    jac_cur.Multiply('N','N',1.0,deriv_asp,xcurr,1.0);

    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_WEG6; ++inode) {
      for (int dim = 0; dim < NUMDIM_WEG6; ++dim) {
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
void DRT::ELEMENTS::So_shw6::soshw6_evaluateT(const Epetra_SerialDenseMatrix jac,
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



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WEG6
