/*!----------------------------------------------------------------------*##
\file so_hex8_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET10
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_tet10.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "so_tet10_evaluate.h"

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




/*----------------------------------------------------------------------***+
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
//int DRT::Elements::So_hex8::Evaluate(ParameterList& params,
int DRT::Elements::So_tet10::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("So_tet10::Evaluate");

//   #ifdef TET_NO_IMPLEMENT //not yet implemented
  // start with "none"
  DRT::Elements::So_tet10::ActionType act = So_tet10::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_tet10::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_tet10::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_tet10::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_tet10::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_tet10::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_tet10::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_tet10::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_tet10::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_tet10::calc_struct_update_istep;
  else dserror("Unknown type of action for So_tet10");

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
      sotet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);//*
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
      sotet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);//*
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
      sotet10_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);//*
    }
    break;

    // evaluate stresses
    case calc_struct_stress: {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      Epetra_SerialDenseMatrix stresses(NUMGPT_SOTET10,NUMSTR_SOTET10);//*
      sotet10_stress(actmat,mydisp,&stresses); //*
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
  

  dserror("warninng:So_tet10::Evaluate not fully implemented");
  return 0;
}










/*----------------------------------------------------------------------**####
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_tet10::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  DSTraceHelper dst("So_tet10::EvaluateNeumann");
   #ifdef TET_NO_IMPLEMENT //not yet implemented

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
    curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **

// ============================================================================
// CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS
// ============================================================================
// pointer to (static) shape function array
//  for each node, evaluated at each gp
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
  for (int i=0; i<NUMNOD_SOH8; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }

  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // get submatrix of deriv at actual gp
    Epetra_SerialDenseMatrix deriv_gp(NUMDIM_SOH8,NUMGPT_SOH8);
    for (int m=0; m<NUMDIM_SOH8; ++m) {
      for (int n=0; n<NUMGPT_SOH8; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_SOH8*gp+m,n);
      }
    }

    // compute the Jacobian matrix
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

    double fac = (*weights)(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
        elevec1[nodid+dim] += (*shapefct)(nodid,gp) * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }

  }/* ==================================================== end of Loop over GP */
  #endif //TET_NO_IMPLEMENT //not yet implemented
   dserror("");
  return 0;
} // DRT::Elements::Shell8::s8_EvaluateNeumann






/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_tet10::sot10_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      struct _MATERIAL*          material)       // element material data
{
  DSTraceHelper dst("So_tet10::sotet10_nlnstiffmass");

/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
** =============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_SOTET10][NUMGPT_SOTET10]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_SOTET10*NUMDIM][NUMNOD_SOTET10]
/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOTET10]
  soh8_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOTET10,NUMDIM_SOTET10);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOTET10,NUMDIM_SOTET10);  // current  coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_10  Y_10  Z_10 ]
    */
  
  for (int i=0; i<NUMNOD_SOTET10; ++i){  					//**+
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOTET10+2];
  }


  
  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET10; ++gp) {

    // get submatrix of deriv at actual gp
    Epetra_SerialDenseMatrix deriv_gp(NUMNOD_SOTET10 ,NUMCOORD_SOTET10);
    /*for (int m=0; m<NUMDIM_SOTET10; ++m) {
      for (int n=0; n<NUMGPT_SOTET10; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_SOTET10*gp+m,n);
      }
    }*/
    
    for (int shape_func=0; shape_func<NUMNOD_SOTET10; ++shape_func) {
      for (int ksi_num=0;ksi_num<NUMCOORD_SOTET10; ++ksi_num) {
        deriv_gp(shape_func,ksi_num)=(*deriv)(NUMNOD_SOTET10*ksi_num+shape_func,gp);
      }
    }

    /* compute the Jacobian matrix which looks like:
    **         [  1        1        1  	    1      ]
    **     J = [ x_,xsi1  x_,ksi2  x_,ksi3  x,ksi4 ]
    **         [ y_,ksi1  y_,ksi2  y_,ksi3  y,ksi4 ]
    * 		     [ z_,ksi1  z_,ksi2  z_,ksi3  z,ksi4 ]
    */
    Epetra_SerialDenseMatrix jac_temp(NUMCOORD_SOTET10-1,NUMCOORD_SOTET10);
    jac_temp.Multiply('T','N',1.0,xrefe,deriv_gp,1.0);
    Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
    for (int i=0; i<4; i++) jac(0,i)=1;
    for (int row=0;row<3;row++)
    {
    	for (int col=0;col<4;col++)
    	jac(row+1,col)=jac_temp(row,col);	
    }
    
    ~jac_temp;    

    // compute determinant of Jacobian by Sarrus' rule	//*################
   /* double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);*/
    double detJ=det_volf(jac);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_SOTET10,NUMNOD_SOTET10);
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

    /* non-linear B-operator (may be so called, meaning
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
    soh8_mat_sel(material,&stress,&cmat,&density,&glstrain, &defgrd, gp);
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
  
    #endif //TET_NO_IMPLEMENT //not yet implemented
  dserror("DRT::Elements::So_tet10::sot10_nlnstiffmass not implemented yet ");
  
  
  return;
} // DRT::Elements::So_tet10::sot10_nlnstiffmass







/*----------------------------------------------------------------------****+
 |  shape functions and derivatives for So_tet10              volf 06/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_tet10::sot10_shapederiv(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseMatrix** deriv,     // pointer to pointer of derivs
      Epetra_SerialDenseVector** weights)   // pointer to pointer of weights
{
  DSTraceHelper dst("So_tet10::sotet10_shapederiv");

  #ifdef TET_NO_IMPLEMENT //not yet implemented
  // static matrix objects, kept in memory
  static Epetra_SerialDenseMatrix  f(NUMNOD_SOTET10,NUMGPT_SOTET10);  // shape functions
  static Epetra_SerialDenseMatrix df(NUMDOF_SOTET10,NUMNOD_SOTET10);  // derivatives
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOTET10);   // weights for each gp
  static bool fdf_eval;                      // flag for re-evaluate everything

  //const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  //Quadrature rule from Carlos A. Felippa: Adv. FEM  §16.4 
  const double gploc_alpha    = (5.0 +3.0*sqrt(5))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 -    sqrt(5))/20.0; 
  
  //const double gpw      = 1.0;              // weight at every gp for linear fct
  const double gpw      = 0.25;              // weight at every gp for linear fct

  if (fdf_eval==true){             // if true f,df already evaluated
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv    = &df;            // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    return;
  } else {
    // (ksi1, ksi2, ksi3 ,ksi4) gp-locations of fully integrated linear 10-node Tet
    const double ksi1[NUMGPT_SOTET10] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
    const double ksi2[NUMGPT_SOTET10] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
    const double ksi3[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
    const double ksi4[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};
    const double w[NUMGPT_SOTET10] = {   gpw,   gpw,   gpw,   gpw};

    // fill up nodal f at each gp 
    for (int i=0; i<NUMGPT_SOTET10; ++i) {
      f(0,i) = ksi1[i] * (2*ksi1[i] -1);
      f(1,i) = ksi2[i] * (2*ksi2[i] -1);
      f(2,i) = ksi3[i] * (2*ksi3[i] -1);
      f(3,i) = ksi4[i] * (2*ksi4[i] -1);
      f(4,i) = 4 * ksi1[i] * ksi2[i];
      f(5,i) = 4 * ksi2[i] * ksi3[i];
      f(6,i) = 4 * ksi3[i] * ksi1[i];
      f(7,i) = 4 * ksi1[i] * ksi4[i];  
      f(8,i) = 4 * ksi2[i] * ksi4[i];
      f(9,i) = 4 * ksi3[i] * ksi4[i];
      weightfactors[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
    } 

    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOTET10; ++i) {
        // df wrt to ksi1 "+0" for each node(0..9) at each gp [i]
        df(NUMCOORD_SOTET10*0+0,i) = 4 * ksi1[i]-1;
   		df(NUMCOORD_SOTET10*1+0,i) = 0;
      	df(NUMCOORD_SOTET10*2+0,i) = 0;
     	df(NUMCOORD_SOTET10*3+0,i) = 0;
    	df(NUMCOORD_SOTET10*4+0,i) = 4 * ksi2[i];
      	df(NUMCOORD_SOTET10*5+0,i) = 0;
      	df(NUMCOORD_SOTET10*6+0,i) = 4 * ksi3[i];
      	df(NUMCOORD_SOTET10*7+0,i) = 4 * ksi4[i];  
      	df(NUMCOORD_SOTET10*8+0,i) = 0;
      	df(NUMCOORD_SOTET10*9+0,i) = 0;

        // df wrt to ksi2 "+1" for each node(0..9) at each gp [i]
        df(NUMCOORD_SOTET10*0+1,i) = 0;
      	df(NUMCOORD_SOTET10*1+1,i) = 4 * ksi2[i] - 1;
      	df(NUMCOORD_SOTET10*2+1,i) = 0;
      	df(NUMCOORD_SOTET10*3+1,i) = 0;
      	df(NUMCOORD_SOTET10*4+1,i) = 4 * ksi1[i];
      	df(NUMCOORD_SOTET10*5+1,i) = 4 * ksi3[i];
      	df(NUMCOORD_SOTET10*6+1,i) = 0;
      	df(NUMCOORD_SOTET10*7+1,i) = 0;  
      	df(NUMCOORD_SOTET10*8+1,i) = 4 * ksi4[i];
      	df(NUMCOORD_SOTET10*9+1,i) = 0;

        // df wrt to ksi3 "+2" for each node(0..9) at each gp [i]
        df(NUMCOORD_SOTET10*0+2,i) = 0;
      	df(NUMCOORD_SOTET10*1+2,i) = 0;
      	df(NUMCOORD_SOTET10*2+2,i) = 4 * ksi3[i] - 1;
      	df(NUMCOORD_SOTET10*3+2,i) = 0;
      	df(NUMCOORD_SOTET10*4+2,i) = 0;
      	df(NUMCOORD_SOTET10*5+2,i) = 4 * ksi2[i];
      	df(NUMCOORD_SOTET10*6+2,i) = 4 * ksi1[i];
      	df(NUMCOORD_SOTET10*7+2,i) = 0;  
      	df(NUMCOORD_SOTET10*8+2,i) = 0;
      	df(NUMCOORD_SOTET10*9+2,i) = 4 * ksi4[i];
        
        // df wrt to ksi4 "+2" for each node(0..9) at each gp [i]
        df(NUMCOORD_SOTET10*0+3,i) = 0;
      	df(NUMCOORD_SOTET10*1+3,i) = 0;
      	df(NUMCOORD_SOTET10*2+3,i) = 0;
      	df(NUMCOORD_SOTET10*3+3,i) = 4 * ksi4[i] - 1;
      	df(NUMCOORD_SOTET10*4+3,i) = 0;
     	df(NUMCOORD_SOTET10*5+3,i) = 0;
      	df(NUMCOORD_SOTET10*6+3,i) = 0;
      	df(NUMCOORD_SOTET10*7+3,i) = 4 * ksi1[i];  
      	df(NUMCOORD_SOTET10*8+3,i) = 4 * ksi2[i];
      	df(NUMCOORD_SOTET10*9+3,i) = 4 * ksi3[i];
    }

    // return adresses of just evaluated matrices
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv = &df;               // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = true;               // now all arrays are filled statically
  }
  
  #endif //TET_NO_IMPLEMENT //not yet implemented
  dserror("So_tet10::sotet10_shapederiv(..) not implemented yet ");
  return;
}  // So_tet10::sotet10_shapederiv

int DRT::Elements::Sotet10Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

void cut_volf(&Epetra_SerialDenseMatrix in_matrix, int A_row,int A_col,int B_row,int B_col,&Epetra_SerialDenseMatrix out_matrix)
{
	out_matrix.Reshape(B_row-A_row,B_col-A_col);
	
	for (int i_row=0;i_row < out_matrix.M();i_row++)
	for (int i_col; i_col < out_matrix.N();i_col++)
	{
		out_matrix(i_row,i_col)= in_matrix(i_row+A_row,i_col+A_col);			
	}	
}


double det_volf(&Epetra_SerialDenseMatrix in_matrix)
{
	Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1;in_matrix.N()-1);
	
	if (in_matrix.N()==1)
	{
		return in_matrix(0,0);	
	}	
	else if (in_matrix.N()==2)
	{
		return 	((in_matrix(0,0)*in_matrix(1,1))-(in_matrix(0,1)*in_matrix(1,0)));
	}
	else if (in_matrix.N()>2)
	{
		double out_det=0;
		for (int i_col=0;i_col < in_natrix.N();i_col++)
		{
			for (c_col=0;c_col < i_col;i_col++)
			{
				for(row=1;row<in_natrix.N();row++)
				temp_matrix(row,c_col)=in_matrix(row,c_col);							
			}
		
			for (c_col=i_col+1;c_col <  in_natrix.N();i_col++)
			{
			for(row=1;row<in_natrix.N();row++)
			temp_matrix(row,c_col-1)=in_matrix(row,c_col);	
			}
		out_det+=det(in_matrix);
		}
		return out_det;
	}	
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET10
