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
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

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

extern "C"
{
#include "../headers/standardtypes.h"
// see if we can avoid this #include "../shell8/shell8.h"
}
#include "../drt_lib/dstrc.H"
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_shw6::soshw6_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
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

    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_WEG6,NUMNOD_WEG6);
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(N_XYZ,deriv_gp);// set X=N_XYZ, B=deriv_gp
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();        
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if ((err != 0) && (err2!=0)) dserror("Inversion of Jacobian failed");

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    Epetra_SerialDenseMatrix defgrd(NUMDIM_WEG6,NUMDIM_WEG6);
    defgrd.Multiply('T','T',1.0,xcurr,N_XYZ,1.0);

    // Right Cauchy-Green tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_WEG6,NUMDIM_WEG6);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,1.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_WEG6);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);


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
    Epetra_SerialDenseMatrix bop(NUMSTR_WEG6,NUMDOF_WEG6);
    for (int i=0; i<NUMNOD_WEG6; ++i) {
      bop(0,NODDOF_WEG6*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_WEG6*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_WEG6*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_WEG6*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_WEG6*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_WEG6*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_WEG6*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_WEG6*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_WEG6*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_WEG6*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_WEG6*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_WEG6*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_WEG6*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_WEG6*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_WEG6*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_WEG6*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_WEG6*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_WEG6*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
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

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_WEG6,NUMDOF_WEG6);
    cb.Multiply('N','N',1.0,cmat,bop,1.0);          // temporary C . B
    (*stiffmatrix).Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
    sfac.Scale(detJ * (*weights)(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
    vector<double> SmB_L(NUMDIM_WEG6);     // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<NUMNOD_WEG6; ++inod){
      SmB_L[0] = sfac(0) * N_XYZ(0,inod) + sfac(3) * N_XYZ(1,inod) + sfac(5) * N_XYZ(2,inod);
      SmB_L[1] = sfac(3) * N_XYZ(0,inod) + sfac(1) * N_XYZ(1,inod) + sfac(4) * N_XYZ(2,inod);
      SmB_L[2] = sfac(5) * N_XYZ(0,inod) + sfac(4) * N_XYZ(1,inod) + sfac(2) * N_XYZ(2,inod);
      for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod){
        double bopstrbop = 0.0;            // intermediate value
        for (int idim=0; idim<NUMDIM_WEG6; ++idim) bopstrbop += N_XYZ(idim,jnod) * SmB_L[idim];
        (*stiffmatrix)(NUMDIM_WEG6*inod+0,NUMDIM_WEG6*jnod+0) += bopstrbop;
        (*stiffmatrix)(NUMDIM_WEG6*inod+1,NUMDIM_WEG6*jnod+1) += bopstrbop;
        (*stiffmatrix)(NUMDIM_WEG6*inod+2,NUMDIM_WEG6*jnod+2) += bopstrbop;
      }
    } // end of integrate `geometric' stiffness ******************************


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
} // DRT::Elements::Shell8::s8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_shw6::soshw6_anssetup(
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
        DRT::Utils::shape_function_3D_deriv1(deriv, r[i], s[i], t[i], wedge6);
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
    Epetra_SerialDenseMatrix deriv_asp(NUMDIM_WEG6,numsp);
    for (int m=0; m<NUMDIM_WEG6; ++m) {
      for (int n=0; n<numsp; ++n) {
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




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WEG6
