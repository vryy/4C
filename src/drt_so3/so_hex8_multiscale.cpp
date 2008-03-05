/*!----------------------------------------------------------------------
\file so_hex8_multiscale.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_utils.H"
#include "Epetra_SerialDenseSolver.h"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize stresses and density (public)                    lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized stress and
// density for multi-scale analyses by averaging the corresponding
// quantities over the initial volume
// currently no EAS is implemented, but this can easily be
// incorporated as soon as the algorithm reliably works


// dynamic homogenization based on averaging the 1st PK stresses
void DRT::ELEMENTS::So_hex8::soh8_homog(ParameterList&  params,
                                        vector<double>& disp,
                                        const double    time,
                                        vector<double>& residual)
{
  // check whether we only have to calculate the initial density
  bool onlydens = params.get<bool>("onlydens", false);

  double homogdens = 0.;

  double P11 = 0.;
  double P12 = 0.;
  double P13 = 0.;
  double P21 = 0.;
  double P22 = 0.;
  double P23 = 0.;
  double P31 = 0.;
  double P32 = 0.;
  double P33 = 0.;

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

    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_SOH8,NUMNOD_SOH8);
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(N_XYZ,deriv_gp);// set X=N_XYZ, B=deriv_gp
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if ((err != 0) && (err2!=0)) dserror("Inversion of Jacobian failed");

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

    double integrationfactor = detJ * (*weights)(gp);

    homogdens += integrationfactor*density;

    // check out whether we have to compute the average of the first
    // Piola Kirchhoff stresses since we prescribe a deformation
    // gradient on the microlevel thus implying that the macroscopic
    // deformation gradient is the average of the microscopic
    // counterpart (in the nonlinear case, not all macroscopic
    // kinematic quantities may be obtained as the volume average of
    // their microstructural counterparts -> a corresponding choice
    // has to be made)

    if (!onlydens) {
      Epetra_SerialDenseMatrix S(3,3);
      S(0,0) = stress(0);
      S(0,1) = stress(3);
      S(0,2) = stress(5);
      S(1,0) = S(0,1);
      S(1,1) = stress(1);
      S(1,2) = stress(4);
      S(2,0) = S(0,2);
      S(2,1) = S(1,2);
      S(2,2) = stress(2);

      Epetra_SerialDenseMatrix P(3,3);
      P.Multiply('N', 'N', 1.0, defgrd, S, 0.);

      P.Scale(integrationfactor);
      P11 += P(0,0);
      P12 += P(0,1);
      P13 += P(0,2);
      P21 += P(1,0);
      P22 += P(1,1);
      P23 += P(1,2);
      P31 += P(2,0);
      P32 += P(2,1);
      P33 += P(2,2);
    }

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  if (!onlydens) {
    double homogP11 = params.get<double>("homogP11", 0.0);
    params.set("homogP11", P11+homogP11);
    double homogP12 = params.get<double>("homogP12", 0.0);
    params.set("homogP12", P12+homogP12);
    double homogP13 = params.get<double>("homogP13", 0.0);
    params.set("homogP13", P13+homogP13);
    double homogP21 = params.get<double>("homogP21", 0.0);
    params.set("homogP21", P21+homogP21);
    double homogP22 = params.get<double>("homogP22", 0.0);
    params.set("homogP22", P22+homogP22);
    double homogP23 = params.get<double>("homogP23", 0.0);
    params.set("homogP23", P23+homogP23);
    double homogP31 = params.get<double>("homogP31", 0.0);
    params.set("homogP31", P31+homogP31);
    double homogP32 = params.get<double>("homogP32", 0.0);
    params.set("homogP32", P32+homogP32);
    double homogP33 = params.get<double>("homogP33", 0.0);
    params.set("homogP33", P33+homogP33);
  }

  return;
}


#endif
#endif
