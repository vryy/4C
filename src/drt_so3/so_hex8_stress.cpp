/*!----------------------------------------------------------------------
\file so_hex8_stress.cpp
\brief Evaluate the element stresses

The element stresses are evaluated at the Gauss points.

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
#include "so_hex8.H"
#include "../discret/drt_discret.H"
#include "../discret/drt_utils.H"
#include "../discret/drt_exporter.H"
#include "../discret/drt_dserror.H"
#include "../discret/linalg_utils.H"
#include "../discret/linalg_serialdensematrix.H"
#include "../discret/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


extern "C" 
{
#include "../headers/standardtypes.h"
// see if we can avoid this #include "../shell8/shell8.h"
}
#include "../discret/dstrc.H"
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  Do stress calculation (private)                            maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_stress(struct _MATERIAL* material, 
                                         vector<double>& disp,
                                         Epetra_SerialDenseMatrix* stresses)
{
  DSTraceHelper dst("So_hex8::soh8_stress");
  
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
    
    xcurr(i,0) = xrefe(i,0) + disp[i*NUMDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NUMDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NUMDOF_SOH8+2];
  }
  
  // Epetra_SerialDenseMatrix stresses(NUMGPT_SOH8,NUMSTR_SOH8);
  
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
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if (err != 0) dserror("Inversion of Jacobian failed");
    
    // (material) deformation gradient F = d xxurr / d xrefe
    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
    defgrd.Multiply('N','N',1.0,N_XYZ,xcurr,1.0);
    
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
    soh8_mat_sel(&stress,&cmat,&density,&glstrain);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc
    
    // safe gausspoint stresses
    for (int i=0; i<NUMSTR_SOH8; ++i) (*stresses)(gp,i) = stress(i);
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
  
  data_.Add("Stresses",stresses);
    
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
