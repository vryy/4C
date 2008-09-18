/*!----------------------------------------------------------------------*
\file so_ctet10_stress.cpp
\brief Evaluate the element stresses

The element stresses are evaluated at the Gauss points.

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#if 0

#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_ctet10.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         vlf 08/07    |
 | vector of material laws                                              |
 | defined in global_control.c											|
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

#if 0
/*----------------------------------------------------------------------*
 |  Do stress calculation (private)                            vlf 05/08|
 |  Consistent computation of stresses, not used because too expensive  |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::so_ctet10_nlnstiffstress(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      struct _MATERIAL*         material,       // element material data
      const bool                cauchy)         // stress output options
{

 if (elestress != NULL){
/* =========================*
 * SUB STRUCTURE for CTET_10
 * =========================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOCTET10+1,NUMDIM_SOCTET10);  // material coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_10  Y_10  Z_10 ]
    */
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOCTET10,NUMDIM_SOCTET10);  // current  coord. of element
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_10  y_10  z_10 ]
    */

  Epetra_SerialDenseMatrix xdisp(NUMNOD_SOCTET10,NUMDIM_SOCTET10);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOCTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0] - Nodes()[0]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1] - Nodes()[0]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2] - Nodes()[0]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOCTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOCTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOCTET10+2];

    xdisp(i,0) = disp[i*NODDOF_SOCTET10+0];
    xdisp(i,1) = disp[i*NODDOF_SOCTET10+1];
    xdisp(i,2) = disp[i*NODDOF_SOCTET10+2];
  }

  //create the midpoint of the tetrahedron as the 11th node of the element
  xrefe(NUMNOD_SOCTET10,0)=ctet10_midpoint(0);
  xrefe(NUMNOD_SOCTET10,1)=ctet10_midpoint(1);
  xrefe(NUMNOD_SOCTET10,2)=ctet10_midpoint(2);

  //build the sub-elements for the whole tetrahedron
  SUB_STRUCTURE sub(xrefe);

  //remove the extra added 11th node
  xrefe.Reshape(NUMNOD_SOCTET10,NUMDIM_SOCTET10);

  //create the integrator from the sub-elements
  L_AJ_integrator L_aj_int(sub);
  L_AJ_stress_integrator stress_int(sub);
  Tet10c_integrator_5point pt5_integr(sub);
  const Integrator_tet10_10node  node10_integr;

  vector<Epetra_SerialDenseMatrix> stress_vec;
  stress_vec.resize(NUMNOD_SOCTET10);
  for (int node = 0; node <NUMNOD_SOCTET10; node++)
  {stress_vec[node].Reshape(3,3);}
  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOCTET10; gp++) {
    /* structure of L_aj_int.deriv_gp, which replaces N_XYZ:
    **                               [   dN_1     dN_1     dN_1   ]
    **                               [  ------   ------   ------  ]
    **                               [    dX       dY       dZ    ]
    **    L_XYZ (is equivalent to)=  [     |        |        |    ]
    **                               [                            ]
    **                               [   dN_10    dN_10    dN_10  ]
    **                               [  -------  -------  ------- ]
    **                               [    dX       dY       dZ    ]
    */

    //                                      d xcurr
    // (material) deformation gradient F = --------- = L_XYZ^T * xcurr
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    /* Here additionally a decompostion of F = L_XYZ^T * (xrefe + (x)disp)
     * Knowing that L_XYZ^T * xrefe = diag(1,1,1) ,
     * we can ensure F = diag(1,1,1) if disp = 0
     */

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOCTET10,NUMDIM_SOCTET10);

    defgrd.Multiply('T','N',1.0,xdisp,L_aj_int.deriv_gp[gp],0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;


    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT

    // Right Cauchy-Green tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOCTET10,NUMDIM_SOCTET10);

    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);

	#ifdef VERBOSE_OUTPUT
	cout << "cauchygreen\n" << cauchygreen;
	getchar();
	#endif //VERBOSE_OUTPUT

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOCTET10);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    #ifdef VERBOSE_OUTPUT
	cout << "glstrain\n" << glstrain;
	#endif //VERBOSE_OUTPUT

	/*----------------------------------------------------------------------*
      the B-operator used is equivalent to the one used in hex8, this needs
      to be checked if it is ok, but from the mathematics point of view, the only
      thing that needed to be changed is tho NUMDOF       vlf 07/07
      ----------------------------------------------------------------------*/
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

    Epetra_SerialDenseMatrix bop(NUMSTR_SOCTET10,NUMDOF_SOCTET10);
    #ifdef VERBOSE_OUTPUT
    cout << bop;
    cout << defgrd;
    cout << N_XYZ;
    #endif //VERBOSE_OUTPUT

    for (int numnode=0; numnode<NUMNOD_SOTET10; numnode++) {
    	for (int numdof=0; numdof<NODDOF_SOTET10; numdof++) {
      	bop(0,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*L_aj_int.deriv_gp[gp](numnode,0);
      	bop(1,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*L_aj_int.deriv_gp[gp](numnode,1);
      	bop(2,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*L_aj_int.deriv_gp[gp](numnode,2);
      	bop(3,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*L_aj_int.deriv_gp[gp](numnode,1) + \
      			    						   defgrd(numdof,1)*L_aj_int.deriv_gp[gp](numnode,0);
      	bop(4,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*L_aj_int.deriv_gp[gp](numnode,2) + \
      										   defgrd(numdof,2)*L_aj_int.deriv_gp[gp](numnode,1);
      	bop(5,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*L_aj_int.deriv_gp[gp](numnode,0) + \
      										   defgrd(numdof,0)*L_aj_int.deriv_gp[gp](numnode,2);
    	}
    }

  	#ifdef VERBOSE_OUTPUT
	cout << "bop\n" << bop;
	#endif //VERBOSE_OUTPUT

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */

    Epetra_SerialDenseMatrix cmat(NUMSTR_SOCTET10,NUMSTR_SOCTET10);
    Epetra_SerialDenseVector stress(NUMSTR_SOCTET10);
    double density;
    so_ctet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);

  	Epetra_SerialDenseVector temp_vec(NUMCOORD_SOCTET10);
  	temp_vec.Multiply('N','N',1.0,L_aj_int.Mbc_inv,pt5_integr.shapefct_gp[gp],0.0);

    Epetra_SerialDenseMatrix PK2_stress(3,3);
    PK2_stress(0,0)=stress(0);
    PK2_stress(1,1)=stress(1);
    PK2_stress(2,2)=stress(2);
    PK2_stress(0,1)= PK2_stress(1,0)=stress(3);
	PK2_stress(1,2)= PK2_stress(2,1)=stress(4);
	PK2_stress(0,2)= PK2_stress(2,0)=stress(5);

	Epetra_SerialDenseMatrix PK1_stress(3,3);

	PK1_stress.Multiply('N','N',1.0,defgrd,PK2_stress,0.0);


    for (int node = 0; node <NUMNOD_SOCTET10; node++)
    {
       double temp_weight = temp_vec.Dot(node10_integr.shapefct_gp_lin[node])*L_aj_int.weights(gp);
       Epetra_SerialDenseMatrix temp_stress(PK1_stress);
       temp_stress.Scale(temp_weight);
       stress_vec[node]+=temp_stress;
      	//cout << "stress Vec " << node  << " weight " << temp_weight << PK2_stress;

    }

 	#ifdef VERBOSE_OUTPUT
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

	if (force != NULL && stiffmatrix != NULL){
	#ifdef VERBOSE_OUTPUT
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
	/**UNCOMMENT NEXT LINE FOR NONLINEAR KINEMATICS	**/
    }

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/


  for (int gp=0; gp<NUMNOD_SOCTET10; gp++) {
    /* structure of L_aj_int.deriv_gp, which replaces N_XYZ:
    **                               [   dN_1     dN_1     dN_1   ]
    **                               [  ------   ------   ------  ]
    **                               [    dX       dY       dZ    ]
    **    L_XYZ (is equivalent to)=  [     |        |        |    ]
    **                               [                            ]
    **                               [   dN_10    dN_10    dN_10  ]
    **                               [  -------  -------  ------- ]
    **                               [    dX       dY       dZ    ]
    */

    //                                      d xcurr
    // (material) deformation gradient F = --------- = L_XYZ^T * xcurr
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    /* Here additionally a decompostion of F = L_XYZ^T * (xrefe + (x)disp)
     * Knowing that L_XYZ^T * xrefe = diag(1,1,1) ,
     * we can ensure F = diag(1,1,1) if disp = 0
     */

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOCTET10,NUMDIM_SOCTET10);

    defgrd.Multiply('T','N',1.0,xdisp,stress_int.deriv_gp[gp],0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;

    //cout << "solving matrix"; getchar();

    Epetra_SerialDenseMatrix PK2_Stress(3,3);
    Epetra_SerialDenseMatrix PK1_Stress(stress_vec[gp]);
    Epetra_SerialDenseSolver solve_for_PK;  // solve A.X=B
    solve_for_PK.SetMatrix(defgrd);            // set A=Mbc)
    /*cout << defgrd;
    cout << PK2_Stress;
    cout << PK1_Stress;
    cout << stress_vec[gp];*/
    solve_for_PK.SetVectors(PK2_Stress,PK1_Stress);// set X=Mbc_inv, B=I_matrix
    //cout << "solving matrix"; getchar();
    solve_for_PK.FactorWithEquilibration(true);
    solve_for_PK.Factor();
    solve_for_PK.Solve();         // Mbc_inv = Mbc^-1.I_matrix
    //if ((err != 0) && (err2!=0)){
    //	dserror("Inversion of Mbc");
    //}
	//cout << "final PK2 stress" ; //<<  PK2_Stress;
	//cout << "matrix solved"; getchar();
    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT

   if (elestress != NULL){
        (*elestress)(gp,0) = PK2_Stress(0,0);
        (*elestress)(gp,1) = PK2_Stress(1,1);
        (*elestress)(gp,2) = PK2_Stress(2,2);
        (*elestress)(gp,3) = PK2_Stress(0,1);
        (*elestress)(gp,4) = PK2_Stress(1,2);
        (*elestress)(gp,5) = PK2_Stress(2,0);
    }

    //    cout << (*elestress);
   /* =========================================================================*/
  }/* ==================================================== end of Loop over NODES */
  }
  return;
} // DRT::ELEMENTS::So_ctet10::soctet10_nlnstress
#endif

#if 1
/*----------------------------------------------------------------------*
 | Do stress calculation (private)                            vlf 05/08 |
 | Simplified evaluation of stresses, by computing them from the local  |
 | assumed deformation gradient                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::so_ctet10_nlnstiffstress(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      struct _MATERIAL*         material,       // element material data
      const bool                cauchy)         // stress output options
{

 if (elestress != NULL){
/* =========================*
 * SUB STRUCTURE for CTET_10
 * =========================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOCTET10+1,NUMDIM_SOCTET10);  // material coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_10  Y_10  Z_10 ]
    */
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOCTET10,NUMDIM_SOCTET10);  // current  coord. of element
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_10  y_10  z_10 ]
    */

  Epetra_SerialDenseMatrix xdisp(NUMNOD_SOCTET10,NUMDIM_SOCTET10);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOCTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0] - Nodes()[0]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1] - Nodes()[0]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2] - Nodes()[0]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOCTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOCTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOCTET10+2];

    xdisp(i,0) = disp[i*NODDOF_SOCTET10+0];
    xdisp(i,1) = disp[i*NODDOF_SOCTET10+1];
    xdisp(i,2) = disp[i*NODDOF_SOCTET10+2];
  }

  //create the midpoint of the tetrahedron as the 11th node of the element
  xrefe(NUMNOD_SOCTET10,0)=ctet10_midpoint(0);
  xrefe(NUMNOD_SOCTET10,1)=ctet10_midpoint(1);
  xrefe(NUMNOD_SOCTET10,2)=ctet10_midpoint(2);

  //build the sub-elements for the whole tetrahedron
  SUB_STRUCTURE sub(xrefe);

  //remove the extra added 11th node
  xrefe.Reshape(NUMNOD_SOCTET10,NUMDIM_SOCTET10);

  //create the integrator from the sub-elements
  L_AJ_stress_integrator stress_int(sub);

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMNOD_SOCTET10; gp++) {
    /* structure of L_aj_int.deriv_gp, which replaces N_XYZ:
    **                               [   dN_1     dN_1     dN_1   ]
    **                               [  ------   ------   ------  ]
    **                               [    dX       dY       dZ    ]
    **    L_XYZ (is equivalent to)=  [     |        |        |    ]
    **                               [                            ]
    **                               [   dN_10    dN_10    dN_10  ]
    **                               [  -------  -------  ------- ]
    **                               [    dX       dY       dZ    ]
    */

    //                                      d xcurr
    // (material) deformation gradient F = --------- = L_XYZ^T * xcurr
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    /* Here additionally a decompostion of F = L_XYZ^T * (xrefe + (x)disp)
     * Knowing that L_XYZ^T * xrefe = diag(1,1,1) ,
     * we can ensure F = diag(1,1,1) if disp = 0
     */

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOCTET10,NUMDIM_SOCTET10);

    defgrd.Multiply('T','N',1.0,xdisp,stress_int.deriv_gp[gp],0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;


    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT

    // Right Cauchy-Green tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOCTET10,NUMDIM_SOCTET10);

    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);

	#ifdef VERBOSE_OUTPUT
	cout << "cauchygreen\n" << cauchygreen;
	getchar();
	#endif //VERBOSE_OUTPUT

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOCTET10);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    #ifdef VERBOSE_OUTPUT
	cout << "glstrain\n" << glstrain;
	#endif //VERBOSE_OUTPUT

	/*----------------------------------------------------------------------*
      the B-operator used is equivalent to the one used in hex8, this needs
      to be checked if it is ok, but from the mathematics point of view, the only
      thing that needed to be changed is tho NUMDOF       vlf 07/07
      ----------------------------------------------------------------------*/
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

    Epetra_SerialDenseMatrix bop(NUMSTR_SOCTET10,NUMDOF_SOCTET10);
    #ifdef VERBOSE_OUTPUT
    cout << bop;
    cout << defgrd;
    cout << N_XYZ;
    #endif //VERBOSE_OUTPUT

    for (int numnode=0; numnode<NUMNOD_SOTET10; numnode++) {
    	for (int numdof=0; numdof<NODDOF_SOTET10; numdof++) {
      	bop(0,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*stress_int.deriv_gp[gp](numnode,0);
      	bop(1,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*stress_int.deriv_gp[gp](numnode,1);
      	bop(2,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*stress_int.deriv_gp[gp](numnode,2);
      	bop(3,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*stress_int.deriv_gp[gp](numnode,1) + \
      			    						   defgrd(numdof,1)*stress_int.deriv_gp[gp](numnode,0);
      	bop(4,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*stress_int.deriv_gp[gp](numnode,2) + \
      										   defgrd(numdof,2)*stress_int.deriv_gp[gp](numnode,1);
      	bop(5,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*stress_int.deriv_gp[gp](numnode,0) + \
      										   defgrd(numdof,0)*stress_int.deriv_gp[gp](numnode,2);
    	}
    }

  	#ifdef VERBOSE_OUTPUT
	cout << "bop\n" << bop;
	#endif //VERBOSE_OUTPUT

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */

    Epetra_SerialDenseMatrix cmat(NUMSTR_SOCTET10,NUMSTR_SOCTET10);
    Epetra_SerialDenseVector stress(NUMSTR_SOCTET10);
    double density;
    so_ctet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);

     if (elestress != NULL){
      for (int i = 0; i < NUMSTR_SOTET10; ++i) {
        (*elestress)(gp,i) = stress(i);
      }
    }

 	#ifdef VERBOSE_OUTPUT
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
     //    cout << (*elestress);
 }
  return;
} // DRT::ELEMENTS::So_ctet10::soctet10_nlnstress
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

#endif
