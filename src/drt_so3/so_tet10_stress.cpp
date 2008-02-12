/*!----------------------------------------------------------------------*
\file so_tet10_stress.cpp
\brief Evaluate the element stresses

The element stresses are evaluated at the Gauss points.

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_integrator.H"
#include "so_tet10.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
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
 |                                                         vlf 08/07    |
 | vector of material laws                                              |
 | defined in global_control.c											|
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------**####
 |  Do stress calculation (private)                            maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_stress(struct _MATERIAL* material,
                                         vector<double>& disp,
                                         Epetra_SerialDenseMatrix* stresses)
{
  DSTraceHelper dst("So_tet10::sotet10_stress");

/* =============================================================================*
 * CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
 * =============================================================================*/
  const static DRT::ELEMENTS::Integrator_tet10_4point tet10_dis;
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOTET10,NUMDIM_SOTET10);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOTET10,NUMDIM_SOTET10);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NUMDOF_SOTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NUMDOF_SOTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NUMDOF_SOTET10+2];
  }

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    
    /* compute the Jacobian matrix which looks like this:
    **         [  1        1        1  	     1      ]
    **   jac = [ X_,xsi1  X_,xsi2  X_,xsi3  X_,xsi4 ]
    **         [ Y_,xsi1  Y_,xsi2  Y_,xsi3  Y_,xsi4 ]
    **		   [ Z_,xsi1  Z_,xsi2  Z_,xsi3  Z_,xsi4 ]
    */
    
    Epetra_SerialDenseMatrix jac_temp(NUMCOORD_SOTET10-1,NUMCOORD_SOTET10);
    Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
    jac_temp.Multiply('T','N',1.0,xrefe,(tet10_dis.deriv_gp)[gp],0.0);
   
    for (int i=0; i<4; i++) jac(0,i)=1;
    for (int row=0;row<3;row++)
    {
    	for (int col=0;col<4;col++)
    	jac(row+1,col)=jac_temp(row,col);	
    }

    #ifdef VERBOSE_OUTPUT
    cout << "jac\n" << jac;
    cout << "xrefe\n" << xrefe;
	#endif //VERBOSE_OUTPUT

    /* compute partial derivatives at gp xsi_1, xsi_2, xsi_3, xsi_4 material coordinates
    ** by solving   Jac . partials = I_aug   for partials
    ** Inverse of Jacobian is therefore not explicitly computed
    */ 
    /* structure of partials:
    **             [  dxsi_1   dxsi_1   dxsi_1  ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **             [     |        |        |    ]
    ** partials =  [                            ]
    **             [  dxsi_4   dxsi_4   dxsi_4  ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    */   
    
    Epetra_SerialDenseMatrix I_aug(NUMCOORD_SOTET10,NUMDIM_SOTET10);
    Epetra_SerialDenseMatrix partials(NUMCOORD_SOTET10,NUMDIM_SOTET10);
    Epetra_SerialDenseMatrix N_XYZ(NUMNOD_SOTET10,NUMDIM_SOTET10);
    I_aug(1,0)=1;
	I_aug(2,1)=1;
	I_aug(3,2)=1;
	
	#ifdef VERBOSE_OUTPUT
	cout << "jac\n" << jac;
	#endif //VERBOSE_OUTPUT
	
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(partials,I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();        
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) && (err2!=0)){
    	dserror("Inversion of Jacobian failed");
    }

    #ifdef VERBOSE_OUTPUT
    cout << "I_aug\n" << I_aug;
    cout << "partials\n" << partials;
    cout << "deriv_gp\n" << deriv_gp;
    #endif //VERBOSE_OUTPUT

    N_XYZ.Multiply('N','N',1.0,(tet10_dis.deriv_gp)[gp],partials,0.0); //N_XYZ = N_xsi_k*partials
      
    /* structure of N_XYZ:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    N_XYZ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_10    dN_10    dN_10  ]
    **             [  -------  -------  ------- ]
    **             [    dX       dY       dZ    ]
    */
    
    #ifdef VERBOSE_OUTPUT
    cout << "N_XYZ\n" << N_XYZ;
    #endif //VERBOSE_OUTPUT
    //                                      d xcurr 
    // (material) deformation gradient F = --------- = xcurr^T * N_XYZ^T
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

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOTET10,NUMDIM_SOTET10);
    defgrd.Multiply('T','N',1.0,xcurr,N_XYZ,0.0);
    
    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT
	
    // Right Cauchy-Green tensor = F^T * F 
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOTET10,NUMDIM_SOTET10);
   
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);
    
	#ifdef VERBOSE_OUTPUT
	cout << "cauchygreen\n" << cauchygreen;
	getchar();
	#endif //VERBOSE_OUTPUT
	
    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOTET10);
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
    
    Epetra_SerialDenseMatrix bop(NUMSTR_SOTET10,NUMDOF_SOTET10);
    #ifdef VERBOSE_OUTPUT
    cout << bop;
    cout << defgrd;
    cout << N_XYZ;
    #endif //VERBOSE_OUTPUT
    
    for (int i=0; i<NUMNOD_SOTET10; i++) {
      bop(0,NODDOF_SOTET10*i+0) = defgrd(0,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET10*i+1) = defgrd(1,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET10*i+2) = defgrd(2,0)*N_XYZ(i,0);
      bop(1,NODDOF_SOTET10*i+0) = defgrd(0,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET10*i+1) = defgrd(1,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET10*i+2) = defgrd(2,1)*N_XYZ(i,1);
      bop(2,NODDOF_SOTET10*i+0) = defgrd(0,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET10*i+1) = defgrd(1,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET10*i+2) = defgrd(2,2)*N_XYZ(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET10*i+0) = defgrd(0,0)*N_XYZ(i,1) + defgrd(0,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET10*i+1) = defgrd(1,0)*N_XYZ(i,1) + defgrd(1,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET10*i+2) = defgrd(2,0)*N_XYZ(i,1) + defgrd(2,1)*N_XYZ(i,0);
      bop(4,NODDOF_SOTET10*i+0) = defgrd(0,1)*N_XYZ(i,2) + defgrd(0,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET10*i+1) = defgrd(1,1)*N_XYZ(i,2) + defgrd(1,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET10*i+2) = defgrd(2,1)*N_XYZ(i,2) + defgrd(2,2)*N_XYZ(i,1);
      bop(5,NODDOF_SOTET10*i+0) = defgrd(0,2)*N_XYZ(i,0) + defgrd(0,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET10*i+1) = defgrd(1,2)*N_XYZ(i,0) + defgrd(1,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET10*i+2) = defgrd(2,2)*N_XYZ(i,0) + defgrd(2,0)*N_XYZ(i,2);
    }
    

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOTET10,NUMSTR_SOTET10);
    Epetra_SerialDenseVector stress(NUMSTR_SOTET10);
    double density;
    so_tet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // safe gausspoint stresses
    for (int i=0; i<NUMSTR_SOTET10; ++i) (*stresses)(gp,i) = stress(i);
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  data_.Add("Stresses",stresses);
} //So_tet10::sotet10_stress


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET
