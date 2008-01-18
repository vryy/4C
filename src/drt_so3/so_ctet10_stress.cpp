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
#ifdef D_SOTET
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
void DRT::Elements::So_ctet10::so_ctet10_stress(struct _MATERIAL* material,
                                         vector<double>& disp,
                                         Epetra_SerialDenseMatrix* stresses)
{
 
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
  for (int i=0; i<NUMNOD_SOCTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0] - Nodes()[0]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1] - Nodes()[0]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2] - Nodes()[0]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOCTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOCTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOCTET10+2];
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

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOCTET10,NUMDIM_SOCTET10);
    defgrd.Multiply('T','N',1.0,xcurr,L_aj_int.deriv_gp[gp],0.0);
    
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
    
    for (int i=0; i<NUMNOD_SOCTET10; i++) {
      bop(0,NODDOF_SOCTET10*i+0) = defgrd(0,0)*L_aj_int.deriv_gp[gp](i,0);
      bop(0,NODDOF_SOCTET10*i+1) = defgrd(1,0)*L_aj_int.deriv_gp[gp](i,0);
      bop(0,NODDOF_SOCTET10*i+2) = defgrd(2,0)*L_aj_int.deriv_gp[gp](i,0);
      bop(1,NODDOF_SOCTET10*i+0) = defgrd(0,1)*L_aj_int.deriv_gp[gp](i,1);
      bop(1,NODDOF_SOCTET10*i+1) = defgrd(1,1)*L_aj_int.deriv_gp[gp](i,1);
      bop(1,NODDOF_SOCTET10*i+2) = defgrd(2,1)*L_aj_int.deriv_gp[gp](i,1);
      bop(2,NODDOF_SOCTET10*i+0) = defgrd(0,2)*L_aj_int.deriv_gp[gp](i,2);
      bop(2,NODDOF_SOCTET10*i+1) = defgrd(1,2)*L_aj_int.deriv_gp[gp](i,2);
      bop(2,NODDOF_SOCTET10*i+2) = defgrd(2,2)*L_aj_int.deriv_gp[gp](i,2);
      /* ~~~ */
      bop(3,NODDOF_SOCTET10*i+0) = defgrd(0,0)*L_aj_int.deriv_gp[gp](i,1) +\
                                   defgrd(0,1)*L_aj_int.deriv_gp[gp](i,0);
      bop(3,NODDOF_SOCTET10*i+1) = defgrd(1,0)*L_aj_int.deriv_gp[gp](i,1) +\
                                   defgrd(1,1)*L_aj_int.deriv_gp[gp](i,0);
      bop(3,NODDOF_SOCTET10*i+2) = defgrd(2,0)*L_aj_int.deriv_gp[gp](i,1) +\
                                   defgrd(2,1)*L_aj_int.deriv_gp[gp](i,0);
      bop(4,NODDOF_SOCTET10*i+0) = defgrd(0,1)*L_aj_int.deriv_gp[gp](i,2) +\
                                   defgrd(0,2)*L_aj_int.deriv_gp[gp](i,1);
      bop(4,NODDOF_SOCTET10*i+1) = defgrd(1,1)*L_aj_int.deriv_gp[gp](i,2) +\
                                   defgrd(1,2)*L_aj_int.deriv_gp[gp](i,1);
      bop(4,NODDOF_SOCTET10*i+2) = defgrd(2,1)*L_aj_int.deriv_gp[gp](i,2) +\
                                   defgrd(2,2)*L_aj_int.deriv_gp[gp](i,1);
      bop(5,NODDOF_SOCTET10*i+0) = defgrd(0,2)*L_aj_int.deriv_gp[gp](i,0) +\
                                   defgrd(0,0)*L_aj_int.deriv_gp[gp](i,2);
      bop(5,NODDOF_SOCTET10*i+1) = defgrd(1,2)*L_aj_int.deriv_gp[gp](i,0) +\
                                   defgrd(1,0)*L_aj_int.deriv_gp[gp](i,2);
      bop(5,NODDOF_SOCTET10*i+2) = defgrd(2,2)*L_aj_int.deriv_gp[gp](i,0) +\
                                   defgrd(2,0)*L_aj_int.deriv_gp[gp](i,2);
    }
    

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOCTET10,NUMSTR_SOCTET10);
    Epetra_SerialDenseVector stress(NUMSTR_SOCTET10);
    double density;
    so_ctet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // safe gausspoint stresses
    for (int i=0; i<NUMSTR_SOCTET10; ++i) (*stresses)(gp,i) = stress(i);
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  data_.Add("Stresses",stresses);
} //So_ctet10::soctet10_stress


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET
