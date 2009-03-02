/*!----------------------------------------------------------------------
\file so_hex27_multiscale.cpp
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
#include "so_hex27.H"
#include "../drt_lib/drt_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex27::soh27_homog(ParameterList&  params)
{
  double homogdens = 0.;

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_27 with 27 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  LINALG::Matrix<NUMNOD_SOH27,NUMGPT_SOH27>* shapefct; //[NUMNOD_SOH27][NUMGPT_SOH27]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  LINALG::Matrix<NUMGPT_SOH27*NUMDIM_SOH27,NUMNOD_SOH27>* deriv;    //[NUMGPT_SOH27*NUMDIM][NUMNOD_SOH27]
/* pointer to (static) weight factors at each gp */
  LINALG::Matrix<NUMGPT_SOH27,1>* weights;  //[NUMGPT_SOH27]
  soh27_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH27,NUMDIM_SOH27> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH27; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH27,NUMGPT_SOH27> deriv_gp;
  for (int gp=0; gp<NUMGPT_SOH27; ++gp) {

    // get submatrix of deriv at actual gp
    for (int m=0; m<NUMDIM_SOH27; ++m) {
      for (int n=0; n<NUMGPT_SOH27; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_SOH27*gp+m,n);
      }
    }

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<NUMDIM_SOH27,NUMDIM_SOH27> jac;
    //jac.Multiply('N','N',1.0,deriv_gp,xrefe,1.0);
    jac.Multiply(deriv_gp, xrefe);

    // compute determinant of Jacobian
    double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // (material) deformation gradient F (=I in reference configuration)
    LINALG::Matrix<NUMDIM_SOH27,NUMDIM_SOH27> defgrd(true);
    for (unsigned int i = 0;i<3; ++i) defgrd(i,i) = 1.;

    // Green-Lagrange strains matrix (=0 in reference configuration)
    LINALG::Matrix<NUMSTR_SOH27,1> glstrain(true);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    LINALG::Matrix<NUMSTR_SOH27,NUMSTR_SOH27> cmat(true);
    LINALG::Matrix<NUMSTR_SOH27,1> stress(true);
    double density;
    soh27_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double integrationfactor = detJ * (*weights)(gp);

    homogdens += integrationfactor*density;

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27::soh27_read_restart_multi(ParameterList& params)
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOH27; ++gp)
  {

    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_readrestart");
  }
  return;
}

#endif
#endif
