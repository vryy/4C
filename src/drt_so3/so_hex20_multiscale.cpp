/*!----------------------------------------------------------------------
\file so_hex20_multiscale.cpp
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
#include "so_hex20.H"
#include "../drt_lib/drt_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex20::soh20_homog(ParameterList&  params)
{
  double homogdens = 0.;

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_20 with 20 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  LINALG::Matrix<NUMNOD_SOH20,NUMGPT_SOH20>* shapefct; //[NUMNOD_SOH20][NUMGPT_SOH20]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  LINALG::Matrix<NUMNOD_SOH20*NUMDIM_SOH20,NUMNOD_SOH20>* deriv;    //[NUMGPT_SOH20*NUMDIM][NUMNOD_SOH20]
/* pointer to (static) weight factors at each gp */
  LINALG::Matrix<NUMGPT_SOH20,1>* weights;  //[NUMGPT_SOH20]
  soh20_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH20; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> deriv_gp;
  for (int gp=0; gp<NUMGPT_SOH20; ++gp) {

    // get submatrix of deriv at actual gp
    for (int m=0; m<NUMDIM_SOH20; ++m) {
      for (int n=0; n<NUMNOD_SOH20; ++n) {
        deriv_gp(m,n)=(*deriv)(NUMDIM_SOH20*gp+m,n);
      }
    }

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> jac;
    //jac.Multiply('N','N',1.0,deriv_gp,xrefe,1.0);
    jac.Multiply(deriv_gp, xrefe);

    // compute determinant of Jacobian
    double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // (material) deformation gradient F (=I in reference configuration)
    LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> defgrd(true);
    for (unsigned int i = 0;i<3; ++i) defgrd(i,i) = 1.;

    // Green-Lagrange strains matrix (=0 in reference configuration)
    LINALG::Matrix<NUMSTR_SOH20,1> glstrain(true);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    LINALG::Matrix<NUMSTR_SOH20,NUMSTR_SOH20> cmat(true);
    LINALG::Matrix<NUMSTR_SOH20,1> stress(true);
    double density;
    soh20_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
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
void DRT::ELEMENTS::So_hex20::soh20_read_restart_multi(ParameterList& params)
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
  {

    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_readrestart");
  }
  return;
}

#endif
#endif
