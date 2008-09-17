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
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex8::soh8_homog(ParameterList&  params)
{
  double homogdens = 0.;

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMGPT_SOH8>* shapefct; //[NUMNOD_SOH8][NUMGPT_SOH8]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8*NUMDIM_SOH8,NUMNOD_SOH8>* deriv;    //[NUMGPT_SOH8*NUMDIM][NUMNOD_SOH8]
/* pointer to (static) weight factors at each gp */
  LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,1>* weights;  //[NUMGPT_SOH8]
  soh8_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMGPT_SOH8> deriv_gp;
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // get submatrix of deriv at actual gp
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
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> jac;
    //jac.Multiply('N','N',1.0,deriv_gp,xrefe,1.0);
    jac.Multiply(deriv_gp, xrefe);

    // compute determinant of Jacobian
    double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // (material) deformation gradient F (=I in reference configuration)
    Epetra_SerialDenseMatrix defgrd_epetra(NUMDIM_SOH8,NUMDIM_SOH8);
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(defgrd_epetra.A(),true);
    for (unsigned int i = 0;i<3; ++i) defgrd(i,i) = 1.;

    // Green-Lagrange strains matrix (=0 in reference configuration)
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    Epetra_SerialDenseMatrix cmat_epetra(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress_epetra(NUMSTR_SOH8);
    double density;
    soh8_mat_sel(&stress_epetra,&cmat_epetra,&density,&glstrain_epetra,&defgrd_epetra,gp,params);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(cmat_epetra.A(),true);
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_SOH8,1> stress(stress_epetra.A(),true);
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
 |  Set EAS internal variables on the microscale (public)       lw 04/08|
 *----------------------------------------------------------------------*/
// the microscale internal EAS data have to be saved separately for every
// macroscopic Gauss point and set before the determination of microscale
// stiffness etc.

void DRT::ELEMENTS::So_hex8::soh8_set_eas_multi(ParameterList&  params)
{
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldalpha", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldfeas", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldKaainv", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldKda", null);

  if (oldalpha==null || oldfeas==null || oldKaainv==null || oldKda==null)
    dserror("Cannot get EAS internal data from parameter list for multi-scale problems");

  data_.Add("alpha", (*oldalpha)[Id()]);
  data_.Add("feas", (*oldfeas)[Id()]);
  data_.Add("invKaa", (*oldKaainv)[Id()]);
  data_.Add("Kda", (*oldKda)[Id()]);

  return;
}


/*----------------------------------------------------------------------*
 |  Initialize EAS internal variables on the microscale         lw 03/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_eas_init_multi(ParameterList&  params)
{
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > lastalpha =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("lastalpha", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldalpha", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldfeas", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldKaainv", null);
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda =
    params.get<RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > >("oldKda", null);

  (*lastalpha)[Id()] = rcp(new Epetra_SerialDenseMatrix(neas_, 1));
  (*oldalpha)[Id()]  = rcp(new Epetra_SerialDenseMatrix(neas_, 1));
  (*oldfeas)[Id()]   = rcp(new Epetra_SerialDenseMatrix(neas_, 1));
  (*oldKaainv)[Id()] = rcp(new Epetra_SerialDenseMatrix(neas_, neas_));
  (*oldKda)[Id()]    = rcp(new Epetra_SerialDenseMatrix(neas_, NUMDOF_SOH8));

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_read_restart_multi(ParameterList& params)
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {

    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., "multi_readrestart");
  }
  return;
}

#endif
#endif
