/*!----------------------------------------------------------------------
\brief Everything concerning EAS technology for so_hex8
\level 1
\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "so_sh8p8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                              maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_easinit()
{
  // all parameters are stored in Epetra_SerialDenseMatrix as only
  // those can be added to DRT::Container

  // EAS enhanced strain parameters at currently investigated load/time step
  Epetra_SerialDenseMatrix alpha(neas_, 1);
  // EAS enhanced strain parameters of last converged load/time step
  Epetra_SerialDenseMatrix alphao(neas_, 1);
  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  Epetra_SerialDenseMatrix feas(neas_, 1);
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Epetra_SerialDenseMatrix invKaa(neas_, neas_);
  // EAS matrix K_{d alpha}
  Epetra_SerialDenseMatrix Kda(neas_, NUMDOF_SOH8);
  // EAS increment over last Newton step
  Epetra_SerialDenseMatrix eas_inc(neas_, 1);

  // save EAS data into element container
  data_.Add("alpha", alpha);
  data_.Add("alphao", alphao);
  data_.Add("feas", feas);
  data_.Add("invKaa", invKaa);
  data_.Add("Kda", Kda);
  data_.Add("eas_inc", eas_inc);

  return;
}

/*----------------------------------------------------------------------*
 |  re-initialize EAS data (private)                           maf 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_reiniteas(const DRT::ELEMENTS::So_hex8::EASType EASType)
{
  switch (EASType)
  {
    case DRT::ELEMENTS::So_hex8::soh8_easfull:
      neas_ = 21;
      break;
    case DRT::ELEMENTS::So_hex8::soh8_easmild:
      neas_ = 9;
      break;
    case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
      neas_ = 7;
      break;
    case DRT::ELEMENTS::So_hex8::soh8_easa:
      neas_ = DRT::ELEMENTS::So_sh8p8::NUMEAS_A_;
      break;
    case DRT::ELEMENTS::So_hex8::soh8_easnone:
      neas_ = 0;
      break;
  }
  eastype_ = EASType;
  if (eastype_ == DRT::ELEMENTS::So_hex8::soh8_easnone) return;
  Epetra_SerialDenseMatrix* alpha = NULL;                         // EAS alphas
  Epetra_SerialDenseMatrix* alphao = NULL;                        // EAS alphas
  Epetra_SerialDenseMatrix* feas = NULL;                          // EAS history
  Epetra_SerialDenseMatrix* Kaainv = NULL;                        // EAS history
  Epetra_SerialDenseMatrix* Kda = NULL;                           // EAS history
  Epetra_SerialDenseMatrix* eas_inc = NULL;                       // EAS history
  alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");    // get alpha of previous iteration
  alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // get alpha of previous iteration
  feas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
  Kaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
  Kda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
  eas_inc = data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc");
  if (!alpha || !Kaainv || !Kda || !feas || !eas_inc) dserror("Missing EAS history-data");

  alpha->Reshape(neas_, 1);
  alphao->Reshape(neas_, 1);
  feas->Reshape(neas_, 1);
  Kaainv->Reshape(neas_, neas_);
  Kda->Reshape(neas_, NUMDOF_SOH8);
  eas_inc->Reshape(neas_, 1);

  return;
}

/*----------------------------------------------------------------------*
 |  setup of constant EAS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_eassetup(
    std::vector<Epetra_SerialDenseMatrix>** M_GP,  // M-matrix evaluated at GPs
    double& detJ0,                                 // det of Jacobian at origin
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
        T0invT,                                                   // maps M(origin) local to global
    const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe) const  // material element coords
{
  // vector of df(origin)
  static double df0_vector[NUMDIM_SOH8 * NUMNOD_SOH8] = {-0.125, -0.125, -0.125, +0.125, -0.125,
      -0.125, +0.125, +0.125, -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, +0.125, +0.125, +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> df0(df0_vector);  // copy

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0;
  jac0.Multiply(df0, xrefe);
  // compute determinant of Jacobian at origin
  detJ0 = jac0.Determinant();

  // first, build T0^T transformation matrix which maps the M-matrix
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  T0invT(0, 0) = jac0(0, 0) * jac0(0, 0);
  T0invT(1, 0) = jac0(1, 0) * jac0(1, 0);
  T0invT(2, 0) = jac0(2, 0) * jac0(2, 0);
  T0invT(3, 0) = 2 * jac0(0, 0) * jac0(1, 0);
  T0invT(4, 0) = 2 * jac0(1, 0) * jac0(2, 0);
  T0invT(5, 0) = 2 * jac0(0, 0) * jac0(2, 0);

  T0invT(0, 1) = jac0(0, 1) * jac0(0, 1);
  T0invT(1, 1) = jac0(1, 1) * jac0(1, 1);
  T0invT(2, 1) = jac0(2, 1) * jac0(2, 1);
  T0invT(3, 1) = 2 * jac0(0, 1) * jac0(1, 1);
  T0invT(4, 1) = 2 * jac0(1, 1) * jac0(2, 1);
  T0invT(5, 1) = 2 * jac0(0, 1) * jac0(2, 1);

  T0invT(0, 2) = jac0(0, 2) * jac0(0, 2);
  T0invT(1, 2) = jac0(1, 2) * jac0(1, 2);
  T0invT(2, 2) = jac0(2, 2) * jac0(2, 2);
  T0invT(3, 2) = 2 * jac0(0, 2) * jac0(1, 2);
  T0invT(4, 2) = 2 * jac0(1, 2) * jac0(2, 2);
  T0invT(5, 2) = 2 * jac0(0, 2) * jac0(2, 2);

  T0invT(0, 3) = jac0(0, 0) * jac0(0, 1);
  T0invT(1, 3) = jac0(1, 0) * jac0(1, 1);
  T0invT(2, 3) = jac0(2, 0) * jac0(2, 1);
  T0invT(3, 3) = jac0(0, 0) * jac0(1, 1) + jac0(1, 0) * jac0(0, 1);
  T0invT(4, 3) = jac0(1, 0) * jac0(2, 1) + jac0(2, 0) * jac0(1, 1);
  T0invT(5, 3) = jac0(0, 0) * jac0(2, 1) + jac0(2, 0) * jac0(0, 1);


  T0invT(0, 4) = jac0(0, 1) * jac0(0, 2);
  T0invT(1, 4) = jac0(1, 1) * jac0(1, 2);
  T0invT(2, 4) = jac0(2, 1) * jac0(2, 2);
  T0invT(3, 4) = jac0(0, 1) * jac0(1, 2) + jac0(1, 1) * jac0(0, 2);
  T0invT(4, 4) = jac0(1, 1) * jac0(2, 2) + jac0(2, 1) * jac0(1, 2);
  T0invT(5, 4) = jac0(0, 1) * jac0(2, 2) + jac0(2, 1) * jac0(0, 2);

  T0invT(0, 5) = jac0(0, 0) * jac0(0, 2);
  T0invT(1, 5) = jac0(1, 0) * jac0(1, 2);
  T0invT(2, 5) = jac0(2, 0) * jac0(2, 2);
  T0invT(3, 5) = jac0(0, 0) * jac0(1, 2) + jac0(1, 0) * jac0(0, 2);
  T0invT(4, 5) = jac0(1, 0) * jac0(2, 2) + jac0(2, 0) * jac0(1, 2);
  T0invT(5, 5) = jac0(0, 0) * jac0(2, 2) + jac0(2, 0) * jac0(0, 2);

  // now evaluate T0^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, 1> solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(T0invT);
  int err2 = solve_for_inverseT0.Factor();
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) || (err2 != 0)) dserror("Inversion of T0inv (Jacobian0) failed");

  // build EAS interpolation matrix M, evaluated at the 8 GPs of so_hex8

  // fill up M at each gp
  if (eastype_ == soh8_easmild)
  {
    // static Epetra_SerialDenseMatrix M_mild(MAT::NUM_STRESS_3D*NUMGPT_SOH8,neas_);
    static std::vector<Epetra_SerialDenseMatrix> M_mild(NUMGPT_SOH8);
    static bool M_mild_eval = false;
    /* easmild is the EAS interpolation of 9 modes, based on
    **            r 0 0   0 0 0 0 0 0
    **            0 s 0   0 0 0 0 0 0
    **    M =     0 0 t   0 0 0 0 0 0
    **            0 0 0   r s 0 0 0 0
    **            0 0 0   0 0 s t 0 0
    **            0 0 0   0 0 0 0 r t
    */
    if (!M_mild_eval)
    {  // if true M already evaluated
      // (r,s,t) gp-locations of fully integrated linear 8-node Hex
      const double* r = soh8_get_coordinate_of_gausspoints(0);
      const double* s = soh8_get_coordinate_of_gausspoints(1);
      const double* t = soh8_get_coordinate_of_gausspoints(2);

      // fill up M at each gp
      for (unsigned i = 0; i < NUMGPT_SOH8; ++i)
      {
        M_mild[i].Shape(MAT::NUM_STRESS_3D, neas_);
        M_mild[i](0, 0) = r[i];
        M_mild[i](1, 1) = s[i];
        M_mild[i](2, 2) = t[i];

        M_mild[i](3, 3) = r[i];
        M_mild[i](3, 4) = s[i];
        M_mild[i](4, 5) = s[i];
        M_mild[i](4, 6) = t[i];
        M_mild[i](5, 7) = r[i];
        M_mild[i](5, 8) = t[i];
      }
      M_mild_eval = true;  // now the array is filled statically
    }

    // return adress of just evaluated matrix
    *M_GP = &M_mild;  // return adress of static object to target of pointer
  }
  else if (eastype_ == soh8_easfull)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_full(NUMGPT_SOH8);
    static bool M_full_eval = false;
    /* easfull is the EAS interpolation of 21 modes, based on
    **            r 0 0   0 0 0 0 0 0   0  0  0  0  0  0   rs rt 0  0  0  0
    **            0 s 0   0 0 0 0 0 0   0  0  0  0  0  0   0  0  rs st 0  0
    **    M =     0 0 t   0 0 0 0 0 0   0  0  0  0  0  0   0  0  0  0  rt st
    **            0 0 0   r s 0 0 0 0   rt st 0  0  0  0   0  0  0  0  0  0
    **            0 0 0   0 0 s t 0 0   0  0  rs rt 0  0   0  0  0  0  0  0
    **            0 0 0   0 0 0 0 r t   0  0  0  0  rs st  0  0  0  0  0  0
    */
    if (!M_full_eval)
    {  // if true M already evaluated
      // (r,s,t) gp-locations of fully integrated linear 8-node Hex
      const double* r = soh8_get_coordinate_of_gausspoints(0);
      const double* s = soh8_get_coordinate_of_gausspoints(1);
      const double* t = soh8_get_coordinate_of_gausspoints(2);

      // fill up M at each gp
      for (unsigned i = 0; i < NUMGPT_SOH8; ++i)
      {
        M_full[i].Shape(MAT::NUM_STRESS_3D, neas_);
        M_full[i](0, 0) = r[i];
        M_full[i](0, 15) = r[i] * s[i];
        M_full[i](0, 16) = r[i] * t[i];
        M_full[i](1, 1) = s[i];
        M_full[i](1, 17) = r[i] * s[i];
        M_full[i](1, 18) = s[i] * t[i];
        M_full[i](2, 2) = t[i];
        M_full[i](2, 19) = r[i] * t[i];
        M_full[i](2, 20) = s[i] * t[i];

        M_full[i](3, 3) = r[i];
        M_full[i](3, 4) = s[i];
        M_full[i](3, 9) = r[i] * t[i];
        M_full[i](3, 10) = s[i] * t[i];
        M_full[i](4, 5) = s[i];
        M_full[i](4, 6) = t[i];
        M_full[i](4, 11) = r[i] * s[i];
        M_full[i](4, 12) = r[i] * t[i];
        M_full[i](5, 7) = r[i];
        M_full[i](5, 8) = t[i];
        M_full[i](5, 13) = r[i] * s[i];
        M_full[i](5, 14) = s[i] * t[i];
      }
      M_full_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    *M_GP = &M_full;  // return adress of static object to target of pointer
  }
  else if (eastype_ == soh8_eassosh8)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_sosh8(NUMGPT_SOH8);
    static bool M_sosh8_eval = false;
    /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
    ** consisting of 7 modes, based on
    **            r 0 0   0 0 0  0
    **            0 s 0   0 0 0  0
    **    M =     0 0 t   0 0 rt st
    **            0 0 0   r s 0  0
    **            0 0 0   0 0 0  0
    **            0 0 0   0 0 0  0
    */
    if (!M_sosh8_eval)
    {  // if true M already evaluated
      // (r,s,t) gp-locations of fully integrated linear 8-node Hex
      const double* r = soh8_get_coordinate_of_gausspoints(0);
      const double* s = soh8_get_coordinate_of_gausspoints(1);
      const double* t = soh8_get_coordinate_of_gausspoints(2);

      // fill up M at each gp
      for (unsigned i = 0; i < NUMGPT_SOH8; ++i)
      {
        M_sosh8[i].Shape(MAT::NUM_STRESS_3D, neas_);
        M_sosh8[i](0, 0) = r[i];
        M_sosh8[i](1, 1) = s[i];
        M_sosh8[i](2, 2) = t[i];
        M_sosh8[i](2, 5) = r[i] * t[i];
        M_sosh8[i](2, 6) = s[i] * t[i];

        M_sosh8[i](3, 3) = r[i];
        M_sosh8[i](3, 4) = s[i];
      }
      M_sosh8_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    *M_GP = &M_sosh8;  // return adress of static object to target of pointer
  }
  else if (eastype_ == soh8_easa)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_sosh8(NUMGPT_SOH8);
    static bool M_sosh8_eval = false;
    /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
    ** consisting of 7 modes, based on
    **            r 0 0   0 0 0  0             // E_rr
    **            0 s 0   0 0 0  0             // E_ss
    **    M =     0 0 t   0 0 rt st            // E_tt
    **            0 0 0   r s 0  0             // E_rs
    **            0 0 0   0 0 0  0             // E_st
    **            0 0 0   0 0 0  0             // E_tr
    */
    if (!M_sosh8_eval)
    {  // if true M already evaluated
      // (r,s,t) gp-locations of fully integrated linear 8-node Hex
      const double* t = soh8_get_coordinate_of_gausspoints(2);

      // fill up M at each gp
      for (unsigned i = 0; i < NUMGPT_SOH8; ++i)
      {
        M_sosh8[i].Shape(MAT::NUM_STRESS_3D, neas_);
        int e = 0;
        M_sosh8[i](2, e++) = t[i] * t[i] * t[i];

        // sosh8
        // M_sosh8[i](0,e++) = r[i];
        // M_sosh8[i](1,e++) = s[i];
        // M_sosh8[i](2,e++) = t[i]; M_sosh8[i](2,e++) = r[i]*t[i]; M_sosh8[i](2,e++) = s[i]*t[i];
        // M_sosh8[i](3,e++) = r[i]; M_sosh8[i](3,e++) = s[i];

        // mild
        // M_sosh8[i](0,e++) = r[i];
        // M_sosh8[i](1,e++) = s[i];
        // M_sosh8[i](2,e++) = t[i];
        // M_sosh8[i](3,e++) = r[i]; M_sosh8[i](3,e++) = s[i];
        // M_sosh8[i](4,e++) = s[i]; M_sosh8[i](4,e++) = t[i];
        // M_sosh8[i](5,e++) = r[i]; M_sosh8[i](5,e++) = t[i];

        if (e != neas_) dserror("Too many/few EAS shape functions: %d != %d", e, neas_);
      }
      M_sosh8_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    *M_GP = &M_sosh8;  // return adress of static object to target of pointer
  }
  else
  {
    dserror("eastype not implemented");
  }
}  // end of soh8_eassetup
