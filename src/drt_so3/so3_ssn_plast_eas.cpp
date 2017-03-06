/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_eas.cpp
\brief Everything concerning EAS technology for so3_ssn_plast
       Mainly copied from so_hex8_eas.cpp. Redundancy needed,
       because of hard coded Gauss point in the so_hex8
       which do not coincide with the intrepid Gauss points
       used in the so3_ssn_plast.
\level 2
\maintainer Alexander Seitz
*/

/*----------------------------------------------------------------------*
 | headers                                                  seitz 04/14 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                           seitz 04/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasInit()
{
  if (eastype_!=soh8p_easnone && eastype_!=soh8p_easmild &&
      eastype_!=soh8p_easfull && eastype_!=soh8p_eassosh8)
    dserror("unknown EAS type for so3_ssn_plast");
  else
    neas_=eastype_;

  if (eastype_!=soh8p_easnone)
  {
    KaaInv_                             = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,neas_,true));
    Kad_                                = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,numdofperelement_,true));
    feas_                               = Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
    alpha_eas_                          = Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
    alpha_eas_last_timestep_            = Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
    alpha_eas_delta_over_last_timestep_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
    alpha_eas_inc_                      = Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
    Kba_                                = Teuchos::rcp(new std::vector<LINALG::SerialDenseMatrix>(numgpt_,LINALG::SerialDenseMatrix(5,neas_,true)));

    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    if (probtype==prb_tsi)
    {
      KaT_=Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,nen_,true));
      KdT_eas_=Teuchos::rcp(new LINALG::Matrix<numdofperelement_,nen_>);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  setup EAS data (private)                                seitz 04/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasSetup(
    std::vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
    double& detJ0,                      // det of Jacobian at origin
    LINALG::Matrix<numstr_,numstr_>& T0invT,   // maps M(origin) local to global
    const LINALG::Matrix<nen_,nsd_>& xrefe)    // material element coords
{
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<nsd_,nen_> df0;
  DRT::UTILS::shape_function_3D_deriv1(df0, 0.0, 0.0, 0.0, DRT::Element::hex8);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> jac0;
  jac0.Multiply(df0,xrefe);
  // compute determinant of Jacobian at origin
  detJ0 = jac0.Determinant();

  // first, build T0^T transformation matrix which maps the M-matrix
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  T0invT(0,0) = jac0(0,0) * jac0(0,0);
  T0invT(1,0) = jac0(1,0) * jac0(1,0);
  T0invT(2,0) = jac0(2,0) * jac0(2,0);
  T0invT(3,0) = 2 * jac0(0,0) * jac0(1,0);
  T0invT(4,0) = 2 * jac0(1,0) * jac0(2,0);
  T0invT(5,0) = 2 * jac0(0,0) * jac0(2,0);

  T0invT(0,1) = jac0(0,1) * jac0(0,1);
  T0invT(1,1) = jac0(1,1) * jac0(1,1);
  T0invT(2,1) = jac0(2,1) * jac0(2,1);
  T0invT(3,1) = 2 * jac0(0,1) * jac0(1,1);
  T0invT(4,1) = 2 * jac0(1,1) * jac0(2,1);
  T0invT(5,1) = 2 * jac0(0,1) * jac0(2,1);

  T0invT(0,2) = jac0(0,2) * jac0(0,2);
  T0invT(1,2) = jac0(1,2) * jac0(1,2);
  T0invT(2,2) = jac0(2,2) * jac0(2,2);
  T0invT(3,2) = 2 * jac0(0,2) * jac0(1,2);
  T0invT(4,2) = 2 * jac0(1,2) * jac0(2,2);
  T0invT(5,2) = 2 * jac0(0,2) * jac0(2,2);

  T0invT(0,3) = jac0(0,0) * jac0(0,1);
  T0invT(1,3) = jac0(1,0) * jac0(1,1);
  T0invT(2,3) = jac0(2,0) * jac0(2,1);
  T0invT(3,3) = jac0(0,0) * jac0(1,1) + jac0(1,0) * jac0(0,1);
  T0invT(4,3) = jac0(1,0) * jac0(2,1) + jac0(2,0) * jac0(1,1);
  T0invT(5,3) = jac0(0,0) * jac0(2,1) + jac0(2,0) * jac0(0,1);


  T0invT(0,4) = jac0(0,1) * jac0(0,2);
  T0invT(1,4) = jac0(1,1) * jac0(1,2);
  T0invT(2,4) = jac0(2,1) * jac0(2,2);
  T0invT(3,4) = jac0(0,1) * jac0(1,2) + jac0(1,1) * jac0(0,2);
  T0invT(4,4) = jac0(1,1) * jac0(2,2) + jac0(2,1) * jac0(1,2);
  T0invT(5,4) = jac0(0,1) * jac0(2,2) + jac0(2,1) * jac0(0,2);

  T0invT(0,5) = jac0(0,0) * jac0(0,2);
  T0invT(1,5) = jac0(1,0) * jac0(1,2);
  T0invT(2,5) = jac0(2,0) * jac0(2,2);
  T0invT(3,5) = jac0(0,0) * jac0(1,2) + jac0(1,0) * jac0(0,2);
  T0invT(4,5) = jac0(1,0) * jac0(2,2) + jac0(2,0) * jac0(1,2);
  T0invT(5,5) = jac0(0,0) * jac0(2,2) + jac0(2,0) * jac0(0,2);

  // now evaluate T0^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<numstr_,numstr_,1> solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(T0invT);
  int err2 = solve_for_inverseT0.Factor();
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) || (err2!=0)) dserror("Inversion of T0inv (Jacobian0) failed");

  // build EAS interpolation matrix M, evaluated at the 8 GPs of so_hex8

  // fill up M at each gp
  if (eastype_ == soh8p_easmild)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_mild(numgpt_);
    static bool M_mild_eval;
    /* easmild is the EAS interpolation of 9 modes, based on
     **            r 0 0   0 0 0 0 0 0
     **            0 s 0   0 0 0 0 0 0
     **    M =     0 0 t   0 0 0 0 0 0
     **            0 0 0   r s 0 0 0 0
     **            0 0 0   0 0 s t 0 0
     **            0 0 0   0 0 0 0 r t
     **
     ** (r,s,t) gp-locations of fully integrated linear 8-node Hex
     */
    if (!M_mild_eval)// if true M already evaluated
    {
      // fill up M at each gp
      for (int i=0; i<numgpt_; ++i) {
        M_mild[i].Shape(numstr_,neas_);
        M_mild[i](0,0) = xsi_.at(i)(0);
        M_mild[i](1,1) = xsi_.at(i)(1);
        M_mild[i](2,2) = xsi_.at(i)(2);

        M_mild[i](3,3) = xsi_.at(i)(0); M_mild[i](3,4) = xsi_.at(i)(1);
        M_mild[i](4,5) = xsi_.at(i)(1); M_mild[i](4,6) = xsi_.at(i)(2);
        M_mild[i](5,7) = xsi_.at(i)(0); M_mild[i](5,8) = xsi_.at(i)(2);
      }
      M_mild_eval = true;  // now the array is filled statically
    }

    // return adress of just evaluated matrix
    *M_GP = &M_mild;       // return adress of static object to target of pointer
  }
  else if (eastype_==soh8p_easfull)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_full(numgpt_);
    static bool M_full_eval;
    /* easfull is the EAS interpolation of 21 modes, based on
    **            r 0 0   0 0 0 0 0 0   0  0  0  0  0  0   rs rt 0  0  0  0
    **            0 s 0   0 0 0 0 0 0   0  0  0  0  0  0   0  0  rs st 0  0
    **    M =     0 0 t   0 0 0 0 0 0   0  0  0  0  0  0   0  0  0  0  rt st
    **            0 0 0   r s 0 0 0 0   rt st 0  0  0  0   0  0  0  0  0  0
    **            0 0 0   0 0 s t 0 0   0  0  rs rt 0  0   0  0  0  0  0  0
    **            0 0 0   0 0 0 0 r t   0  0  0  0  rs st  0  0  0  0  0  0
    **
    **             (r,s,t) gp-locations of fully integrated linear 8-node Hex
    */
    if (!M_full_eval)// if true M already evaluated
    {
      // fill up M at each gp
      for (int i=0; i<numgpt_; ++i) {
        M_full[i].Shape(numstr_,neas_);
        M_full[i](0,0) = xsi_.at(i)(0);        M_full[i](0,15) = xsi_.at(i)(0)*xsi_.at(i)(1); M_full[i](0,16) = xsi_.at(i)(0)*xsi_.at(i)(2);
        M_full[i](1,1) = xsi_.at(i)(1);        M_full[i](1,17) = xsi_.at(i)(0)*xsi_.at(i)(1); M_full[i](1,18) = xsi_.at(i)(1)*xsi_.at(i)(2);
        M_full[i](2,2) = xsi_.at(i)(2);        M_full[i](2,19) = xsi_.at(i)(0)*xsi_.at(i)(2); M_full[i](2,20) = xsi_.at(i)(1)*xsi_.at(i)(2);

        M_full[i](3,3) = xsi_.at(i)(0); M_full[i](3,4) = xsi_.at(i)(1);   M_full[i](3, 9) = xsi_.at(i)(0)*xsi_.at(i)(2); M_full[i](3,10) = xsi_.at(i)(1)*xsi_.at(i)(2);
        M_full[i](4,5) = xsi_.at(i)(1); M_full[i](4,6) = xsi_.at(i)(2);   M_full[i](4,11) = xsi_.at(i)(0)*xsi_.at(i)(1); M_full[i](4,12) = xsi_.at(i)(0)*xsi_.at(i)(2);
        M_full[i](5,7) = xsi_.at(i)(0); M_full[i](5,8) = xsi_.at(i)(2);   M_full[i](5,13) = xsi_.at(i)(0)*xsi_.at(i)(1); M_full[i](5,14) = xsi_.at(i)(1)*xsi_.at(i)(2);
      }
      M_full_eval = true;  // now the array is filled statically
    }
      // return adress of just evaluated matrix
      *M_GP = &M_full;            // return adress of static object to target of pointer
  }
  else if (eastype_==soh8p_eassosh8)
  {
    static std::vector<Epetra_SerialDenseMatrix> M_sosh8(numgpt_);
    static bool M_sosh8_eval;
    /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
     ** consisting of 7 modes, based on
     **            r 0 0   0 0 0  0
     **            0 s 0   0 0 0  0
     **    M =     0 0 t   0 0 rt st
     **            0 0 0   r s 0  0
     **            0 0 0   0 0 0  0
     **            0 0 0   0 0 0  0
     */
    if (!M_sosh8_eval) // if true M already evaluated
    {
      // fill up M at each gp
      for (int i=0; i<numgpt_; ++i) {
      M_sosh8[i].Shape(numstr_,neas_);
      M_sosh8[i](0,0) = xsi_.at(i)(0);
      M_sosh8[i](1,1) = xsi_.at(i)(1);
      M_sosh8[i](2,2) = xsi_.at(i)(2); M_sosh8[i](2,5) = xsi_.at(i)(0)*xsi_.at(i)(2); M_sosh8[i](2,6) = xsi_.at(i)(1)*xsi_.at(i)(2);

      M_sosh8[i](3,3) = xsi_.at(i)(0); M_sosh8[i](3,4) = xsi_.at(i)(1);
    }

      M_sosh8_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    *M_GP = &M_sosh8;            // return adress of static object to target of pointer

  }
  else
    dserror("this EAS type not yet implemented");


  return;
}

/*----------------------------------------------------------------------*
 |  Defgrd consistent with enhanced GL strain (private)     seitz 04/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::CalcConsistentDefgrd(LINALG::Matrix<3,3> defgrd_disp,
    LINALG::Matrix<6,1> glstrain_mod,
    LINALG::Matrix<3,3>& defgrd_mod)
{
  LINALG::Matrix<3,3> R;      // rotation tensor
  LINALG::Matrix<3,3> U_mod;  // modified right stretch tensor
  LINALG::Matrix<3,3> U_disp; // displacement-based right stretch tensor
  LINALG::Matrix<3,3> EW;     // temporarily store eigenvalues
  LINALG::Matrix<3,3> tmp;    // temporary matrix for matrix matrix matrix products
  LINALG::Matrix<3,3> tmp2;    // temporary matrix for matrix matrix matrix products

  // ******************************************************************
  // calculate modified right stretch tensor
  // ******************************************************************
  for (int i=0; i<3; i++)
    U_mod(i,i) = 2.*glstrain_mod(i) + 1.;
  U_mod(0,1) = glstrain_mod(3);
  U_mod(1,0) = glstrain_mod(3);
  U_mod(1,2) = glstrain_mod(4);
  U_mod(2,1) = glstrain_mod(4);
  U_mod(0,2) = glstrain_mod(5);
  U_mod(2,0) = glstrain_mod(5);

  LINALG::SYEV(U_mod,EW,U_mod);
  for (int i=0; i<3; ++i)
    EW(i,i) = sqrt(EW(i,i));
  tmp.Multiply(U_mod,EW);
  tmp2.MultiplyNT(tmp,U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(defgrd_disp,defgrd_disp);

  LINALG::SYEV(U_disp,EW,U_disp);
  for (int i=0; i<3; ++i)
    EW(i,i) = sqrt(EW(i,i));
  tmp.Multiply(U_disp,EW);
  tmp2.MultiplyNT(tmp,U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(defgrd_disp,U_disp);
  defgrd_mod.Multiply(R,U_mod);

  // you're done here
  return;

}


template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex18>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>;
