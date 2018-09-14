/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_eas.cpp
\brief Everything concerning EAS technology for so3_ssn_plast
       Mainly copied from so_hex8_eas.cpp. Redundancy needed,
       because of hard coded Gauss point in the so_hex8
       which do not coincide with the intrepid Gauss points
       used in the so3_ssn_plast.
\level 2
\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------*
 | headers                                                  seitz 04/14 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"
#include "../../linalg/linalg_utils.H"
#include "../../linalg/linalg_fixedsizematrix.H"
#include "../../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                           seitz 04/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasInit()
{
  switch (eastype_)
  {
    case soh8p_easnone:
      neas_ = PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easnone>::neas;
      break;
    case soh8p_easmild:
      neas_ = PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas;
      break;
    case soh8p_easfull:
      neas_ = PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas;
      break;
    case soh8p_eassosh8:
      neas_ = PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas;
      break;
    case soh18p_eassosh18:
      neas_ = PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas;
      break;
    default:
      dserror("unknown EAS type");
  }

  if (eastype_ != soh8p_easnone)
  {
    KaaInv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_, neas_, true));
    Kad_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_, numdofperelement_, true));
    feas_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_last_timestep_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_delta_over_last_timestep_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    alpha_eas_inc_ = Teuchos::rcp(new LINALG::SerialDenseVector(neas_, true));
    Kba_ = Teuchos::rcp(new std::vector<LINALG::SerialDenseMatrix>(
        numgpt_, LINALG::SerialDenseMatrix(5, neas_, true)));

    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    if (probtype == prb_tsi)
    {
      KaT_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_, nen_, true));
      KdT_eas_ = Teuchos::rcp(new LINALG::Matrix<numdofperelement_, nen_>);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  setup EAS data (private)                                seitz 04/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasSetup()
{
  /* evaluation of EAS variables (which are constant for the following):
  ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
  ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
  ** -> T0^{-T}
  */

  // first, build T0^T transformation matrix which maps the M-matrix
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  SetT0invT()(0, 0) = Jac_0()(0, 0) * Jac_0()(0, 0);
  SetT0invT()(1, 0) = Jac_0()(1, 0) * Jac_0()(1, 0);
  SetT0invT()(2, 0) = Jac_0()(2, 0) * Jac_0()(2, 0);
  SetT0invT()(3, 0) = 2. * Jac_0()(0, 0) * Jac_0()(1, 0);
  SetT0invT()(4, 0) = 2. * Jac_0()(1, 0) * Jac_0()(2, 0);
  SetT0invT()(5, 0) = 2. * Jac_0()(0, 0) * Jac_0()(2, 0);

  SetT0invT()(0, 1) = Jac_0()(0, 1) * Jac_0()(0, 1);
  SetT0invT()(1, 1) = Jac_0()(1, 1) * Jac_0()(1, 1);
  SetT0invT()(2, 1) = Jac_0()(2, 1) * Jac_0()(2, 1);
  SetT0invT()(3, 1) = 2. * Jac_0()(0, 1) * Jac_0()(1, 1);
  SetT0invT()(4, 1) = 2. * Jac_0()(1, 1) * Jac_0()(2, 1);
  SetT0invT()(5, 1) = 2. * Jac_0()(0, 1) * Jac_0()(2, 1);

  SetT0invT()(0, 2) = Jac_0()(0, 2) * Jac_0()(0, 2);
  SetT0invT()(1, 2) = Jac_0()(1, 2) * Jac_0()(1, 2);
  SetT0invT()(2, 2) = Jac_0()(2, 2) * Jac_0()(2, 2);
  SetT0invT()(3, 2) = 2. * Jac_0()(0, 2) * Jac_0()(1, 2);
  SetT0invT()(4, 2) = 2. * Jac_0()(1, 2) * Jac_0()(2, 2);
  SetT0invT()(5, 2) = 2. * Jac_0()(0, 2) * Jac_0()(2, 2);

  SetT0invT()(0, 3) = Jac_0()(0, 0) * Jac_0()(0, 1);
  SetT0invT()(1, 3) = Jac_0()(1, 0) * Jac_0()(1, 1);
  SetT0invT()(2, 3) = Jac_0()(2, 0) * Jac_0()(2, 1);
  SetT0invT()(3, 3) = Jac_0()(0, 0) * Jac_0()(1, 1) + Jac_0()(1, 0) * Jac_0()(0, 1);
  SetT0invT()(4, 3) = Jac_0()(1, 0) * Jac_0()(2, 1) + Jac_0()(2, 0) * Jac_0()(1, 1);
  SetT0invT()(5, 3) = Jac_0()(0, 0) * Jac_0()(2, 1) + Jac_0()(2, 0) * Jac_0()(0, 1);


  SetT0invT()(0, 4) = Jac_0()(0, 1) * Jac_0()(0, 2);
  SetT0invT()(1, 4) = Jac_0()(1, 1) * Jac_0()(1, 2);
  SetT0invT()(2, 4) = Jac_0()(2, 1) * Jac_0()(2, 2);
  SetT0invT()(3, 4) = Jac_0()(0, 1) * Jac_0()(1, 2) + Jac_0()(1, 1) * Jac_0()(0, 2);
  SetT0invT()(4, 4) = Jac_0()(1, 1) * Jac_0()(2, 2) + Jac_0()(2, 1) * Jac_0()(1, 2);
  SetT0invT()(5, 4) = Jac_0()(0, 1) * Jac_0()(2, 2) + Jac_0()(2, 1) * Jac_0()(0, 2);

  SetT0invT()(0, 5) = Jac_0()(0, 0) * Jac_0()(0, 2);
  SetT0invT()(1, 5) = Jac_0()(1, 0) * Jac_0()(1, 2);
  SetT0invT()(2, 5) = Jac_0()(2, 0) * Jac_0()(2, 2);
  SetT0invT()(3, 5) = Jac_0()(0, 0) * Jac_0()(1, 2) + Jac_0()(1, 0) * Jac_0()(0, 2);
  SetT0invT()(4, 5) = Jac_0()(1, 0) * Jac_0()(2, 2) + Jac_0()(2, 0) * Jac_0()(1, 2);
  SetT0invT()(5, 5) = Jac_0()(0, 0) * Jac_0()(2, 2) + Jac_0()(2, 0) * Jac_0()(0, 2);

  // now evaluate T0^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<numstr_, numstr_, 1> solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(SetT0invT());
  int err2 = solve_for_inverseT0.Factor();
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) || (err2 != 0)) dserror("Inversion of T0inv (Jacobian0) failed");

  // reset EAS matrices
  KaaInv_->Shape(neas_, neas_);
  Kad_->Shape(neas_, numdofperelement_);
  if (KaT_ != Teuchos::null) KaT_->Shape(neas_, nen_);
  feas_->Size(neas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Defgrd consistent with enhanced GL strain (private)     seitz 04/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::CalcConsistentDefgrd()
{
  static LINALG::Matrix<numstr_, 1> glstrain_mod(false);
  glstrain_mod(0) = 0.5 * (RCG()(0, 0) - 1.0);
  glstrain_mod(1) = 0.5 * (RCG()(1, 1) - 1.0);
  glstrain_mod(2) = 0.5 * (RCG()(2, 2) - 1.0);
  glstrain_mod(3) = RCG()(0, 1);
  glstrain_mod(4) = RCG()(1, 2);
  glstrain_mod(5) = RCG()(2, 0);

  LINALG::Matrix<3, 3> R;       // rotation tensor
  LINALG::Matrix<3, 3> U_mod;   // modified right stretch tensor
  LINALG::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
  LINALG::Matrix<3, 3> EW;      // temporarily store eigenvalues
  LINALG::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
  LINALG::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

  // ******************************************************************
  // calculate modified right stretch tensor
  // ******************************************************************
  for (int i = 0; i < 3; i++) U_mod(i, i) = 2. * glstrain_mod(i) + 1.;
  U_mod(0, 1) = glstrain_mod(3);
  U_mod(1, 0) = glstrain_mod(3);
  U_mod(1, 2) = glstrain_mod(4);
  U_mod(2, 1) = glstrain_mod(4);
  U_mod(0, 2) = glstrain_mod(5);
  U_mod(2, 0) = glstrain_mod(5);

  LINALG::SYEV(U_mod, EW, U_mod);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_mod, EW);
  tmp2.MultiplyNT(tmp, U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(Defgrd(), Defgrd());

  LINALG::SYEV(U_disp, EW, U_disp);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_disp, EW);
  tmp2.MultiplyNT(tmp, U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(Defgrd(), U_disp);
  SetDefgrdMod().Multiply(R, U_mod);

  // you're done here
  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasShape(const int gp)
{
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;  // EAS matrix M at all GPs
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
    if (!M_mild_eval)  // if true M already evaluated
    {
      // fill up M at each gp
      for (int i = 0; i < numgpt_; ++i)
      {
        M_mild[i].Shape(numstr_, neas_);
        M_mild[i](0, 0) = xsi_.at(i)(0);
        M_mild[i](1, 1) = xsi_.at(i)(1);
        M_mild[i](2, 2) = xsi_.at(i)(2);

        M_mild[i](3, 3) = xsi_.at(i)(0);
        M_mild[i](3, 4) = xsi_.at(i)(1);
        M_mild[i](4, 5) = xsi_.at(i)(1);
        M_mild[i](4, 6) = xsi_.at(i)(2);
        M_mild[i](5, 7) = xsi_.at(i)(0);
        M_mild[i](5, 8) = xsi_.at(i)(2);
      }
      M_mild_eval = true;  // now the array is filled statically
    }

    // return adress of just evaluated matrix
    M_GP = &M_mild;  // return adress of static object to target of pointer
  }
  else if (eastype_ == soh8p_easfull)
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
    if (!M_full_eval)  // if true M already evaluated
    {
      // fill up M at each gp
      for (int i = 0; i < numgpt_; ++i)
      {
        M_full[i].Shape(numstr_, neas_);
        M_full[i](0, 0) = xsi_.at(i)(0);
        M_full[i](0, 15) = xsi_.at(i)(0) * xsi_.at(i)(1);
        M_full[i](0, 16) = xsi_.at(i)(0) * xsi_.at(i)(2);
        M_full[i](1, 1) = xsi_.at(i)(1);
        M_full[i](1, 17) = xsi_.at(i)(0) * xsi_.at(i)(1);
        M_full[i](1, 18) = xsi_.at(i)(1) * xsi_.at(i)(2);
        M_full[i](2, 2) = xsi_.at(i)(2);
        M_full[i](2, 19) = xsi_.at(i)(0) * xsi_.at(i)(2);
        M_full[i](2, 20) = xsi_.at(i)(1) * xsi_.at(i)(2);

        M_full[i](3, 3) = xsi_.at(i)(0);
        M_full[i](3, 4) = xsi_.at(i)(1);
        M_full[i](3, 9) = xsi_.at(i)(0) * xsi_.at(i)(2);
        M_full[i](3, 10) = xsi_.at(i)(1) * xsi_.at(i)(2);
        M_full[i](4, 5) = xsi_.at(i)(1);
        M_full[i](4, 6) = xsi_.at(i)(2);
        M_full[i](4, 11) = xsi_.at(i)(0) * xsi_.at(i)(1);
        M_full[i](4, 12) = xsi_.at(i)(0) * xsi_.at(i)(2);
        M_full[i](5, 7) = xsi_.at(i)(0);
        M_full[i](5, 8) = xsi_.at(i)(2);
        M_full[i](5, 13) = xsi_.at(i)(0) * xsi_.at(i)(1);
        M_full[i](5, 14) = xsi_.at(i)(1) * xsi_.at(i)(2);
      }
      M_full_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    M_GP = &M_full;  // return adress of static object to target of pointer
  }
  else if (eastype_ == soh8p_eassosh8)
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
    if (!M_sosh8_eval)  // if true M already evaluated
    {
      // fill up M at each gp
      for (int i = 0; i < numgpt_; ++i)
      {
        M_sosh8[i].Shape(numstr_, neas_);
        M_sosh8[i](0, 0) = xsi_.at(i)(0);
        M_sosh8[i](1, 1) = xsi_.at(i)(1);
        M_sosh8[i](2, 2) = xsi_.at(i)(2);
        M_sosh8[i](2, 5) = xsi_.at(i)(0) * xsi_.at(i)(2);
        M_sosh8[i](2, 6) = xsi_.at(i)(1) * xsi_.at(i)(2);

        M_sosh8[i](3, 3) = xsi_.at(i)(0);
        M_sosh8[i](3, 4) = xsi_.at(i)(1);
      }

      M_sosh8_eval = true;  // now the array is filled statically
    }
    // return adress of just evaluated matrix
    M_GP = &M_sosh8;  // return adress of static object to target of pointer
  }
  else
    dserror("this EAS type not yet implemented");

  // transform EAS shape functions from parameter space to actual space
  SetM_eas().LightShape(numstr_, neas_);
  SetM_eas().Zero();
  switch (eastype_)
  {
    case soh8p_easfull:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
          SetM_eas().A(), DetJac_0() / DetJ(), T0invT().A(), (M_GP->at(gp)).A());
      break;
    case soh8p_easmild:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
          SetM_eas().A(), DetJac_0() / DetJ(), T0invT().A(), (M_GP->at(gp)).A());
      break;
    case soh8p_eassosh8:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
          SetM_eas().A(), DetJac_0() / DetJ(), T0invT().A(), (M_GP->at(gp)).A());
      break;
    case soh8p_easnone:
      break;
    default:
      dserror("Don't know what to do with EAS type %d", eastype_);
      break;
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::EasEnhanceStrains()
{
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  static LINALG::Matrix<numstr_, 1> total_glstrain(false);
  total_glstrain(0) = 0.5 * (RCG()(0, 0) - 1.0);
  total_glstrain(1) = 0.5 * (RCG()(1, 1) - 1.0);
  total_glstrain(2) = 0.5 * (RCG()(2, 2) - 1.0);
  total_glstrain(3) = RCG()(0, 1);
  total_glstrain(4) = RCG()(1, 2);
  total_glstrain(5) = RCG()(2, 0);
  // add enhanced strains = M . alpha to GL strains to "unlock" element
  switch (eastype_)
  {
    case soh8p_easfull:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, M_eas().A(), alpha_eas_->A());
      break;
    case soh8p_easmild:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, M_eas().A(), alpha_eas_->A());
      break;
    case soh8p_eassosh8:
      LINALG::DENSEFUNCTIONS::multiply<double, numstr_,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, M_eas().A(), alpha_eas_->A());
      break;
    case soh8p_easnone:
      break;
    default:
      dserror("Don't know what to do with EAS type %d", eastype_);
      break;
  }

  for (int i = 0; i < nsd_; ++i) SetRCG()(i, i) = 2. * total_glstrain(i) + 1.;
  SetRCG()(0, 1) = SetRCG()(1, 0) = total_glstrain(3);
  SetRCG()(2, 1) = SetRCG()(1, 2) = total_glstrain(4);
  SetRCG()(0, 2) = SetRCG()(2, 0) = total_glstrain(5);

  // calculate deformation gradient consistent with modified GL strain tensor
  CalcConsistentDefgrd();
}

template class DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex18>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::nurbs27>;
