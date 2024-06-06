/*----------------------------------------------------------------------*/
/*! \file
\brief Everything concerning EAS technology for so3_ssn_plast
       Mainly copied from so_hex8_eas.cpp. Redundancy needed,
       because of hard coded Gauss point in the so_hex8
       which do not coincide with the intrepid Gauss points
       used in the so3_ssn_plast.
\level 2
*/

/*----------------------------------------------------------------------*
 | headers                                                  seitz 04/14 |
 *----------------------------------------------------------------------*/
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_so3_plast_ssn.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                           seitz 04/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::eas_init()
{
  switch (eastype_)
  {
    case soh8p_easnone:
      neas_ = PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easnone>::neas;
      break;
    case soh8p_easmild:
      neas_ = PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas;
      break;
    case soh8p_easfull:
      neas_ = PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas;
      break;
    case soh8p_eassosh8:
      neas_ = PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas;
      break;
    case soh18p_eassosh18:
      neas_ = PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas;
      break;
    default:
      FOUR_C_THROW("unknown EAS type");
  }

  if (eastype_ != soh8p_easnone)
  {
    KaaInv_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, neas_, true));
    Kad_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, numdofperelement_, true));
    feas_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_last_timestep_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_delta_over_last_timestep_ =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    alpha_eas_inc_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(neas_, true));
    Kba_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseMatrix>(
        numgpt_, Core::LinAlg::SerialDenseMatrix(5, neas_, true)));

    Core::ProblemType probtype = Global::Problem::Instance()->GetProblemType();
    if (probtype == Core::ProblemType::tsi)
    {
      KaT_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(neas_, nen_, true));
      KdT_eas_ = Teuchos::rcp(new Core::LinAlg::Matrix<numdofperelement_, nen_>);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  setup EAS data (private)                                seitz 04/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::eas_setup()
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
  set_t0inv_t()(0, 0) = jac_0()(0, 0) * jac_0()(0, 0);
  set_t0inv_t()(1, 0) = jac_0()(1, 0) * jac_0()(1, 0);
  set_t0inv_t()(2, 0) = jac_0()(2, 0) * jac_0()(2, 0);
  set_t0inv_t()(3, 0) = 2. * jac_0()(0, 0) * jac_0()(1, 0);
  set_t0inv_t()(4, 0) = 2. * jac_0()(1, 0) * jac_0()(2, 0);
  set_t0inv_t()(5, 0) = 2. * jac_0()(0, 0) * jac_0()(2, 0);

  set_t0inv_t()(0, 1) = jac_0()(0, 1) * jac_0()(0, 1);
  set_t0inv_t()(1, 1) = jac_0()(1, 1) * jac_0()(1, 1);
  set_t0inv_t()(2, 1) = jac_0()(2, 1) * jac_0()(2, 1);
  set_t0inv_t()(3, 1) = 2. * jac_0()(0, 1) * jac_0()(1, 1);
  set_t0inv_t()(4, 1) = 2. * jac_0()(1, 1) * jac_0()(2, 1);
  set_t0inv_t()(5, 1) = 2. * jac_0()(0, 1) * jac_0()(2, 1);

  set_t0inv_t()(0, 2) = jac_0()(0, 2) * jac_0()(0, 2);
  set_t0inv_t()(1, 2) = jac_0()(1, 2) * jac_0()(1, 2);
  set_t0inv_t()(2, 2) = jac_0()(2, 2) * jac_0()(2, 2);
  set_t0inv_t()(3, 2) = 2. * jac_0()(0, 2) * jac_0()(1, 2);
  set_t0inv_t()(4, 2) = 2. * jac_0()(1, 2) * jac_0()(2, 2);
  set_t0inv_t()(5, 2) = 2. * jac_0()(0, 2) * jac_0()(2, 2);

  set_t0inv_t()(0, 3) = jac_0()(0, 0) * jac_0()(0, 1);
  set_t0inv_t()(1, 3) = jac_0()(1, 0) * jac_0()(1, 1);
  set_t0inv_t()(2, 3) = jac_0()(2, 0) * jac_0()(2, 1);
  set_t0inv_t()(3, 3) = jac_0()(0, 0) * jac_0()(1, 1) + jac_0()(1, 0) * jac_0()(0, 1);
  set_t0inv_t()(4, 3) = jac_0()(1, 0) * jac_0()(2, 1) + jac_0()(2, 0) * jac_0()(1, 1);
  set_t0inv_t()(5, 3) = jac_0()(0, 0) * jac_0()(2, 1) + jac_0()(2, 0) * jac_0()(0, 1);


  set_t0inv_t()(0, 4) = jac_0()(0, 1) * jac_0()(0, 2);
  set_t0inv_t()(1, 4) = jac_0()(1, 1) * jac_0()(1, 2);
  set_t0inv_t()(2, 4) = jac_0()(2, 1) * jac_0()(2, 2);
  set_t0inv_t()(3, 4) = jac_0()(0, 1) * jac_0()(1, 2) + jac_0()(1, 1) * jac_0()(0, 2);
  set_t0inv_t()(4, 4) = jac_0()(1, 1) * jac_0()(2, 2) + jac_0()(2, 1) * jac_0()(1, 2);
  set_t0inv_t()(5, 4) = jac_0()(0, 1) * jac_0()(2, 2) + jac_0()(2, 1) * jac_0()(0, 2);

  set_t0inv_t()(0, 5) = jac_0()(0, 0) * jac_0()(0, 2);
  set_t0inv_t()(1, 5) = jac_0()(1, 0) * jac_0()(1, 2);
  set_t0inv_t()(2, 5) = jac_0()(2, 0) * jac_0()(2, 2);
  set_t0inv_t()(3, 5) = jac_0()(0, 0) * jac_0()(1, 2) + jac_0()(1, 0) * jac_0()(0, 2);
  set_t0inv_t()(4, 5) = jac_0()(1, 0) * jac_0()(2, 2) + jac_0()(2, 0) * jac_0()(1, 2);
  set_t0inv_t()(5, 5) = jac_0()(0, 0) * jac_0()(2, 2) + jac_0()(2, 0) * jac_0()(0, 2);

  // now evaluate T0^{-T} with solver
  Core::LinAlg::FixedSizeSerialDenseSolver<numstr_, numstr_, 1> solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(set_t0inv_t());
  int err2 = solve_for_inverseT0.Factor();
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) || (err2 != 0)) FOUR_C_THROW("Inversion of T0inv (Jacobian0) failed");

  // reset EAS matrices
  KaaInv_->shape(neas_, neas_);
  Kad_->shape(neas_, numdofperelement_);
  if (KaT_ != Teuchos::null) KaT_->shape(neas_, nen_);
  feas_->size(neas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Defgrd consistent with enhanced GL strain (private)     seitz 04/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::calc_consistent_defgrd()
{
  static Core::LinAlg::Matrix<numstr_, 1> glstrain_mod(false);
  glstrain_mod(0) = 0.5 * (rcg()(0, 0) - 1.0);
  glstrain_mod(1) = 0.5 * (rcg()(1, 1) - 1.0);
  glstrain_mod(2) = 0.5 * (rcg()(2, 2) - 1.0);
  glstrain_mod(3) = rcg()(0, 1);
  glstrain_mod(4) = rcg()(1, 2);
  glstrain_mod(5) = rcg()(2, 0);

  Core::LinAlg::Matrix<3, 3> R;       // rotation tensor
  Core::LinAlg::Matrix<3, 3> U_mod;   // modified right stretch tensor
  Core::LinAlg::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
  Core::LinAlg::Matrix<3, 3> EW;      // temporarily store eigenvalues
  Core::LinAlg::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
  Core::LinAlg::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

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

  Core::LinAlg::SYEV(U_mod, EW, U_mod);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_mod, EW);
  tmp2.MultiplyNT(tmp, U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(defgrd(), defgrd());

  Core::LinAlg::SYEV(U_disp, EW, U_disp);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_disp, EW);
  tmp2.MultiplyNT(tmp, U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(defgrd(), U_disp);
  set_defgrd_mod().Multiply(R, U_mod);

  // you're done here
  return;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::eas_shape(const int gp)
{
  std::vector<Core::LinAlg::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
  // build EAS interpolation matrix M, evaluated at the 8 GPs of so_hex8

  // fill up M at each gp
  if (eastype_ == soh8p_easmild)
  {
    static std::vector<Core::LinAlg::SerialDenseMatrix> M_mild(numgpt_);
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
        M_mild[i].shape(numstr_, neas_);
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
    static std::vector<Core::LinAlg::SerialDenseMatrix> M_full(numgpt_);
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
        M_full[i].shape(numstr_, neas_);
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
    static std::vector<Core::LinAlg::SerialDenseMatrix> M_sosh8(numgpt_);
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
        M_sosh8[i].shape(numstr_, neas_);
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
    FOUR_C_THROW("this EAS type not yet implemented");

  // transform EAS shape functions from parameter space to actual space
  set_m_eas().shape(numstr_, neas_);
  set_m_eas().putScalar(0.0);
  switch (eastype_)
  {
    case soh8p_easfull:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
          set_m_eas().values(), det_jac_0() / det_j(), t0inv_t().A(), (M_GP->at(gp)).values());
      break;
    case soh8p_easmild:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
          set_m_eas().values(), det_jac_0() / det_j(), t0inv_t().A(), (M_GP->at(gp)).values());
      break;
    case soh8p_eassosh8:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
          set_m_eas().values(), det_jac_0() / det_j(), t0inv_t().A(), (M_GP->at(gp)).values());
      break;
    case soh8p_easnone:
      break;
    default:
      FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
      break;
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::eas_enhance_strains()
{
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  static Core::LinAlg::Matrix<numstr_, 1> total_glstrain(false);
  total_glstrain(0) = 0.5 * (rcg()(0, 0) - 1.0);
  total_glstrain(1) = 0.5 * (rcg()(1, 1) - 1.0);
  total_glstrain(2) = 0.5 * (rcg()(2, 2) - 1.0);
  total_glstrain(3) = rcg()(0, 1);
  total_glstrain(4) = rcg()(1, 2);
  total_glstrain(5) = rcg()(2, 0);
  // add enhanced strains = M . alpha to GL strains to "unlock" element
  switch (eastype_)
  {
    case soh8p_easfull:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, m_eas().values(), alpha_eas_->values());
      break;
    case soh8p_easmild:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, m_eas().values(), alpha_eas_->values());
      break;
    case soh8p_eassosh8:
      Core::LinAlg::DenseFunctions::multiply<double, numstr_,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, 1>(
          1.0, total_glstrain.A(), 1.0, m_eas().values(), alpha_eas_->values());
      break;
    case soh8p_easnone:
      break;
    default:
      FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
      break;
  }

  for (int i = 0; i < nsd_; ++i) set_rcg()(i, i) = 2. * total_glstrain(i) + 1.;
  set_rcg()(0, 1) = set_rcg()(1, 0) = total_glstrain(3);
  set_rcg()(2, 1) = set_rcg()(1, 2) = total_glstrain(4);
  set_rcg()(0, 2) = set_rcg()(2, 0) = total_glstrain(5);

  // calculate deformation gradient consistent with modified GL strain tensor
  calc_consistent_defgrd();
}

template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
