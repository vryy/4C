/*----------------------------------------------------------------------*/
/*!
\file elast_coupvarga.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupVarga MUE 200.0

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupvarga.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupVarga::CoupVarga(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE")),
  beta_(matdata->GetDouble("BETA"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupVarga::CoupVarga()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupVarga::CoupVarga(MAT::ELASTIC::PAR::CoupVarga* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupVarga::AddShearMod(
  bool& haveshearmod,  ///< non-zero shear modulus was added
  double& shearmod  ///< variable to add upon
  ) const
{
  // indeed, a shear modulus is provided
  haveshearmod = haveshearmod or true;

  // material parameters for isochoric part
  shearmod += params_->mue_;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupVarga::AddCoefficientsStretchesPrincipal(
  LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
  LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
  const LINALG::Matrix<3,1>& prstr  ///< principal stretches, [lambda_1, lambda_2, lambda_3]
  )
{
#if 0
  // here, the isochoric neo-Hooke material in principal stretches is 
  // implemented to verify implementation
  
  // energy density
  //   Psi = mu [ \bar{lam}_1^2 + \bar{lam}_2^2 + \bar{lam}_3^2 - 3 ]
  //       = mu [ (J^{-1/3} lam_1)^2 + (J^{-1/3} lam_2)^2 + (J^{-1/3} lam_3)^2 - 3 ]
  //       = mu [ J^{-2/3} (lam_1^2 + lam_2^2 + lam_3^2) - 3 ]
  //       = mu [ J^{-2/3} I_C - 3 ]
  // with
  //   mu ... shear modulus
  //   lam_\alpha ... principal stretches, \alpha=1,2,3
  //   \bar{lam}_\alpha = J^{-1/3} lam_\alpha ... modified principal stretches
  //   I_C ... 1st invariant of right Cauchy-Green 2-tensor
  //   J = lam_1*lam_2*lam_3 ... determinant of deformation gradient

  // shear modulus
  const double mue = params_->mue_;
  // determinant of deformation gradient
  const double jac = prstr(0)*prstr(1)*prstr(2);
  const double jac23 = std::pow(jac,-2.0/3.0);
  const double jac53 = std::pow(jac,-5.0/3.0);
  const double jac83 = std::pow(jac,-8.0/3.0);
  // 1st invariant of right Cauchy-Green tensor
  const double fst = prstr(0)*prstr(0) + prstr(1)*prstr(1) + prstr(2)*prstr(2);
  
  // first derivatives
  gamma(0)  // ,0
    += mue*jac23*prstr(0)
    + (-1.0/3.0)*mue*fst*jac53*prstr(1)*prstr(2);
  gamma(1)  // ,1
    += mue*jac23*prstr(1)
    + (-1.0/3.0)*mue*fst*jac53*prstr(0)*prstr(2);
  gamma(2)  // ,2
    += mue*jac23*prstr(2)
    + (-1.0/3.0)*mue*fst*jac53*prstr(0)*prstr(1);

  // second derivatives
  delta(0)  // ,00
    += mue*jac23 
    + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(1)*prstr(2)
    + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(1)*prstr(2)
    + (5.0/9.0)*mue*fst*jac83*prstr(1)*prstr(2)*prstr(1)*prstr(2);
  delta(1)  // ,11
    += mue*jac23 
    + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(2)
    + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(2)
    + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(2)*prstr(0)*prstr(2);
  delta(2)  // ,22
    += mue*jac23 
    + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(1)
    + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(1)
    + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(1)*prstr(0)*prstr(1);
  delta(3)  // ,01
    += (-2.0/3.0)*mue*prstr(0)*jac53*prstr(0)*prstr(2)
    + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(1)*prstr(2)
    + (5.0/9.0)*mue*fst*jac83*prstr(1)*prstr(2)*prstr(0)*prstr(2)
    + (-1.0/3.0)*mue*fst*jac53*prstr(2);
  delta(4)  // ,12
    += (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(1)
    + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(2)
    + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(2)*prstr(0)*prstr(1)
    + (-1.0/3.0)*mue*fst*jac53*prstr(0);
  delta(5)  // ,20
    += (-2.0/3.0)*mue*prstr(2)*jac53*prstr(1)*prstr(2)
    + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(0)*prstr(1)
    + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(1)*prstr(1)*prstr(2)
    + (-1.0/3.0)*mue*fst*jac53*prstr(1);

#else
  // parameters
  const double alpha = 2.0*params_->mue_ - params_->beta_;
  const double beta = params_->beta_;

  // first derivatives
  // \frac{\partial Psi}{\partial \lambda_1}
  gamma(0) += alpha - beta/(prstr(0)*prstr(0));
  // \frac{\partial Psi}{\partial \lambda_2}
  gamma(1) += alpha - beta/(prstr(1)*prstr(1));
  // \frac{\partial Psi}{\partial \lambda_3}
  gamma(2) += alpha - beta/(prstr(2)*prstr(2));

  // second derivatives
  // \frac{\partial^2 Psi}{\partial\lambda_1^2}
  delta(0) += 2.0*beta/(prstr(0)*prstr(0)*prstr(0));
  // \frac{\partial^2 Psi}{\partial\lambda_2^2}
  delta(1) += 2.0*beta/(prstr(1)*prstr(1)*prstr(1));
  // \frac{\partial^2 Psi}{\partial\lambda_3^2}
  delta(2) += 2.0*beta/(prstr(2)*prstr(2)*prstr(2));
  // \frac{\partial^2 Psi}{\partial\lambda_1 \partial\lambda_2}
  delta(3) += 0.0;
  // \frac{\partial^2 Psi}{\partial\lambda_2 \partial\lambda_3}
  delta(4) += 0.0;
  // \frac{\partial^2 Psi}{\partial\lambda_3 \partial\lambda_1}
  delta(5) += 0.0;
#endif

  // done
  return;
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
