/*----------------------------------------------------------------------*/
/*!
\file elast_isovarga.cpp
\brief


the input line should read
  MAT 1 ELAST_IsoVarga MUE 200.0

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
#include "elast_isovarga.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVarga::IsoVarga(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVarga::IsoVarga()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVarga::IsoVarga(MAT::ELASTIC::PAR::IsoVarga* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVarga::AddShearMod(
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
void MAT::ELASTIC::IsoVarga::AddCoefficientsStretchesModified(
  LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
  LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
  const LINALG::Matrix<3,1>& modstr  ///< modified principal stretches, [lambda_1, lambda_2, lambda_3]
  )
{
#if 0
  // here, the isochoric neo-Hooke material in principal stretches is 
  // implemented to verify implementation
  
  // energy density
  //   Psi = mu/2 [ \bar{lam}_1^2 + \bar{lam}_2^2 + \bar{lam}_3^2 - 3 ]
  //       = mu/2 [ (J^{-1/3} lam_1)^2 + (J^{-1/3} lam_2)^2 + (J^{-1/3} lam_3)^2 - 3 ]
  //       = mu/2 [ J^{-2/3} (lam_1^2 + lam_2^2 + lam_3^2) - 3 ]
  //       = mu/2 [ J^{-2/3} I_C - 3 ]
  // with
  //   mu ... shear modulus
  //   lam_\alpha ... principal stretches, \alpha=1,2,3
  //   \bar{lam}_\alpha = J^{-1/3} lam_\alpha ... modified principal stretches
  //   I_C ... 1st invariant of right Cauchy-Green 2-tensor
  //   J = lam_1*lam_2*lam_3 ... determinant of deformation gradient

  // shear modulus
  const double mue = params_->mue_;
  
  // first derivatives
  gamma(0)  // ,0
    += mue*modstr(0);
  gamma(1)  // ,1
    += mue*modstr(1);
  gamma(2)  // ,2
    += mue*modstr(2);

  // second derivatives
  delta(0)  // ,00
    += mue;
  delta(1)  // ,11
    += mue;
  delta(2)  // ,22
    += mue;
  delta(3)  // ,01
    += 0.0;
  delta(4)  // ,12
    += 0.0;
  delta(5)  // ,20
    += 0.0;

#else
  // first derivatives
  gamma(0) += 2.0*params_->mue_;  // \frac{\partial Psi}{\partial \bar{\lambda}_1}
  gamma(1) += 2.0*params_->mue_;  // \frac{\partial Psi}{\partial \bar{\lambda}_2}
  gamma(2) += 2.0*params_->mue_;  // \frac{\partial Psi}{\partial \bar{\lambda}_3}

  // second derivatives
  delta(0) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_1^2}
  delta(1) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_2^2}
  delta(2) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_3^2}
  delta(3) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_1 \partial\lambda_2}
  delta(4) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_2 \partial\lambda_3}
  delta(5) += 0.0;  // \frac{\partial^2 Psi}{\partial\lambda_3 \partial\lambda_1}
#endif

  // done
  return;
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
