/*----------------------------------------------------------------------*/
/*!
\file elast_isoneohooke.H
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOTWO_H_
#define ELAST_COUPANISOEXPOTWO_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a neo-Hooke material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 CoupAnisoExpoTwo MUE 100
    class CoupAnisoExpoTwo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAnisoExpoTwo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoExpoTwo() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double k1_;
      double k2_;
      double k3_;
      double k4_;

      //@}

    };  // class CoupAnisoExpoTwo

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Coupled anisotropic exponential fiber function, implemented for two possible fiber families
  /// [1]
  /// This is a hyperelastic, anisotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac {k_1}{2 k_2} \left(e^{k_2 (IV_{\boldsymbol C}-1)^2 }-1 \right)+\frac {k_3}{2 k_4} \left(e^{k_4 (VI_{\boldsymbol C}-1)^2 }-1 \right).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, TC Gasser, A viscoelastic model for fiber-reinforced composites at finite strains: ... 2000
  /// </ul>
  ///
  /// \author rausch,tk
  /// \date 05/09
  class CoupAnisoExpoTwo : public Summand
  {

  public:

    /// empty constructor
    CoupAnisoExpoTwo();

    /// constructor with given material parameters
    CoupAnisoExpoTwo(MAT::ELASTIC::PAR::CoupAnisoExpoTwo* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoexpotwo; }

    //@}

    /// Add principal coeffiencts.
    /// Nothing to do for this (decoupled) material!
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;

      return;
    };

    /// Add modified coeffiencts.
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add anisotropic principal coeffiencts.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    );

    /// Add anisotropic modified coeffiencts.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsModifiedAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// @name Access methods
    //@{
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }
    double K3() const { return params_->k3_; }
    double K4() const { return params_->k4_; }

    //@}

    /// @name Write methods
    //@{
    void SetK1(const double k1) { params_->k1_ = k1; }
    void SetK2(const double k2) { params_->k2_ = k2; }
    void SetK3(const double k3) { params_->k3_ = k3; }
    void SetK4(const double k4) { params_->k4_ = k4; }
    //@}

    private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupAnisoExpoTwo* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // ISONEOHOOKE_H_
