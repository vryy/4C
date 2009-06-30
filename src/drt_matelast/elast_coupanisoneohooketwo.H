/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooketwo.H
\brief
This file contains the routines required for the anisotropic material
using two  neo-Hooke materials in both directions according to

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISONEOHOOKETWO_H_
#define ELAST_COUPANISONEOHOOKETWO_H_
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
    /// material parameters for anisochoric contribution of a
    /// neo-Hooke material with two fiber directions
    ///
    ///
    /// <h3>Input line</h3>
    /// MAT 1 CoupAnisoNeoHookeTwo C1 100 C2 100
    class CoupAnisoNeoHookeTwo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAnisoNeoHookeTwo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAnisoNeoHookeTwo() { ; }

      /// @name material parameters
      //@{

      /// fiber params
      double c1_;
      double c2_;

      //@}

    };  // class CoupAnisoNeoHookeTwo

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Coupled anisotropic exponential fiber function, implemented for two possible fiber families
  /// [1]
  /// This is a hyperelastic, anisotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = c_1 \left(IV_{\boldsymbol C}-1\right)+c_2 \left(VI_{\boldsymbol C}-1\right).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, TC Gasser, A viscoelastic model for fiber-reinforced composites at finite strains: ... 2000
  /// </ul>
  ///
  /// \author rausch
  /// \date 05/09
  class CoupAnisoNeoHookeTwo : public Summand
  {

  public:

    /// empty constructor
    CoupAnisoNeoHookeTwo();

    /// constructor with given material parameters
    CoupAnisoNeoHookeTwo(MAT::ELASTIC::PAR::CoupAnisoNeoHookeTwo* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupanisoneohooketwo; }

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
    double C1() const { return params_->c1_; }
    double C2() const { return params_->c2_; }

    //@}

    /// @name Write methods
    //@{
    void SetC1(const double c1) { params_->c1_ = c1; }
    void SetC2(const double c2) { params_->c2_ = c2; }
    //@}

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupAnisoNeoHookeTwo* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // COUPANISONEOHOOKETWO_H_
