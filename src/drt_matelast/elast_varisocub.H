/*----------------------------------------------------------------------*/
/*!
\file elast_varisocub.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a cubic material.

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VARISOCUB_H_
#define ELAST_IOSCUB_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a cubic material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_VarIsoCub FRAC 0.5 C 1
    class VarIsoCub : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VarIsoCub(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VarIsoCub() { ; }

      /// @name material parameters
      //@{

      /// Fraction parameter
      double frac_;
      /// Shear modulus
      double c_;

      //@}

    };  // class VarIsoCub

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric cubic material
  ///
  /// This is the isochoric part of a third order hyperelastic, isotropic
  /// material depending on the modified trace of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = frac* C (\overline{I}_{\boldsymbol{C}}-3)^3.
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 11/09
  class VarIsoCub : public Summand
  {

  public:

    /// empty constructor
    VarIsoCub();

    /// constructor with given material parameters
    VarIsoCub(MAT::ELASTIC::PAR::VarIsoCub* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_varisocub; }

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
    );

    /// Add anisotropic principal coeffiencts.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

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
    double Frac() const { return params_->frac_; }
    double C() const { return params_->c_; }

    //@}

    /// @name Write methods
    //@{
    void SetFrac(const double frac) { params_->frac_ = frac; }
    void SetC(const double c) { params_->c_ = c; }
    //@}

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::VarIsoCub* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // VARISOCUB_H_
