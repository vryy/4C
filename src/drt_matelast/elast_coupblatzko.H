/*----------------------------------------------------------------------*/
/*!
\file elast_coupblatzko.H
\brief
This file contains the routines required for Blatz and Ko material model
according to Holzapfel, "Nonlinear solid mechanics", 2001.

<pre>
Maintainer: Sophie Rausch
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPBLATZKO_H_
#define ELAST_COUPBLATZKO_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
namespace ELAST
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for Blatz and Ko material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupBlatzKo MUE 1.044E7 NUE 0.3 F 0.5
    class CoupBlatzKo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupBlatzKo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupBlatzKo() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      const double mue_;
      /// Possion's ratio
      const double nue_;
      /// interpolation parameter
      const double f_;

      //@}

    };  // class CoupBlatzKo

  }  // namespace PAR
  /*----------------------------------------------------------------------*/
  /// Coupled Blatz and Ko material according to [1].
  ///
  /// Strain energy function is given by
  ///
  ///
  ///\f[
  ///   \Psi= f \frac {\mu} 2 \left[ (I_{\boldsymbol C}-3)+\frac 1
  ///         {\beta} ( III_{\boldsymbol C}^{-\beta} -1) \right]
  ///         +(1-f) \frac {\mu} 2 \left[\left( \frac {II_{\boldsymbol
  ///         C}}{III_{\boldsymbol C}}-3 \right) + \frac 1 {\beta}
  ///         (III_{\boldsymbol C}^{\beta}-1)\right]
  ///\f]
  ///
  /// \f[
  ///    with \ \beta=\nue-0.5
  /// \f]
  ///
  ///
  /// [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
  ///     Journal of Applied Mechanics, 2000
  /// \author rausch
  /// \date 04/09
  class CoupBlatzKo : public Summand
  {

  public:

    /// empty constructor
    CoupBlatzKo();

    /// constructor with given material parameters
    CoupBlatzKo(MAT::ELAST::PAR::CoupBlatzKo* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupblatzko; }
    //@}
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// Add modified coeffiencts.
    /// Nothing to do for this (coupled) material!
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


  private:

    /// my material parameters
    MAT::ELAST::PAR::CoupBlatzKo* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // COUPBLATZKO_H_
