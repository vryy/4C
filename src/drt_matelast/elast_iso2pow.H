/*----------------------------------------------------------------------*/
/*!
\file elast_iso2pow.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a general power-type material.

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISO2POW_H_
#define ELAST_ISO2POW_H_
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
    /// material parameters for isochoric contribution of a general power material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_Iso2Pow C 1 D 1
    class Iso2Pow : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      Iso2Pow(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Iso2Pow() { ; }

      /// @name material parameters
      //@{

      /// material parameters
      double c_;
      int d_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetC(const double c) { c_ = c; }
      void SetD(const int d) { d_ = d; }
      //@}

    };  // class Iso2Pow

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric general power material
  ///
  /// This is the isochoric part of a variable order hyperelastic, isotropic
  /// material depending on the seccond invariant of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C (\overline{II}_{\boldsymbol{C}}-3)^D.
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 03/11
  class Iso2Pow : public Summand
  {

  public:

    /// empty constructor
    Iso2Pow();

    /// constructor with given material parameters
    Iso2Pow(MAT::ELASTIC::PAR::Iso2Pow* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_iso2pow; }

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
    double C() const { return params_->c_; }
    int D() const { return params_->d_; }

    //@}

    /// Indicator for anisotropic summand
    virtual void IsAnisotropic(
      bool& isaniso ///< global indicator for anisotropy
      )
    {
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::Iso2Pow* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // ISO2POW_H_
