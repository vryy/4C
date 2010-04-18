/*----------------------------------------------------------------------*/
/*!
\file elast_isoexpo.H
\brief

This file contains the routines required to calculate the isochoric contribution
of an exponential material

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOEXPO_H_
#define ELAST_ISOEXPO_H_
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
    /// material parameters for isochoric contribution of a exponential type material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoExpo K1 1. K2 1.
    class IsoExpo : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoExpo(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoExpo() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double k1_;
      double k2_;

      //@}

    };  // class IsoExpo

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric exponential material
  ///
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac{k_1}{2k_2} (e^{(\overline{I}_{\boldsymbol{C}}-3)^2}-1).
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 08/09
  class IsoExpo : public Summand
  {

  public:

    /// empty constructor
    IsoExpo();

    /// constructor with given material parameters
    IsoExpo(MAT::ELASTIC::PAR::IsoExpo* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoexpo; }

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
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }

    //@}

    /// @name Write methods
    //@{
    void SetK1(const double k1) { params_->k1_ = k1; }
    void SetK2(const double k2) { params_->k2_ = k2; }
    //@}
  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoExpo* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif
