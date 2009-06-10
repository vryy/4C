/*----------------------------------------------------------------------*/
/*!
\file elast_vologden.H
\brief
This file contains the routines required for the volumetic function
suggested by Ogden, see Doll, S. and Schweizerhof, K.
On the Development of Volumetric Strain Energy Functions
Journal of Applied Mechanics, 2000

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLOGDEN_H_
#define ELAST_VOLOGDEN_H_
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
    /// material parameters for volumetric contribution
    /// \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_VolOgden KAPPA 100 BETA -2
    class VolOgden : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VolOgden(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VolOgden() { ; }

      /// @name material parameters
      //@{

      /// Dilatation modulus
      const double kappa_;
      /// empiric constant
      const double beta_;

      //@}

    };  // class VolOgden

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Volumetric Ogden material according to [1].
  ///
  /// Strain energy function is given by
  ///
  ///
  ///\f[
  ///   \Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)
  ///\f]
  ///
  ///
  /// [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
  ///     Journal of Applied Mechanics, 2000
  /// \author rausch
  /// \date 04/09
  class VolOgden : public Summand
  {

  public:

    /// empty constructor
    VolOgden();

    /// constructor with given material parameters
    VolOgden(MAT::ELASTIC::PAR::VolOgden* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_vologden; }

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



  private:

    /// my material parameters
   MAT::ELASTIC::PAR::VolOgden* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // VOLOGDEN_H_
