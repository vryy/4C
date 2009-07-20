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
#ifndef ELAST_ISONEOHOOKE_H_
#define ELAST_ISONEOHOOKE_H_
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
    /// MAT 1 ELAST_IsoNeoHooke MUE 100
    class IsoNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double mue_;

      //@}

    };  // class IsoNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric neo-Hooke material according to [1]
  ///
  /// This is a hyperelastic, isotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac \mu 2 (\overline{I}_{\boldsymbol{C}}-3).
  ///\f]
  ///
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk
  /// \date 04/09
  class IsoNeoHooke : public Summand
  {

  public:

    /// empty constructor
    IsoNeoHooke();

    /// constructor with given material parameters
    IsoNeoHooke(MAT::ELASTIC::PAR::IsoNeoHooke* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoneohooke; }

    /// add shear modulus equivalent
    virtual void AddShearMod(
      bool& haveshearmod,  ///< non-zero shear modulus was added
      double& shearmod  ///< variable to add upon
      ) const;

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
    double Mue() const { return params_->mue_; }

    //@}

    /// @name Write methods
    //@{
    void SetMue(const double mue) { params_->mue_ = mue; }
    //@}
  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // ISONEOHOOKE_H_
