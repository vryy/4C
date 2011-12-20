/*----------------------------------------------------------------------*/
/*!
\file elast_varisoneohooke.H
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
#ifndef ELAST_VARISONEOHOOKE_H_
#define ELAST_VARISONEOHOOKE_H_
#ifdef CCADISCRET

// Linalg matrices used instead; included in elast_summand
//#include <Epetra_SerialDenseVector.h>
//#include <Epetra_SerialDenseMatrix.h>
// included in elast_summand
//#include "../drt_inpar/inpar_material.H"
#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

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
    /// MAT 1 ELAST_VarIsoNeoHooke FRAC 0.5 MUE 100
    class VarIsoNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VarIsoNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VarIsoNeoHooke() { ; }

      /// @name material parameters
      //@{


      /// Fraction of continuen still there
      double frac_;

      /// Shear modulus
      double mue_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class VarIsoNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric neo-Hooke material according to [1]
  ///
  /// added the factor frac to account the fact that after collagease
  /// elastase a certain proportion of the corresponding fibers is
  /// destroid and has to be excluded from the SEF
  ///
  /// This is a hyperelastic, isotropic material
  /// of the most simple kind.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = frac \frac \mu 2 (\overline{I}_{\boldsymbol{C}}-3).
  //         = frac \frac{\mu}{2} (J^{-2/3}{I}_{\boldsymbol{C}} - 3)
  ///\f]
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch
  /// \date 04/09
  class VarIsoNeoHooke : public Summand
  {

  public:

    /// empty constructor
    VarIsoNeoHooke();

    /// constructor with given material parameters
    VarIsoNeoHooke(MAT::ELASTIC::PAR::VarIsoNeoHooke* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_varisoneohooke; }

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

    /// @name Access methods
    //@{
    double Frac() const { return params_->frac_; }
    double Mue()  const { return params_->mue_; }

    //@}

    /// @name Write methods
    //@{
    void SetFrac(const double frac) { params_->frac_ = frac; }
    void SetMue(const double mue)   { params_->mue_ = mue; }
    //@}
  private:

    /// my material parameters
    MAT::ELASTIC::PAR::VarIsoNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // VARISONEOHOOKE_H_
