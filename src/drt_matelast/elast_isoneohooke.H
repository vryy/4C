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

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

      /// @name Write methods
      //@{
      void SetMue(const double mue) { mue_ = mue; }
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
  //         = \frac{\mu}{2} (J^{-2/3}{I}_{\boldsymbol{C}} - 3)
  ///\f]
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

    /// Add modified coeffiencts.
    virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );


    /// @name Access methods
    //@{
    double Mue() const { return params_->mue_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco ///< global indicator for isotropic splitted viscos formulation
        )
    {
      isomod = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISONEOHOOKE_H_
