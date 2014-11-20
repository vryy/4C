/*----------------------------------------------------------------------*/
/*!
\file elast_isotestmaterial.H

\brief
This file contains the routines required to calculate the isochoric contribution
of a Material, which is not realistic, but contains all possible derivations of invariants.
With this material in combination with volsussmannbathe, it is possible to test all
isochoric parts of the Elasthyper-Toolbox.

<pre>
Maintainer: Anna Birzle
            birzle@lnm.mw.tum.de
            089/289 15255
</pre>


*----------------------------------------------------------------------*/
#ifndef ELAST_ISOTESTMATERIAL_H_
#define ELAST_ISOTESTMATERIAL_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of test material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoTestMaterial C1 100 C2 50
    class IsoTestMaterial : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoTestMaterial(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoTestMaterial() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double c1_;
      double c2_;

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
      void SetC1(const double c1) { c1_ = c1; }
      void SetC2(const double c2) { c2_ = c2; }
      //@}

    };  // class IsoTestMaterial

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric Material to test the elasthyper-Toolbox
  ///
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + 0.5 C1 (\overline{I}_{\boldsymbol{C}}-3)^2
  ///        + C2 (\overline{II}_{\boldsymbol{C}}-3)  + 0.5 C2 (\overline{II}_{\boldsymbol{C}}-3)^2
  ///        + D (\overline{I}_{\boldsymbol{C}}-3) (\overline{II}_{\boldsymbol{C}}-3).
  ///\f]
  ///
  /// with D = C1 + 2 C2
  ///
  /// \author birzle
  /// \date 11/2014
  class IsoTestMaterial : public Summand
  {

  public:

    /// constructor with given material parameters
    IsoTestMaterial(MAT::ELASTIC::PAR::IsoTestMaterial* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isotestmaterial; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );

//    /// Add derivations with respect to modified invariants.
//    virtual void AddDerivationsModified(
//        LINALG::Matrix<3,1>& dPmodI, ///< first derivative with respect to modified invariants
//        LINALG::Matrix<6,1>& ddmodPII, ///< second derivative with respect to modified invariants
//        const LINALG::Matrix<3,1>& modinv ///< modified invariants of right Cauchy-Green tensor
//    );

    /// @name Access methods
    //@{
    double C1() const { return params_->c1_; }
    double C2() const { return params_->c2_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco, ///< global indicator for isotropic splitted viscos formulation
        bool& viscogenmax  ///< global indicator for viscous contribution according the SLS-Model
        )
    {
      isomod = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoTestMaterial* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif
