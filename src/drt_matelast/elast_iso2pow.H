/*----------------------------------------------------------------------*/
/*!
\file elast_iso2pow.H
\brief This file contains the routines required to calculate the isochoric contribution
of a general power-type material.

\level 1

<pre>
\maintainer Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISO2POW_H_
#define ELAST_ISO2POW_H_

#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"

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

    /// constructor with given material parameters
    Iso2Pow(MAT::ELASTIC::PAR::Iso2Pow* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_iso2pow; }

    //@}

    // add strain energy
    virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy function
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const LINALG::Matrix<6,1> glstrain,///< Green-Lagrange strain
        const int eleGID                   ///< element GID
        );

    // Add derivatives with respect to modified invariants.
    virtual void AddDerivativesModified(
        LINALG::Matrix<3,1>& dPmodI,       ///< first derivative with respect to modified invariants
        LINALG::Matrix<6,1>& ddmodPII,     ///< second derivative with respect to modified invariants
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const int eleGID                   ///< element GID
        );
    /// @name Access methods
    //@{
    double C() const { return params_->c_; }
    int D() const { return params_->d_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< general indicator, if one viscoelastic formulation is used
      )
    {
      isomod = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::Iso2Pow* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISO2POW_H_
