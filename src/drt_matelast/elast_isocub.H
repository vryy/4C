/*----------------------------------------------------------------------*/
/*!
\file elast_isocub.H
\brief This file contains the routines required to calculate the isochoric contribution
of a cubic material.

\level 1

<pre>
\maintainer Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOCUB_H_
#define ELAST_ISOCUB_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

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
    /// MAT 1 ELAST_IsoCub C 1
    class IsoCub : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoCub(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoCub() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double c_;

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
      //@}

    };  // class IsoCub

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric cubic material
  ///
  /// This is the isochoric part of a third order hyperelastic, isotropic
  /// material depending on the modified trace of the right Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C (\overline{I}_{\boldsymbol{C}}-3)^3.
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 11/09
  class IsoCub : public Summand
  {

  public:

    /// constructor with given material parameters
    IsoCub(MAT::ELASTIC::PAR::IsoCub* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isocub; }

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
    MAT::ELASTIC::PAR::IsoCub* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISOCUB_H_
