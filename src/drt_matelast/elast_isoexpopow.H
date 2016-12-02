/*----------------------------------------------------------------------*/
/*!
\file elast_isoexpopow.H

\brief
This file contains the routines required to calculate the isochoric contribution
of an exponential material.
The input line should read
MAT 1 ELAST_IsoExpoPow K1 5000. K2 5. C 1.
C = Exponent D

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOEXPOPOW_H_
#define ELAST_ISOEXPOPOW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

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
    /// MAT 1 ELAST_IsoExpoPow K1 5000. K2 5. C 1.
    /// C = Exponent D
    class IsoExpoPow : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoExpoPow(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoExpoPow() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double k1_;
      double k2_;
      /// exponent
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
      void SetK1(const double k1) { k1_ = k1; }
      void SetK2(const double k2) { k2_ = k2; }
      //@}

    };  // class IsoExpoPow

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric exponential material
  ///
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = \frac{k_1}{2k_2} (e^{k_2 (\overline{I}_{\boldsymbol{C}}-3)^d}-1).
  ///\f]
  ///
  ///
  /// \author rausch
  /// \date 08/09
  class IsoExpoPow : public Summand
  {

  public:

    /// constructor with given material parameters
    IsoExpoPow(MAT::ELASTIC::PAR::IsoExpoPow* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoexpopow; }

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
        LINALG::Matrix<3,1>& dPmodI, ///< first derivative with respect to modified invariants
        LINALG::Matrix<6,1>& ddmodPII, ///< second derivative with respect to modified invariants
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const int eleGID                   ///< element GID
        );

    /// @name Access methods
    //@{
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }

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
    MAT::ELASTIC::PAR::IsoExpoPow* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif
