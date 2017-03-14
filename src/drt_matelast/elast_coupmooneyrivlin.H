/*----------------------------------------------------------------------*/
/*!
\file elast_coupmooneyrivlin.H

\brief
This file contains the routines required to calculate the contribution
of a CoupMooneyRivlin material material.
The input line should read
  MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPMOONEYRIVLIN_H_
#define ELAST_COUPMOONEYRIVLIN_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a CoupMooneyRivlinan material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1
    class CoupMooneyRivlin : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupMooneyRivlin(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupMooneyRivlin() { ; }

      /// @name material parameters
      //@{

      /// material parameters


      double c1_;
      double c2_;
      double c3_;


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
      void SetC3(const double c3) { c3_ = c3; }
      //@}

    };  // class CoupMooneyRivlin

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// CoupMooneyRivlinan material
  ///
  ///
  ///
  ///
  /// Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.
  ///
  /// The underlying strain-energy function is (expressed in invariants I1 and I2):
  ///
  /// \f[
  ///     \Psi= c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + c_3 * (J - 1)^2
  ///\f]
  ///
  /// For references see Holzapfel p. 245,
  ///
  /// Parameters are \f$ c_1, c_2\f$ and \f$c_3\f$ as penalty factor to enforce incompressibility
  /// (shear modulus \f$\mu = (c_1 + c_2) / 2)\f$
  ///
  ///
  ///
  /// \author rausch
  /// \date 03/12


  class CoupMooneyRivlin : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupMooneyRivlin(MAT::ELASTIC::PAR::CoupMooneyRivlin* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupmooneyrivlin; }

    //@}

    // add strain energy
    virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy function
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const LINALG::Matrix<6,1>& glstrain,///< Green-Lagrange strain
        const int eleGID                   ///< element GID
        );

    virtual void AddDerivativesPrincipal(
        LINALG::Matrix<3,1>& dPI, ///< first derivative with respect to invariants
        LINALG::Matrix<6,1>& ddPII, ///< second derivative with respect to invariants
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const int eleGID                  ///< element GID
        );

    /// add the derivatives of a coupled strain energy functions associated with a purely isochoric deformation
    virtual void AddCoupDerivVol(const double j, double& dPj1, double& dPj2, double& dPj3);

    /// @name Access methods
    //@{
    double C1() const { return params_->c1_; }
    double C2() const { return params_->c2_; }
    double C3() const { return params_->c3_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
        )
    {
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupMooneyRivlin* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPMOONEYRIVLIN_H_
