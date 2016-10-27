/*----------------------------------------------------------------------*/
/*!
\file elast_coupneohooke.H

\brief bla bla bla

\level 1

\maintainer Sophie Rausch

\brief This file contains the routines required to calculate the isochoric
contribution of a CoupNeoHookean material material.
*----------------------------------------------------------------------*/
#ifndef ELAST_COUPNEOHOOKE_H_
#define ELAST_COUPNEOHOOKE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a CoupNeoHookean material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupNeoHooke YOUNG 1 NUE 1
    class CoupNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// material parameters


      /// Young's modulus
      double youngs_;
      /// Possion's ratio
      double nue_;

      /// nue \(1-2 nue)
      double beta_;
      /// Shear modulus / 2
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
      void SetBeta(const double beta) { beta_ = beta; }
      //@}

    };  // class CoupNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// CoupNeoHookean material
  ///
  /// This is the summand of a hyperelastic, isotropic CoupNeoHookean material
  /// depending on the first and the third invariant of the right Cauchy-Green tensor.
  /// The formulation is based on [1] p. 247,248 and 263
  ///
  /// The implemented material is the coupled form of the compressible NeoHook
  /// model. The Parameters read in are the Young's modulus and the Poisson's
  /// ratio.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = c(I_{\boldsymbol{C}}-3)+\frac {c}{\beta} (J^{-2\beta}-1)
  ///\f]
  ///
  /// with
  ///\f[
  ///   \beta = \frac {\nu}{1-2\nu}
  ///\f]
  ///
  /// \f$ c=\frac {\mu}{2} = \frac {\text{Young's modulus}}{4(1+\nu)} \f$
  /// and \f$\mu\f$ and \f$\nu\f$ denoting the shear modulus
  /// and the Poisson's ratio, respectively.
  ///
  /// [1] Holzapfel, G. A., Nonlinear Solid Mechanics, 2002
  ///
  /// \author rausch
  /// \date 03/11
  class CoupNeoHooke : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupNeoHooke(MAT::ELASTIC::PAR::CoupNeoHooke* params);


    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupneohooke; }

    /// add shear modulus equivalent
    virtual void AddShearMod(
      bool& haveshearmod,  ///< non-zero shear modulus was added
      double& shearmod  ///< variable to add upon
      ) const;

    /// add young's modulus equivalent
    virtual void AddYoungsMod(double& young, double& shear, double& bulk) {young+=YOUNGS();};

    //@}

    // add strain energy
    virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy function
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const LINALG::Matrix<6,1> glstrain,///< Green-Lagrange strain
        const int eleGID                   ///< element GID
        );

    virtual void AddDerivativesPrincipal(
        LINALG::Matrix<3,1>& dPI, ///< first derivative with respect to invariants
        LINALG::Matrix<6,1>& ddPII, ///< second derivative with respect to invariants
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const int eleGID                  ///< element GID
        );

    /// @name Access methods
    //@{
    double C() const { return params_->c_; }
    double BETA() const { return params_->beta_; }
    double NUE() const { return params_->nue_; }
    double YOUNGS() const { return params_->youngs_; }
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
    MAT::ELASTIC::PAR::CoupNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // COUPNEOHOOKE_H_
