/*----------------------------------------------------------------------*/
/*!
\file elast_couplogmixneohooke.H

\brief
This file contains the routines required for mixed-logarithmic neo-Hooke material
according to example 1 in Stupkiewicz "An ALE formulation for implicit time integration of quasi-
steady state wear problems" Comput. Methods Appl. Mech. Engrg. 260, 2013.
The input line should read either
  MAT 1 ELAST_CoupLogMixNeoHooke MODE YN C1 1.044E7 C2 0.3
  ( C1 = Young, c2 = nue)
or
  MAT 1 ELAST_CoupLogMixNeoHooke MODE Lame C1 1. C2 1.

\level 1

<pre>
\maintainer Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_LOGMIXNEOHOOKE_H_
#define ELAST_LOGMIXNEOHOOKE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for mixed logarithmic neo-Hooke material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupLogMixNeoHooke MODE YN C1 1.0 C2 0.3
    /// MAT 1 ELAST_CoupLogMixNeoHooke MODE Lame C1 1.0 C2 1.0
    class CoupLogMixNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupLogMixNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupLogMixNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double mue_;
      /// Lame's constant
      double lambda_;

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
      void SetMue    (const double mue    ) { mue_     =  mue    ; }
      void SetLambda (const double lambda ) { lambda_  =  lambda ; }
      //@}

    };  // class CoupLogNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Logarithmic mixed neo-Hooke material:
  ///
  /// This is a compressible, hyperelastic, isotropic material.
  ///
  /// The material strain energy density function is
  ///\f[
  /// \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
  ///       - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
  ///      + \frac{\lambda}{2} \big( \sqrt{I\!I\!I_{\boldsymbol{C}}} -1 \big)^2
  ///\f]
  ///

  class CoupLogMixNeoHooke : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupLogMixNeoHooke(MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_couplogmixneohooke; }

    /// add shear modulus equivalent
    virtual void AddShearMod(
      bool& haveshearmod,  ///< non-zero shear modulus was added
      double& shearmod  ///< variable to add upon
      ) const;

    //@}

    /// Main material call to determine  2nd PK stress and material constitutive tensor
    ///
    /// The material strain energy density function is shown above.
    ///
    /// Here is \f$I_{\boldsymbol{C}}\f$ the first principal invariant
    /// of the right Cauchy--Green strain tensor \f$\boldsymbol{C}\f$
    /// and \f$I\!I\!I_{\boldsymbol{C}}\f$ its third.
    /// The isochoric part is proportional to \f$\mu\f$,
    /// whereas the \f$\lambda\f$-proportional part constitutes the volumetric contribution.
    ///
    /// The parameters are the shear modulus. In #params_->parmode_==0 it is set directly,
    /// in #params_->parmode_==1 it is computed by
    ///\f[
    /// \mu = \frac{E}{2(1+\nu)}
    ///\f]
    /// and Lame's coefficient, again for #params_->parmode_==1, otherwise set by user,
    ///\f[
    /// \lambda = \left\{\begin{array}{ll}
    ///      \frac{E \nu}{(1+\nu) (1-2\nu)} & \text{if $\nu \neq 1/2$}
    ///   \\ 0                              & \text{else}
    /// \end{array}\right.
    ///\f]
    ///

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
    virtual void AddCoupDerivVol(const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
    {dserror("not implemented");}

    /// @name Access methods
    //@{
    double Mue() const { return params_->mue_; }
    double Lambda() const { return params_->lambda_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc,    ///< global indicator for isotropic principal formulation
        bool& isomod,      ///< global indicator for isotropic splitted formulation
        bool& anisoprinc,  ///< global indicator for anisotropic principal formulation
        bool& anisomod,    ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
      )
    {
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_LOGMIXNEOHOOKE_H_
