/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.
The input line should read either
  MAT 1 ELAST_CoupLogNeoHooke YOUNG 1.044E7 NUE 0.3
or
  MAT 1 ELAST_CoupLogNeoHooke MUE 1. LAMBDA 1.

\level 1

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_LOGNEOHOOKE_H_
#define ELAST_LOGNEOHOOKE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for logarithmic neo-Hooke material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_CoupLogNeoHooke MODE YN C1 1.0 C2 0.3
      /// MAT 1 ELAST_CoupLogNeoHooke MODE Lame C1 1.0 C2 1.0
      class CoupLogNeoHooke : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupLogNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~CoupLogNeoHooke() { ; }

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
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        void SetMue(const double mue) { mue_ = mue; }
        void SetLambda(const double lambda) { lambda_ = lambda; }
        //@}

      };  // class CoupLogNeoHooke

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Logarithmic neo-Hooke material according to [1]
    ///
    /// This is a compressible, hyperelastic, isotropic material
    /// of the most simple kind.
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
    ///       - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
    ///      + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
    ///\f]
    /// which was taken from [1]. More details at #AddCoefficientsPrincipal()
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] J Bonet and RD Wood, "Nonlinear continuum mechanics for finite
    ///          element analysis", Cambridge, 1997.
    /// <li> [2] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
    /// </ul>
    ///
    /// \author bborn
    /// \date 04/09
    class CoupLogNeoHooke : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupLogNeoHooke(MAT::ELASTIC::PAR::CoupLogNeoHooke* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_couplogneohooke;
      }

      /// add shear modulus equivalent
      virtual void AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
          double& shearmod                          ///< variable to add upon
          ) const;

      //@}

      /// Main material call to determine  2nd PK stress and material constitutive tensor
      ///
      /// The material strain energy density function is
      ///\f[
      /// \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
      ///       - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
      ///      + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
      ///\f]
      /// which was taken from [1].
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
      /// The 2nd Piola--Kirchhoff stress is
      ///\f[
      /// \boldsymbol{S} = \mu \big( \boldsymbol{1} - \boldsymbol{C}^{-1} \big)
      ///                + \lambda \, \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \,  \boldsymbol{C}^{-1}
      ///\f]
      ///
      /// The material constitutive 4-tensor
      /// \f$\boldsymbol{C}_\text{m}=C_{IJKL}\boldsymbol{E}^I\otimes\boldsymbol{E}^J\otimes\boldsymbol{E}^K\otimes\boldsymbol{E}^L\f$
      /// is determined by
      ///\f[
      /// C_{IJKL}
      /// = \lambda (\boldsymbol{C}^{-1})_{IJ} \, (\boldsymbol{C}^{-1})_{KL}
      /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
      ///   \, \frac{1}{2} \big( \boldsymbol{C}^{-1})_{IK} \, (\boldsymbol{C}^{-1})_{JL}
      ///                        + (\boldsymbol{C}^{-1})_{IL} \, (\boldsymbol{C}^{-1})_{JK} \big)
      ///\f]
      /// or
      ///\f[
      /// \boldsymbol{C}_\text{m}
      /// = \lambda \boldsymbol{C}^{-1} \otimes  \boldsymbol{C}^{-1}
      /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
      ///   \, \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
      ///\f]

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          int eleGID                             ///< element GID
      );

      virtual void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,                             ///< Gauss point
          int eleGID                          ///< element GID
      );

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      virtual void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
      {
        dserror("not implemented");
      }

      /// @name Access methods
      //@{
      double Mue() const { return params_->mue_; }
      double Lambda() const { return params_->lambda_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
      )
      {
        isoprinc = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupLogNeoHooke* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // LOGNEOHOOKE_H_
