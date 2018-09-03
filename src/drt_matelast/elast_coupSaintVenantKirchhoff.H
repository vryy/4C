/*----------------------------------------------------------------------*/
/*!
\file elast_coupSaintVenantKirchhoff.H

\brief
This file contains the routines required for Saint-Venant-Kirchhoff material
The input line should read
  MAT 1 ELAST_CoupSVK YOUNG 1.044E7 NUE 0.3

\level 1

<pre>
\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_SVK_H_
#define ELAST_SVK_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      class CoupSVK : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupSVK(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~CoupSVK() { ; }

        /// @name material parameters
        //@{

        /// Shear modulus
        double mue_;
        /// Lame's constant
        double lambda_;

        //@}

        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        void SetLambda(const double lambda) { lambda_ = lambda; }
        void SetMue(const double mue) { mue_ = mue; }
        //@}
      };

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Saint Venant Kirchhoff - Material
    ///
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = \Big(\frac{1}{4}\mu + \frac{1}{8}\lambda \Big) \, I_{1}^2
    ///       - \big( \frac{3}{4}\lambda + \frac{1}{2}\mu \big)  \, I_{1} - \frac{1}{2}\mu \, I_{2}
    ///      + \frac{9}{8}\lambda + \frac{3}{4}\mu
    ///\f]
    /// More details at #AddCoefficientsPrincipal()
    ///
    ///
    /// \date 12/14


    class CoupSVK : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupSVK(MAT::ELASTIC::PAR::CoupSVK* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupSVK; }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          const int eleGID                       ///< element GID
      );

      virtual void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const int eleGID                    ///< element GID
      );

      virtual void AddThirdDerivativesPrincipalIso(
          LINALG::Matrix<10, 1>& dddPIII_iso,     ///< third derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv_iso,  ///< principal isotropic invariants
          const int eleGID);                      ///< element GID

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      virtual void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4);

      /// @name Access methods
      //@{
      double LAMBDA() const { return params_->lambda_; }
      double MUE() const { return params_->mue_; }
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


      /// a young's modulus equivalent
      virtual void AddYoungsMod(double& young, double& shear, double& bulk)
      {
        young += 9. * params_->mue_ * (3. * params_->lambda_ + 2. * params_->mue_) /
                 (params_->lambda_ + params_->mue_);
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupSVK* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ELAST_SVK_H_
