/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required for the volumetic function
suggested by Ogden, see Doll, S. and Schweizerhof, K.
On the Development of Volumetric Strain Energy Functions
Journal of Applied Mechanics, 2000
The input line should read
  MAT 1 ELAST_VolOgden KAPPA 100 BETA -2

\level 1

\maintainer Fabian Braeu

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLOGDEN_H_
#define ELAST_VOLOGDEN_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for volumetric contribution
      /// \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_VolOgden KAPPA 100 BETA -2
      class VolOgden : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        VolOgden(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~VolOgden() { ; }

        /// @name material parameters
        //@{

        /// Dilatation modulus
        double kappa_;
        /// empiric constant
        double beta_;

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
        void SetKappa(const double kappa) { kappa_ = kappa; }
        void SetBeta(const double beta) { beta_ = beta; }
        //@}


      };  // class VolOgden

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Volumetric Ogden material according to [1].
    ///
    /// Strain energy function is given by
    ///
    ///
    ///\f[
    ///   \Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)
    ///\f]
    ///
    ///
    /// [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
    ///     Journal of Applied Mechanics, 2000
    /// \author rausch
    /// \date 04/09
    class VolOgden : public Summand
    {
     public:
      /// constructor with given material parameters
      VolOgden(MAT::ELASTIC::PAR::VolOgden* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_vologden; }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          const int eleGID                       ///< element GID
      );

      // Add derivatives with respect to modified invariants.
      virtual void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const int eleGID                     ///< element GID
      );

      /// Add third derivative w.r.t. J
      virtual void Add3rdVolDeriv(const LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3);

      /// @name Access methods
      //@{
      double Kappa() const { return params_->kappa_; }
      double Beta() const { return params_->beta_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
      )
      {
        isomod = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::VolOgden* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // VOLOGDEN_H_
