/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a coupled Saint-Venant-Kirchhoff material

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPSAINTVENANTKIRCHHOFF_HPP
#define FOUR_C_MATELAST_COUPSAINTVENANTKIRCHHOFF_HPP

#include "4C_config.hpp"

#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      class CoupSVK : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupSVK(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double mue_;
        /// Lame's constant
        double lambda_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };
    }  // namespace PAR

    /*!
     * @brief Saint Venant Kirchhoff - Material
     *
     *
     * The material strain energy density function is
     * \f[
     * \Psi = \Big(\frac{1}{4}\mu + \frac{1}{8}\lambda \Big) \, I_{1}^2
     *       - \big( \frac{3}{4}\lambda + \frac{1}{2}\mu \big)  \, I_{1} - \frac{1}{2}\mu \, I_{2}
     *      + \frac{9}{8}\lambda + \frac{3}{4}\mu
     * \f]
     * More details at #AddCoefficientsPrincipal()
     */
    class CoupSVK : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupSVK(MAT::ELASTIC::PAR::CoupSVK* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_coupSVK;
      }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy function
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                      ///< Gauss point
          int eleGID                                   ///< element GID
          ) override;

      void AddDerivativesPrincipal(
          CORE::LINALG::Matrix<3, 1>& dPI,    ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<6, 1>& ddPII,  ///< second derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,     ///< Gauss point
          int eleGID  ///< element GID
          ) override;

      void AddThirdDerivativesPrincipalIso(
          CORE::LINALG::Matrix<10, 1>&
              dddPIII_iso,  ///< third derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>& prinv_iso,  ///< principal isotropic invariants
          int gp,                                       ///< Gauss point
          int eleGID) override;                         ///< element GID

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        isoprinc = true;
        return;
      };


      /// a young's modulus equivalent
      void AddYoungsMod(double& young, double& shear, double& bulk) override
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

FOUR_C_NAMESPACE_CLOSE

#endif
