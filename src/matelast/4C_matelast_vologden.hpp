/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the volumetic contribution suggested by Ogden, see "Doll, S. and
Schweizerhof, K. On the Development of Volumetric Strain Energy Functions Journal of Applied
Mechanics, 2000"

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_VOLOGDEN_HPP
#define FOUR_C_MATELAST_VOLOGDEN_HPP

#include "4C_config.hpp"

#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace Elastic
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for volumetric contribution
       * \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_VolOgden KAPPA 100 BETA -2
       */
      class VolOgden : public Core::Mat::PAR::Parameter
      {
       public:
        /// standard constructor
        VolOgden(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Dilatation modulus
        double kappa_;
        /// empiric constant
        double beta_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<Core::Mat::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "Mat::Elastic::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class VolOgden
    }     // namespace PAR

    /*!
     * @brief Volumetric Ogden material according to [1].
     *
     *  Strain energy function is given by
     * \f[
     *    \Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)
     * \f]
     *
     * [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
     *      Journal of Applied Mechanics, 2000
     */
    class VolOgden : public Summand
    {
     public:
      /// constructor with given material parameters
      VolOgden(Mat::Elastic::PAR::VolOgden* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType MaterialType() const override
      {
        return Core::Materials::mes_vologden;
      }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy function
          const Core::LinAlg::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const Core::LinAlg::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const Core::LinAlg::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                      ///< Gauss point
          const int eleGID                             ///< element GID
          ) override;

      // Add derivatives with respect to modified invariants.
      void add_derivatives_modified(
          Core::LinAlg::Matrix<3, 1>&
              dPmodI,  ///< first derivative with respect to modified invariants
          Core::LinAlg::Matrix<6, 1>&
              ddPmodII,  ///< second derivative with respect to modified invariants
          const Core::LinAlg::Matrix<3, 1>&
              modinv,       ///< modified invariants of right Cauchy-Green tensor
          int gp,           ///< Gauss point
          const int eleGID  ///< element GID
          ) override;

      /// Add third derivative w.r.t. J
      void Add3rdVolDeriv(const Core::LinAlg::Matrix<3, 1>& modinv, double& d3PsiVolDJ3) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        isomod = true;
        return;
      };


     private:
      /// my material parameters
      Mat::Elastic::PAR::VolOgden* params_;
    };

  }  // namespace Elastic
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
