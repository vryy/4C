/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the volumetric SussmanBathe material according to "Doll, S. and
Schweizerhof, K. On the Development of Volumetric Strain Energy Functions Journal of Applied
Mechanics, 2000"

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_VOLSUSSMANBATHE_HPP
#define FOUR_C_MATELAST_VOLSUSSMANBATHE_HPP

#include "4C_config.hpp"

#include "4C_mat_par_parameter.hpp"
#include "4C_matelast_summand.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for volumetric contribution \f$\Psi=\kappa(J-1)^2\f$
       *
       *  <h3>Input line</h3>
       *  MAT 1 ELAST_VolSussmanBathe KAPPA 10000
       */
      class VolSussmanBathe : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        VolSussmanBathe(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Dilatation modulus
        double kappa_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class VolSussmanBathe

    }  // namespace PAR

    /*!
     * @brief Volumetric SussmanBathe material according to [1].
     *
     * Strain energy function is given by
     * \f[
     *    \Psi = \frac \kappa 2 (J-1)^2
     * \f]
     *
     *  [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
     *      Journal of Applied Mechanics, 2000
     */
    class VolSussmanBathe : public Summand
    {
     public:
      /// constructor with given material parameters
      VolSussmanBathe(MAT::ELASTIC::PAR::VolSussmanBathe* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_volsussmanbathe;
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


      // Add derivatives with respect to modified invariants.
      void AddDerivativesModified(
          CORE::LINALG::Matrix<3, 1>&
              dPmodI,  ///< first derivative with respect to modified invariants
          CORE::LINALG::Matrix<6, 1>&
              ddPmodII,  ///< second derivative with respect to modified invariants
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Add third derivative w.r.t. J
      void Add3rdVolDeriv(const CORE::LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3) override;

      /// @name Access methods
      //@{
      double Kappa() const { return params_->kappa_; }
      //@}

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
      MAT::ELASTIC::PAR::VolSussmanBathe* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
