/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric contribution of a Mooney-Rivlin-type material

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_ISOMOONEYRIVLIN_HPP
#define FOUR_C_MATELAST_ISOMOONEYRIVLIN_HPP

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
       * @brief material parameters for isochoric contribution of a Mooney-Rivlin type material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_IsoMooneyRivlin C1 100 C2 50
       */
      class IsoMooneyRivlin : public Core::Mat::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoMooneyRivlin(const Core::Mat::PAR::Parameter::Data& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double c1_;
        double c2_;

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
      };  // class IsoMooneyRivlin
    }     // namespace PAR

    /*!
     * @brief Isochoric Mooney-Rivlin material according to [1]
     *
     * Strain energy function is given by
     * \f[
     *    \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{II}_{\boldsymbol{C}}-3).
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * </ul>
     */
    class IsoMooneyRivlin : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoMooneyRivlin(Mat::Elastic::PAR::IsoMooneyRivlin* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType MaterialType() const override
      {
        return Core::Materials::mes_isomooneyrivlin;
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
          int eleGID                                   ///< element GID
          ) override;

      // Add derivatives with respect to modified invariants.
      void add_derivatives_modified(
          Core::LinAlg::Matrix<3, 1>&
              dPmodI,  ///< first derivative with respect to modified invariants
          Core::LinAlg::Matrix<6, 1>&
              ddPmodII,  ///< second derivative with respect to modified invariants
          const Core::LinAlg::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

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
      Mat::Elastic::PAR::IsoMooneyRivlin* params_;
    };

  }  // namespace Elastic
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
