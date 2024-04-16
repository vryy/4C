/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isochoric Neo Hooke material

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_ISONEOHOOKE_HPP
#define FOUR_C_MATELAST_ISONEOHOOKE_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_matelast_summand.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a neo-Hooke material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_IsoNeoHooke MUE 100
       */
      class IsoNeoHooke : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoNeoHooke(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double mue_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class IsoNeoHooke
    }     // namespace PAR

    /*!
     * @brief Isochoric neo-Hooke material according to [1]
     *
     * This is a hyperelastic, isotropic material
     * of the most simple kind.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = \frac{\mu}{2} (\overline{I}_{\boldsymbol{C}}-3) = \frac{\mu}{2}
     *   (J^{-2/3}{I}_{\boldsymbol{C}}-3)
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * </ul>
     */
    class IsoNeoHooke : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoNeoHooke(MAT::ELASTIC::PAR::IsoNeoHooke* params);


      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_isoneohooke; }

      /// add shear modulus equivalent
      void AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
          double& shearmod                  ///< variable to add upon
      ) const override;

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

      /// @name Access methods
      //@{
      double Mue() const { return params_->mue_; }
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
      MAT::ELASTIC::PAR::IsoNeoHooke* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
