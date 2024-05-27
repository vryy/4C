/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a volumetric general power-type material in terms of the Jacobi
determinant

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUP3POW_HPP
#define FOUR_C_MATELAST_COUP3POW_HPP

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
      /*!
       * @brief material parameters for volumtric contribution of a general power material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_Coup3Pow C 1 D 1
       */
      class Coup3Pow : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup3Pow(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c_;
        int d_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<CORE::MAT::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class Coup3Pow
    }     // namespace PAR

    /*!
     * @brief Volumtric general power material
     *
     * This is a summand of a variable order hyperelastic, volumetric
     * material depending on the first invariant of the right Cauchy-Green tensor.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = C (I_3^(1/3)-1)^D.
     * \f]
     */
    class Coup3Pow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup3Pow(MAT::ELASTIC::PAR::Coup3Pow* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_coup3pow;
      }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1> glstrain,  ///< Green-Lagrange strain
          const int eleGID                            ///< element GID
      );

      void add_derivatives_principal(
          CORE::LINALG::Matrix<3, 1>& dPI,    ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<6, 1>& ddPII,  ///< second derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,     ///< Gauss point
          int eleGID  ///< element GID
          ) override;

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


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::Coup3Pow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
