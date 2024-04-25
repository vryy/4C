/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isotropic coupled material depending on the first Cauchy-Green
invariant and the Jacobi determinant

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUP13APOW_HPP
#define FOUR_C_MATELAST_COUP13APOW_HPP

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
       * @brief material parameters for coupled contribution of a multiplicative coupled I1 and I3
       * power material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_Coup13aPow C 1 D 1 A 1
       */
      class Coup13aPow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup13aPow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c_;  // constant
        int d_;     // exponent of all
        double a_;  // negative exponent of I3

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
      };  // class Coup13aPow

    }  // namespace PAR

    /*!
     * @brief Coupled general power material
     *
     * This is a summand of a variable order hyperelastic, coupled
     * material depending on the first and third invariant of the right
     * Cauchy-Green tensor. The third invariant has a additional negative exponent.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = c (I_{\boldsymbol{C}}*(III_{\boldsymbol{C}}^(-a))-3)^d.
     * \f]
     */
    class Coup13aPow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup13aPow(MAT::ELASTIC::PAR::Coup13aPow* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coup13apow; }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy function
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1>&
              glstrain,  ///< Green-Lagrange strain in strain like Voigt notation
          int gp,        ///< Gauss point
          int eleGID     ///< element GID
          ) override;

      void AddDerivativesPrincipal(
          CORE::LINALG::Matrix<3, 1>& dPI,    ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<6, 1>& ddPII,  ///< second derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,     ///< Gauss point
          int eleGID  ///< element GID
          ) override;

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4) override
      {
        FOUR_C_THROW("not implemented");
      }

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
      MAT::ELASTIC::PAR::Coup13aPow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
