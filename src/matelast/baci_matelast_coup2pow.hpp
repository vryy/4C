/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isotropic general power-type material in terms of the second
Cauchy-Green invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUP2POW_HPP
#define FOUR_C_MATELAST_COUP2POW_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_matelast_summand.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a general power material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_Coup2Pow C 1 D 1
       */
      class Coup2Pow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup2Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c_;
        int d_;

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
      };  // class Coup1Pow

    }  // namespace PAR

    /*!
     * @brief Isochoric general power material
     *
     * This is a summand of variable order hyperelastic, isotropic
     * material depending on the seccond invariant of the right Cauchy-Green tensor.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = C (II_{\boldsymbol{C}}-3)^D.
     * \f]
     */
    class Coup2Pow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup2Pow(MAT::ELASTIC::PAR::Coup2Pow* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coup2pow; }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy function
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                      ///< Gauss point
          const int eleGID                             ///< element GID
          ) override;


      void AddDerivativesPrincipal(
          CORE::LINALG::Matrix<3, 1>& dPI,    ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<6, 1>& ddPII,  ///< second derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>&
              prinv,        ///< principal invariants of right Cauchy-Green tensor
          int gp,           ///< Gauss point
          const int eleGID  ///< element GID
          ) override;

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4) override
      {
        dserror("not implemented");
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
      MAT::ELASTIC::PAR::Coup2Pow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
