/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isotropic exponential coupled material according to
Weickenmeier_2014

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPEXPPOL_HPP
#define FOUR_C_MATELAST_COUPEXPPOL_HPP

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
       * @brief material parameters for compressible soft tissue material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupExpPol A 600. B 2. C 5.
       */
      class CoupExpPol : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupExpPol(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{
        double a_;
        /// constant for linear part of I_1
        double b_;
        /// constant for linear part of J
        double c_;
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
      };  // class CoupExpPol

    }  // namespace PAR

    /*!
     * @brief coupled hyperelastic, compressible, isotropic material according to [1] with linear
     * summands in exponent The material strain energy density function is
     * \f[
     *    \Psi = a \exp[ b(I_1- 3) - (2b + c)ln{J} + c(J-1) ] - a
     * \f]
     *
     * More details at #AddCoefficientsPrincipal()
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] Weickenmeier, Jabareen "Elastic-viscoplastic modeling of soft biological
     *           tissues using a mixed finite element formulation based on the relative
     *           deformation gradient", 2014
     * </ul>
     */
    class CoupExpPol : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupExpPol(MAT::ELASTIC::PAR::CoupExpPol* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coupexppol; }

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

      // add first and second derivative w.r.t. principal invariants
      void AddDerivativesPrincipal(
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
      MAT::ELASTIC::PAR::CoupExpPol* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
