/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a coupled Mooney Rivlin material

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPMOONEYRIVLIN_HPP
#define FOUR_C_MATELAST_COUPMOONEYRIVLIN_HPP

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
       * @brief material parameters for coupled contribution of a CoupMooneyRivlinan material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1
       */
      class CoupMooneyRivlin : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupMooneyRivlin(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c1_;
        double c2_;
        double c3_;
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
      };  // class CoupMooneyRivlin

    }  // namespace PAR

    /*!
     * @brief CoupMooneyRivlinan material
     *
     * Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.
     *
     * The underlying strain-energy function is (expressed in invariants I1 and I2):
     * \f[
     *     \Psi= c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + c_3 * (J - 1)^2
     * \f]
     *
     * For references see Holzapfel p. 245
     *
     * Parameters are \f$ c_1, c_2\f$ and \f$c_3\f$ as penalty factor to enforce incompressibility
     * (shear modulus \f$\mu = (c_1 + c_2) / 2)\f$
     */
    class CoupMooneyRivlin : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupMooneyRivlin(MAT::ELASTIC::PAR::CoupMooneyRivlin* params);


      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_coupmooneyrivlin;
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

      void add_derivatives_principal(
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

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupMooneyRivlin* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
