/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric contribution of a material to test the isochoric
parts of the Elasthyper-Toolbox.

\level 1
*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_ISOTESTMATERIAL_H
#define BACI_MATELAST_ISOTESTMATERIAL_H

#include "baci_mat_par_parameter.H"
#include "baci_matelast_summand.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of test material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_IsoTestMaterial C1 100 C2 50
       */
      class IsoTestMaterial : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoTestMaterial(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double c1_;
        double c2_;

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
      };  // class IsoTestMaterial
    }     // namespace PAR

    /*!
     * @brief Isochoric Material to test the elasthyper-Toolbox
     *
     * This material is not realistic, but contains all possible derivatives of invariants. With
     * this material in combination with volsussmanbathe, it is possible to test all isochoric
     * parts of the Elasthyper-Toolbox.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + 0.5 C1 (\overline{I}_{\boldsymbol{C}}-3)^2
     *        + C2 (\overline{II}_{\boldsymbol{C}}-3)  + 0.5 C2
     *        (\overline{II}_{\boldsymbol{C}}-3)^2
     *        + D (\overline{I}_{\boldsymbol{C}}-3) (\overline{II}_{\boldsymbol{C}}-3).
     * \f]
     *
     * with D = C1 + 2 C2
     */
    class IsoTestMaterial : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoTestMaterial(MAT::ELASTIC::PAR::IsoTestMaterial* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_isotestmaterial;
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
          const int eleGID                             ///< element GID
          ) override;

      void AddDerivativesModified(
          CORE::LINALG::Matrix<3, 1>&
              dPmodI,  ///< first derivative with respect to modified invariants
          CORE::LINALG::Matrix<6, 1>&
              ddPmodII,  ///< second derivative with respect to modified invariants
          const CORE::LINALG::Matrix<3, 1>&
              modinv,       ///< modified invariants of right Cauchy-Green tensor
          int gp,           ///< Gauss point
          const int eleGID  ///< element GID
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
      MAT::ELASTIC::PAR::IsoTestMaterial* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
